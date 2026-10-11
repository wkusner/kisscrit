#!/bin/bash
# =====================================================================================================
#  kisscrit completeness run: one-command setup and run for a Mac (Intel or Apple Silicon), from Terminal
#
#  Usage (one line, from Terminal):
#      bash ~/Downloads/mac_run.sh
#
#  What it does (every step is skipped if already done, so the same line also RESUMES a stopped run):
#    1. installs Miniforge (a private Python, in ~/miniforge3; no admin password, no Homebrew)
#    2. creates a Python environment "kisscrit" with numpy, scipy, networkx, mpmath, numba, git, gh
#    3. logs you in to GitHub once (opens your browser; you type an 8-character code)
#    4. gets the code: clones github.com/<you>/kisscrit, or, if that repository does not exist yet and
#       ~/Downloads/kisscrit-update.bundle is present, creates it as a PRIVATE repo from the bundle
#    5. runs the self-test (1-3 minutes)
#    6. phase A: generates candidate graphs (4 triangulations in parallel), commits and pushes them
#    7. phase B: claims shards one at a time, certifies them on all cores, commits and pushes each one
#    8. prints the summary and pushes it
#
#  After the self-test the run is handed to launchd (a LaunchAgent, org.kisscrit.run): it keeps running
#  with the Terminal window closed, restarts by itself after a reboot or log-in, and keeps the Mac awake
#  (caffeinate; on a laptop keep it plugged in).  Progress:  tail -f ~/kisscrit-run/logs/run.log
#  Stop for good:   launchctl bootout gui/$(id -u)/org.kisscrit.run
#  Interactive run without launchd:   AUTOSTART=0 bash mac_run.sh
#
#  Settings you can override, e.g.   REPO=critconfig WORKERS=3 bash ~/Downloads/mac_run.sh
#    REPO     repository name on your GitHub account                (default: kisscrit)
#    ROUTE    cond (icosahedral-type hulls only) or full             (default: cond)
#    SHARDS   number of certify shards, fixed for the whole run      (default: 40)
#    WORKERS  parallel processes                                     (default: number of cores)
#    BUDGET   max boxes per candidate                                (default: 300000)
#    BUNDLE   path of the git bundle used to seed a new repository   (default: ~/Downloads/kisscrit-update.bundle)
# =====================================================================================================
set -e
REPO="${REPO:-kisscrit}"
ROUTE="${ROUTE:-cond}"
SHARDS="${SHARDS:-40}"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu)}"
BUDGET="${BUDGET:-300000}"
BUNDLE="${BUNDLE:-$HOME/Downloads/kisscrit-update.bundle}"
WORK="$HOME/kisscrit-run"
MF="$HOME/miniforge3"
NAME="$(scutil --get ComputerName 2>/dev/null | tr ' ' '-' || hostname -s)"

AUTOSTART="${AUTOSTART:-1}"
DETACHED=0; [ "$1" = "--detached" ] && DETACHED=1
LABEL="org.kisscrit.run"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
mkdir -p "$WORK/logs"

# one copy at a time
LOCK="$WORK/.lock"
if mkdir "$LOCK" 2>/dev/null; then echo $$ > "$LOCK/pid"
elif kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; then
  echo "A run is already going (pid $(cat "$LOCK/pid")).  Progress: tail -f $WORK/logs/run.log"; exit 0
else echo $$ > "$LOCK/pid"; fi
trap 'rm -rf "$LOCK"' EXIT

# keep the Mac awake for as long as this script runs (-s: also on a laptop, while on power)
caffeinate -ims -w $$ &

if [ -z "$SKIP_SETUP" ]; then   # (SKIP_SETUP=1 is for testing only)
# ---------------------------------------------------------------- 1. Miniforge
if [ ! -x "$MF/bin/conda" ]; then
  say "Installing Miniforge (private Python) into $MF"
  ARCH="$(uname -m)"   # x86_64 on an Intel iMac, arm64 on Apple Silicon
  curl -fsSL -o "$WORK/miniforge.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-${ARCH}.sh"
  bash "$WORK/miniforge.sh" -b -p "$MF"
fi
# shellcheck disable=SC1091
source "$MF/etc/profile.d/conda.sh"

# ---------------------------------------------------------------- 2. environment
if ! conda env list | grep -q '^kisscrit '; then
  say "Creating the Python environment (a few minutes)"
  conda create -y -q -n kisscrit -c conda-forge python=3.11 numpy scipy networkx mpmath numba git gh
fi
conda activate kisscrit

# ---------------------------------------------------------------- 3. GitHub login
if ! gh auth status >/dev/null 2>&1; then
  [ "$DETACHED" = 1 ] && { echo "Not logged in to GitHub: run  bash $WORK/mac_run.sh  in Terminal once."; exit 1; }
  say "Logging in to GitHub: a code is shown below; your browser opens; paste the code there"
  gh auth login --hostname github.com --git-protocol https --web
fi
gh auth setup-git
USER_LOGIN="$(gh api user --jq .login)"
USER_ID="$(gh api user --jq .id)"
git config --global user.name  >/dev/null || git config --global user.name "$USER_LOGIN"
git config --global user.email >/dev/null || git config --global user.email "${USER_ID}+${USER_LOGIN}@users.noreply.github.com"
git config --global pull.rebase true

# ---------------------------------------------------------------- 4. the code
SRC="$WORK/$REPO"
if [ ! -d "$SRC/.git" ]; then
  if gh repo view "$USER_LOGIN/$REPO" >/dev/null 2>&1; then
    say "Cloning $USER_LOGIN/$REPO"
    gh repo clone "$USER_LOGIN/$REPO" "$SRC"
  elif [ ! -f "$BUNDLE" ] && grep -q '^__BUNDLE_BASE64__$' "$0"; then
    # self-contained copy of this script: unpack the embedded bundle
    BUNDLE="$WORK/kisscrit-update.bundle"
    sed -n '/^__BUNDLE_BASE64__$/,$p' "$0" | tail -n +2 \
      | python -c "import sys, base64; sys.stdout.buffer.write(base64.b64decode(sys.stdin.read()))" > "$BUNDLE"
  fi
  if gh repo view "$USER_LOGIN/$REPO" >/dev/null 2>&1; then :
  elif [ -f "$BUNDLE" ]; then
    say "Creating private repository $USER_LOGIN/$REPO from $BUNDLE"
    git clone -q "$BUNDLE" "$SRC"
    (cd "$SRC" && git checkout -q -B main && git remote remove origin \
       && gh repo create "$USER_LOGIN/$REPO" --private --source=. --remote=origin --push)
  else
    echo "No repository $USER_LOGIN/$REPO on GitHub and no bundle at $BUNDLE."
    echo "Download kisscrit-update.bundle into ~/Downloads (or set REPO=<name> of a repo that has the code) and rerun."
    exit 1
  fi
fi
fi
SRC="${SRC_OVERRIDE:-$WORK/$REPO}"
cd "$SRC"
git pull -q || true
# the code may live in a subfolder of the repository (e.g. critconfig/kisscrit)
RS="$(find . -path ./.git -prune -o -name run_shard.py -path '*scripts/completeness/*' -print | head -1)"
if [ -z "$RS" ]; then echo "This repository has no scripts/completeness/run_shard.py: push the new code first."; exit 1; fi
ROOT="$(cd "$(dirname "$RS")/../.." && pwd)"
cd "$ROOT"
echo "code at $ROOT   machine name: $NAME   workers: $WORKERS"

next_shard() {  # prints a newly claimed shard; returns 1 when none are left, 2 on a git error
  local out rc
  if out=$(python scripts/completeness/claim.py --phase "$1" --route "$ROUTE" --of "$SHARDS" --name "$NAME"); then
    echo "$out"; return 0
  else
    rc=$?
    [ $rc -eq 1 ] && return 1
    return 2
  fi
}

push() {  # commit result files and push, retrying if another machine pushed first
  git add -A results/completeness >/dev/null 2>&1 || true
  git commit -q -m "$1" >/dev/null 2>&1 || true
  for i in 1 2 3 4 5; do
    git pull -q --rebase && git push -q && return 0
    sleep $((RANDOM % 20 + 5))
  done
  echo "push failed; results are committed locally and will be pushed next time"
}

# ---------------------------------------------------------------- 5. self-test
say "Self-test"
python scripts/completeness/selftest.py | tee "$WORK/logs/selftest.log" | tail -3

# ---------------------------------------------------------------- hand over to launchd
if [ "$DETACHED" = 0 ] && [ "$AUTOSTART" = 1 ]; then
  cp "$0" "$WORK/mac_run.sh" 2>/dev/null || true
  mkdir -p "$HOME/Library/LaunchAgents"
  cat > "$PLIST" <<PL
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>Label</key><string>$LABEL</string>
  <key>ProgramArguments</key><array><string>/bin/bash</string><string>$WORK/mac_run.sh</string><string>--detached</string></array>
  <key>EnvironmentVariables</key><dict>
    <key>REPO</key><string>$REPO</string><key>ROUTE</key><string>$ROUTE</string>
    <key>SHARDS</key><string>$SHARDS</string><key>WORKERS</key><string>$WORKERS</string>
    <key>BUDGET</key><string>$BUDGET</string><key>BUNDLE</key><string>$BUNDLE</string>
  </dict>
  <key>RunAtLoad</key><true/>
  <key>ProcessType</key><string>Background</string>
  <key>Nice</key><integer>5</integer>
  <key>StandardOutPath</key><string>$WORK/logs/run.log</string>
  <key>StandardErrorPath</key><string>$WORK/logs/run.log</string>
</dict></plist>
PL
  launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
  rm -rf "$LOCK"; trap - EXIT
  launchctl bootstrap "gui/$(id -u)" "$PLIST"
  say "The run continues in the background (launchd job $LABEL). You can close this window."
  echo "Progress:        tail -f $WORK/logs/run.log"
  echo "Summary so far:  cd $ROOT && git pull -q && $MF/envs/kisscrit/bin/python scripts/completeness/merge.py --route $ROUTE"
  echo "Stop for good:   launchctl bootout gui/\$(id -u)/$LABEL"
  exit 0
fi

# ---------------------------------------------------------------- 6. phase A: candidates
say "Phase A: candidate generation ($ROUTE route)"
if [ "$ROUTE" = "cond" ]; then
  for k in 0 1 2 3; do
    python scripts/completeness/run_shard.py gen --route cond --shard $k --of 4 > "$WORK/logs/gen-$k.log" 2>&1 &
  done
  wait
  tail -n 2 "$WORK"/logs/gen-*.log
  push "gen $ROUTE ($NAME)"
else
  while true; do
    K=$(next_shard gen) || { rc=$?; [ $rc -eq 1 ] && break; echo "Could not claim a shard (git error above). Check the network and rerun the same command."; exit 1; }
    python scripts/completeness/run_shard.py gen --route "$ROUTE" --shard "$K" --of "$SHARDS" | tee "$WORK/logs/gen-$K.log"
    push "gen $ROUTE shard $K/$SHARDS ($NAME)"
  done
fi

# ---------------------------------------------------------------- 7. phase B: certification
say "Phase B: certification, $SHARDS shards, claimed one at a time"
# first finish any shard this machine claimed earlier but did not complete (interrupted run)
for f in results/completeness/claims/certify-"$ROUTE"-*-of-"$SHARDS".txt; do
  [ -f "$f" ] || continue
  case "$(cat "$f")" in "$NAME "*) ;; *) continue ;; esac
  K="$(basename "$f" .txt | sed -E 's/^certify-[a-z]+-([0-9]+)-of-.*/\1/')"
  [ -f "results/completeness/cert/done-$ROUTE-$K-of-$SHARDS.txt" ] && continue
  echo "--- resuming shard $K of $SHARDS ($(date))"
  python scripts/completeness/run_shard.py certify --route "$ROUTE" --shard "$K" --of "$SHARDS" \
         --workers "$WORKERS" --budget "$BUDGET" 2>&1 | tee -a "$WORK/logs/cert-$K.log"
  push "certify $ROUTE shard $K/$SHARDS ($NAME)"
done
while true; do
  K=$(next_shard certify) || { rc=$?; [ $rc -eq 1 ] && break; echo "Could not claim a shard (git error above). Check the network and rerun the same command."; exit 1; }
  echo "--- shard $K of $SHARDS ($(date))"
  python scripts/completeness/run_shard.py certify --route "$ROUTE" --shard "$K" --of "$SHARDS" \
         --workers "$WORKERS" --budget "$BUDGET" 2>&1 | tee "$WORK/logs/cert-$K.log"
  push "certify $ROUTE shard $K/$SHARDS ($NAME)"
done

# ---------------------------------------------------------------- 8. summary
say "Summary"
python scripts/completeness/merge.py --route "$ROUTE" | tee "$WORK/logs/summary.log"
git add -A results/completeness >/dev/null 2>&1 || true
git commit -q -m "summary $ROUTE ($NAME)" >/dev/null 2>&1 || true
git pull -q --rebase && git push -q || true
say "Done. All shards are claimed. If other machines are still working, rerun later for a fresh summary."
[ "$DETACHED" = 1 ] && launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
