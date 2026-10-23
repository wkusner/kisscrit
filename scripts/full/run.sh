#!/usr/bin/env bash
# run.sh: the kisscrit completeness run on a Linux server without root (checklist F8; replaces mac_run.sh for
# the new run).  Status: untrusted (run tooling, no part of the proof).  The search writes certificates; the
# checker named in the run record decides; results are labelled "checker verdict at commit X; not certified".
#
# Version 1 (2026-10-22, F8 milestone M1).  Linux first; macOS (launchd, Keychain) is milestone M4 and refused
# here.  GitHub pushes (M2) and Zenodo (M3) are not in this version: results are committed locally on the run
# branch run/<route>-<run-id> and never pushed.
#
# Subcommands (each also takes --dry-run: print every command that changes something, run nothing):
#   setup      [--commit REF] [--repo OWNER/NAME | --repo-url URL] [--archive DIR] [--skip-preflight]
#              Miniforge in ~/miniforge3, the conda env "kisscrit-run", gh login (device code), the code in
#              ~/kisscrit-run/repo (gh clone) or from a `git archive` tree with a COMMIT file; with --commit,
#              then the self-test.  No root, nothing outside $HOME.
#   preflight  [--need-disk-gb N]    runs scripts/full/kisscrit-preflight.sh (unchanged); stops on any FAIL
#   selftest   --commit REF [--pytest-subset "FILES"]   C build and make test, pytest, selftest.py, peak RSS
#   start      --route test|cond|full --commit REF [options]   or just `start` to resume the current run
#              options: --shards K (40) --budget B (300000) --workers W (nproc) --workers-per-shard w (2)
#              --only-shards LIST --fresh --checker SCRIPT --survival auto|systemd|tmux|foreground
#              --run-id ID --new --need-disk-gb N (50) --prepare-only --allow-partial-selftest
#   status     [--final]
#   stop       [--reason TEXT]
#   check-tokens               whether the Zenodo token files exist and their mode (never their contents)
#   zenodo-sandbox-test        milestone M3, not in this version
# REF is a tag or a full 40-character commit hash; branch names and short hashes are refused.
# Everything lives under ~/kisscrit-run (scratch, logs, trees); never /tmp.  Help: bash run.sh help

set -uo pipefail
export LC_ALL=C

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
HERE="$(dirname "$SELF")"
RUN_HOME="${KISSCRIT_RUN_HOME:-$HOME/kisscrit-run}"
MF="${KISSCRIT_MINIFORGE:-$HOME/miniforge3}"
ENV_NAME=kisscrit-run
ENV_DIR="$MF/envs/$ENV_NAME"
TOKDIR="$HOME/.config/kisscrit"
UNIT=kisscrit-run.service
UNIT_FILE="$HOME/.config/systemd/user/$UNIT"
TMUX_SESSION=kisscrit-run
TMUX_SOCK="${KISSCRIT_TMUX_SOCKET:-kisscrit}"
CRON_TAG="# kisscrit-run @reboot line (scripts/full/run.sh)"
CRONTAB="${KISSCRIT_CRONTAB:-crontab}"
SYSTEMCTL="${KISSCRIT_SYSTEMCTL:-systemctl}"
LOGINCTL="${KISSCRIT_LOGINCTL:-loginctl}"
DEFAULT_REPO=wkusner/kisscrit
DEFAULT_CHECKER=scripts/completeness/check_certificates.py
CONDA_PKGS=(python=3.12 numpy scipy numba python-flint libflint c-compiler make networkx mpmath sympy pytest git gh)
DRY=0

# ------------------------------------------------------------------------------------------------- helpers
say()  { printf '%s\n' "$*"; }
note() { printf '  %s\n' "$*"; }
die()  { printf 'run.sh: %s\n' "$*" >&2; exit 1; }
utc()  { date -u '+%Y-%m-%dT%H:%M:%SZ'; }
# have CMD: CMD is on PATH (KISSCRIT_HIDE="curl ..." pretends some are missing; used by the tests only)
have() { case " ${KISSCRIT_HIDE:-} " in *" $1 "*) return 1 ;; esac; command -v "$1" >/dev/null 2>&1; }
# run CMD...: run it, or print it under --dry-run.  Every command that changes state goes through run.
run()  { if [ "$DRY" = 1 ]; then printf '[dry-run]'; printf ' %q' "$@"; printf '\n'; else "$@"; fi; }
# write_file PATH: stdin to PATH (printed instead under --dry-run)
write_file() {
  if [ "$DRY" = 1 ]; then printf '[dry-run] write %q:\n' "$1"; sed 's/^/[dry-run]   /'
  else mkdir -p "$(dirname "$1")" && cat > "$1.tmp.$$" && mv -f "$1.tmp.$$" "$1"; fi
}
os_name() { uname -s; }
nproc_() { nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1; }
mem_avail_kib() { awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo 0; }
free_gb() { local p=$1; while [ ! -e "$p" ] && [ "$p" != / ]; do p=$(dirname "$p"); done; df -Pk "$p" | awk 'NR==2{printf "%d", $4/1048576}'; }
# the Python of the run: KISSCRIT_PYTHON (tests: the system Python), else the conda env's
run_python() {
  if [ -n "${KISSCRIT_PYTHON:-}" ]; then printf '%s' "$KISSCRIT_PYTHON"
  elif [ -x "$ENV_DIR/bin/python" ]; then printf '%s' "$ENV_DIR/bin/python"
  else return 1; fi
}
HPY="$(run_python 2>/dev/null || command -v python3 || true)"
# the helper: next to this file, or (for a lone copy such as ~/run.sh) in the current run's tree, the clone, or a
# self-test tree; only reading and writing the run's JSON files depends on it before the tree's own copy takes over
HELPER="$HERE/run_state.py"
if [ ! -f "$HELPER" ]; then
  for _h in "$RUN_HOME/runs/$(cat "$RUN_HOME/current" 2>/dev/null)/tree/scripts/full/run_state.py" \
            "$RUN_HOME/repo/scripts/full/run_state.py" "$RUN_HOME"/selftest/*/tree/scripts/full/run_state.py; do
    [ -f "$_h" ] && { HELPER=$_h; break; }
  done
fi
st() { [ -f "$HELPER" ] || die "missing $HELPER (run.sh needs its helper next to it; use the copy inside a run tree)"
       "$HPY" -I "$HELPER" "$@"; }
# st_w: a helper call that changes state (printed under --dry-run)
st_w() { if [ "$DRY" = 1 ]; then run "$HPY" -I "$HELPER" "$@"; else st "$@"; fi; }
[ -d "$ENV_DIR/bin" ] && [ -z "${KISSCRIT_PYTHON:-}" ] && PATH="$ENV_DIR/bin:$PATH"

# ------------------------------------------------------------------------------------------------- sources
# The code comes from a gh clone in $RUN_HOME/repo, or from a tree sent with git archive plus a COMMIT file.
# $RUN_HOME/source.env records which (written by setup).
SRC_KIND=""; SRC_DIR=""
load_source() {
  [ -f "$RUN_HOME/source.env" ] || die "no code source yet: run  bash run.sh setup  first (or setup --archive DIR)"
  # shellcheck disable=SC1091
  . "$RUN_HOME/source.env"; SRC_KIND=$KC_SRC_KIND; SRC_DIR=$KC_SRC_DIR
}
# resolve_commit REF: the full hash of a tag or a full hash; anything else is refused
resolve_commit() {
  local ref=$1 h
  if [ "$SRC_KIND" = archive ]; then
    h=$(tr -d ' \n' < "$SRC_DIR/COMMIT")
    [[ $h =~ ^[0-9a-f]{40}$ ]] || die "$SRC_DIR/COMMIT does not hold a full 40-character commit hash"
    [ "$ref" = "$h" ] || die "--commit $ref does not match the archive's COMMIT file ($h)"
    printf '%s' "$h"; return
  fi
  if [ "$DRY" = 1 ]; then printf '[dry-run] git -C %q fetch -q --tags origin\n' "$SRC_DIR" >&2
  else git -C "$SRC_DIR" fetch -q --tags origin 2>/dev/null || note "(could not fetch from origin; using the clone as it is)" >&2; fi
  if [[ $ref =~ ^[0-9a-f]{40}$ ]]; then
    git -C "$SRC_DIR" cat-file -e "$ref^{commit}" 2>/dev/null || die "commit $ref is not in $SRC_DIR (pushed? fetched?)"
    printf '%s' "$ref"
  elif h=$(git -C "$SRC_DIR" rev-parse -q --verify "refs/tags/$ref^{commit}" 2>/dev/null); then
    printf '%s' "$h"
  else
    die "--commit must be a tag or a full 40-character commit hash, not '$ref' (branch names and short hashes are refused)"
  fi
}
# code_dirty DIR: the git work tree has changes outside results/ (code must be exactly the pinned commit)
code_dirty() { [ -n "$(git -C "$1" status --porcelain --untracked-files=normal -- . ':(exclude)results' 2>/dev/null)" ]; }
git_id() {  # git -c options giving an identity even when none is configured
  local n e; n=$(git config --global user.name 2>/dev/null || true); e=$(git config --global user.email 2>/dev/null || true)
  printf '%s\n' -c "user.name=${n:-kisscrit-run}" -c "user.email=${e:-kisscrit-run@$(hostname -s 2>/dev/null || echo localhost)}"
}
gitc() { local -a id; mapfile -t id < <(git_id); run git -C "$TREE" "${id[@]}" "$@"; }

# ------------------------------------------------------------------------------------------------- current run
CUR_FILE="$RUN_HOME/current"
RUN_ID=""; RDIR=""; CFG=""; STATE=""; TREE=""; KC_TREE=""
set_run() { RUN_ID=$1; RDIR="$RUN_HOME/runs/$RUN_ID"; CFG="$RDIR/config.json"; STATE="$RDIR/state/state.json"; }
load_current() {
  [ -f "$CUR_FILE" ] || return 1
  set_run "$(cat "$CUR_FILE")"; [ -f "$CFG" ] || return 1
  eval "$(st env "$CFG")"; TREE=$KC_TREE; return 0
}
REC() { printf '%s' "$TREE/results/full/run-record.json"; }
lock_free() { ( exec 9>>"$RUN_HOME/run.lock"; flock -n 9 ) 2>/dev/null; }
stopped_flag() { printf '%s' "$RDIR/state/stopped"; }

# ------------------------------------------------------------------------------------------------- preflight
cmd_preflight() {
  local need=50 pf rep fails
  while [ $# -gt 0 ]; do case $1 in --need-disk-gb) need=$2; shift ;; *) die "preflight: unknown option $1" ;; esac; shift; done
  pf="${KISSCRIT_PREFLIGHT:-}"
  for c in "$HERE/kisscrit-preflight.sh" "$SRC_DIR/scripts/full/kisscrit-preflight.sh" "$HOME/kisscrit-preflight.sh"; do
    [ -z "$pf" ] && [ -f "$c" ] && pf=$c
  done
  [ -n "$pf" ] || die "kisscrit-preflight.sh not found next to run.sh or in ~; copy scripts/full/kisscrit-preflight.sh to the server (docs/setup-guide.md A4)"
  say "== preflight: bash $pf --need-disk-gb $need --results-dir $RUN_HOME"
  run bash "$pf" --need-disk-gb "$need" --results-dir "$RUN_HOME"
  [ "$DRY" = 1 ] && return 0
  rep="$HOME/kisscrit-preflight-report.txt"
  [ -f "$rep" ] || die "the preflight wrote no report ($rep)"
  fails=$(grep -F '[FAIL]' "$rep" || true)
  if [ -n "$fails" ]; then
    say ""; say "Preflight found problems; the run cannot start until they are fixed:"
    printf '%s\n' "$fails" | sed 's/^ *\[FAIL\] /  - /'
    printf '%s\n' "$fails" | grep -q 'no response' && say "  A host that does not answer cannot be reached from this machine: ask IT to allow outbound HTTPS to it (docs/setup-guide.md A4)."
    printf '%s\n' "$fails" | grep -q 'free where the results go' && say "  The disk is too small for the run: free space under $RUN_HOME, or ask for space (docs/setup-guide.md A4)."
    exit 3
  fi
  say "preflight: no failures ($(grep -c '\[warn\]' "$rep") warnings; report $rep)"
  say "open files: soft $(ulimit -Sn), hard $(ulimit -Hn); the run raises the soft limit to $(nofile_target)"
}
nofile_target() { local h; h=$(ulimit -Hn); if [ "$h" = unlimited ] || [ "$h" -ge 65536 ]; then echo 65536; else echo "$h"; fi; }

# ------------------------------------------------------------------------------------------------- tokens
# token_mode FILE: octal permission bits; token_ok FILE: owned by me, no group or other bits
token_mode() { stat -c %a "$1" 2>/dev/null || stat -f %Lp "$1"; }
token_ok() { local m; m=$(token_mode "$1"); [ -O "$1" ] && [ $((8#$m & 8#077)) -eq 0 ]; }
cmd_check_tokens() {
  local bad=0 f m
  say "Zenodo tokens in $TOKDIR (contents are never read here; only existence, owner and mode):"
  if [ -d "$TOKDIR" ]; then note "directory: mode $(token_mode "$TOKDIR")"; else note "directory: missing (docs/setup-guide.md A5)"; fi
  for f in zenodo.token zenodo-sandbox.token; do
    if [ ! -e "$TOKDIR/$f" ]; then note "$f: missing"
    else
      m=$(token_mode "$TOKDIR/$f")
      if token_ok "$TOKDIR/$f"; then note "$f: present, mode $m (ok)"
      else note "$f: present, mode $m: REFUSED (readable by others or not yours; run chmod 600 $TOKDIR/$f)"; bad=1; fi
    fi
  done
  return $bad
}
cmd_zenodo_sandbox_test() { die "zenodo-sandbox-test is milestone M3 of F8 and is not in this version of run.sh"; }

# ------------------------------------------------------------------------------------------------- setup
download() {  # download URL OUT: curl, or wget when curl is missing (the server has only wget)
  if have curl; then run curl -fsSL --retry 3 -o "$2" "$1"
  elif have wget; then run wget -q --tries=3 -O "$2" "$1"
  else die "neither curl nor wget is installed"; fi
}
cmd_setup() {
  local commit="" repo=$DEFAULT_REPO archive="" skip_pf=0 url="${KISSCRIT_REPO_URL:-}"
  while [ $# -gt 0 ]; do
    case $1 in
      --commit) commit=$2; shift ;; --repo) repo=$2; shift ;; --archive) archive=$2; shift ;;
      --repo-url) url=$2; shift ;;   # clone with plain git from a URL instead of gh (tests: a local bare repo)
      --skip-preflight) skip_pf=1 ;; *) die "setup: unknown option $1" ;;
    esac; shift
  done
  [ "$(os_name)" = Linux ] || die "this version supports Linux only; macOS (launchd, Keychain) is milestone M4"
  run mkdir -p "$RUN_HOME/logs" "$RUN_HOME/tmp" "$RUN_HOME/dl"
  [ $skip_pf = 1 ] || cmd_preflight
  # 1. Miniforge and the environment (skipped when KISSCRIT_PYTHON names a Python: tests only)
  if [ -n "${KISSCRIT_PYTHON:-}" ]; then say "== KISSCRIT_PYTHON=$KISSCRIT_PYTHON: no Miniforge, no conda env (test mode)"
  else
    if [ ! -x "$MF/bin/conda" ]; then
      say "== installing Miniforge into $MF (no root)"
      download "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" "$RUN_HOME/dl/miniforge.sh"
      run bash "$RUN_HOME/dl/miniforge.sh" -b -p "$MF"
    fi
    if [ ! -x "$ENV_DIR/bin/python" ]; then
      say "== creating the conda environment $ENV_NAME (a few minutes; log $RUN_HOME/logs/conda.log)"
      if [ "$DRY" = 1 ]; then run "$MF/bin/conda" create -y -q -n "$ENV_NAME" -c conda-forge --override-channels "${CONDA_PKGS[@]}"
      else "$MF/bin/conda" create -y -q -n "$ENV_NAME" -c conda-forge --override-channels "${CONDA_PKGS[@]}" > "$RUN_HOME/logs/conda.log" 2>&1 \
             || die "conda environment failed (see $RUN_HOME/logs/conda.log)"; fi
    fi
    PATH="$ENV_DIR/bin:$PATH"; HPY="$ENV_DIR/bin/python"
  fi
  # 2. the code
  if [ -n "$archive" ]; then
    archive=$(cd "$archive" 2>/dev/null && pwd) || die "--archive: no directory $archive"
    [ -f "$archive/COMMIT" ] || die "--archive $archive has no COMMIT file (send it with the git archive one-liner, docs/setup-guide.md B1)"
    printf 'KC_SRC_KIND=archive\nKC_SRC_DIR=%q\n' "$archive" | write_file "$RUN_HOME/source.env"
    say "== code: the archive tree $archive (commit $(cat "$archive/COMMIT")); no GitHub pushes in this mode"
  else
    if [ -z "$url" ]; then
      have gh || die "gh not found (it comes with the conda env)"
      if ! gh auth status --hostname github.com >/dev/null 2>&1; then
        [ -t 0 ] || die "not logged in to GitHub: run  bash run.sh setup  in an interactive ssh session once"
        say "== GitHub login: gh prints a one-time code; on your laptop open https://github.com/login/device and type it"
        run gh auth login --hostname github.com --git-protocol https --web
      fi
      run gh auth setup-git
      if [ "$DRY" = 0 ] && ! git config --global user.name >/dev/null; then
        git config --global user.name "$(gh api user --jq .login)"
        git config --global user.email "$(gh api user --jq .id)+$(gh api user --jq .login)@users.noreply.github.com"
      fi
    fi
    if [ ! -d "$RUN_HOME/repo/.git" ]; then
      say "== cloning ${url:-$repo} into $RUN_HOME/repo"
      if [ -n "$url" ]; then run git clone -q "$url" "$RUN_HOME/repo"; else run gh repo clone "$repo" "$RUN_HOME/repo" -- -q; fi
    fi
    printf 'KC_SRC_KIND=gh\nKC_SRC_DIR=%q\n' "$RUN_HOME/repo" | write_file "$RUN_HOME/source.env"
  fi
  say "setup: $(nproc_) cores, $(( $(mem_avail_kib) / 1048576 )) GB memory available; survival: $(survival_auto)"
  if [ -n "$commit" ]; then cmd_selftest --commit "$commit"
  else say "next: bash run.sh selftest --commit <tag or full hash>"; fi
}

# ------------------------------------------------------------------------------------------------- selftest
# make_tree COMMIT DIR: the files of COMMIT in DIR (no .git): git archive from the clone, or a copy of the archive
make_tree() {
  run mkdir -p "$2"
  if [ "$SRC_KIND" = archive ]; then run cp -a "$SRC_DIR/." "$2/"
  elif [ "$DRY" = 1 ]; then run git -C "$SRC_DIR" archive "$1"; say "[dry-run]   | tar -x -C $2"
  else git -C "$SRC_DIR" archive "$1" | tar -x -C "$2"; fi
}
cc_for_env() {  # the C compiler and FLINT prefix of the run: the conda env's, or the system's in test mode
  if [ -z "${KISSCRIT_PYTHON:-}" ] && compgen -G "$ENV_DIR/bin/*-gcc" >/dev/null; then
    local -a ccs=("$ENV_DIR"/bin/*-gcc); CC=${ccs[0]}; FLINT_PREFIX=$ENV_DIR
  else CC=$(command -v cc || command -v gcc || echo cc); FLINT_PREFIX=""; fi
}
cmd_selftest() {
  local ref="" subset="" commit sdir tree rec log py ok=1 rc out kib1 kib2 np nf ns
  while [ $# -gt 0 ]; do
    case $1 in --commit) ref=$2; shift ;; --pytest-subset) subset=$2; shift ;; *) die "selftest: unknown option $1" ;; esac; shift
  done
  [ -n "$ref" ] || die "selftest needs --commit <tag or full hash>"
  load_source; commit=$(resolve_commit "$ref") || exit 1
  sdir="$RUN_HOME/selftest/$commit"; tree="$sdir/tree"; rec="$RUN_HOME/selftest/$commit.json"
  if [ -z "${KC_INTREE:-}" ]; then
    run rm -rf "$tree"; make_tree "$commit" "$tree"
    if [ "$DRY" = 0 ]; then
      [ -f "$tree/scripts/full/run.sh" ] || die "commit $commit has no scripts/full/run.sh: pin a commit at or after F8"
      KC_INTREE=1 exec bash "$tree/scripts/full/run.sh" selftest --commit "$commit" ${subset:+--pytest-subset "$subset"}
    fi
  fi
  py=$(run_python) || die "no Python for the run: run setup first"
  log="$RUN_HOME/logs/selftest-${commit:0:12}.log"; run mkdir -p "$RUN_HOME/logs" "$RUN_HOME/tmp"
  [ "$DRY" = 1 ] && log=/dev/stdout     # a dry run prints the commands that would write to the log
  export TMPDIR="$RUN_HOME/tmp" NUMBA_CACHE_DIR="$RUN_HOME/cache/numba" PYTHONPYCACHEPREFIX="$RUN_HOME/cache/pyc"
  say "== self-test of $commit in $tree (log $log)"
  [ "$DRY" = 0 ] && { : > "$log"; st jset "$rec" commit "\"$commit\""; st jset "$rec" started_utc "\"$(utc)\""; st jset "$rec" passed false; }
  # 1. the C kernel: build with the run's compiler and the repo flags, then the smoke test
  cc_for_env
  if run make -C "$tree/csrc" "CC=$CC" ${FLINT_PREFIX:+"FLINT_PREFIX=$FLINT_PREFIX"} all >>"$log" 2>&1 \
     && run make -C "$tree/csrc" "CC=$CC" ${FLINT_PREFIX:+"FLINT_PREFIX=$FLINT_PREFIX"} test >>"$log" 2>&1 \
     && { [ "$DRY" = 1 ] || grep -q 'all containment tests passed' "$log"; }
  then say "  ok   C kernel built ($CC) and make test passed"; [ "$DRY" = 0 ] && st jset "$rec" parts.c_kernel '"passed"'
  else say "  FAIL C kernel build or make test (see $log)"; ok=0; [ "$DRY" = 0 ] && st jset "$rec" parts.c_kernel '"failed"'; fi
  # 2. unit tests: tests/test_ball_kill.py in its own process (it sets flint.ctx.dps globally, which changes the
  #    bits of later C/Arb comparisons in the same process; seen 2026-10-21 at cc375bc)
  [ "$DRY" = 0 ] && { mkdir -p "$RUN_HOME/bin"; printf '#!/bin/sh\nexec "%s" "$@"\n' "$CC" > "$RUN_HOME/bin/cc"; cp "$RUN_HOME/bin/cc" "$RUN_HOME/bin/gcc"; chmod +x "$RUN_HOME/bin/cc" "$RUN_HOME/bin/gcc"; }
  local -a ev=("PATH=$RUN_HOME/bin:$(dirname "$py"):$PATH" "KISSCRIT_REQUIRE_ALL=1")   # tests/conftest.py: a skipped test (except DOM1_REGEN) fails the self-test
  [ -n "$FLINT_PREFIX" ] && ev+=("CPATH=$FLINT_PREFIX/include" "LIBRARY_PATH=$FLINT_PREFIX/lib" "LD_LIBRARY_PATH=$FLINT_PREFIX/lib")
  # pytest groups, one process each: "a b c" is split into arguments on spaces
  local -a groups
  if [ -n "$subset" ]; then groups=("$subset")
  else groups=("tests --ignore=tests/test_ball_kill.py"); [ -f "$tree/tests/test_ball_kill.py" ] && groups+=("tests/test_ball_kill.py"); fi
  np=0; nf=0; ns=0
  local g; local -a gargs
  for g in "${groups[@]}"; do
    read -r -a gargs <<< "$g"
    if [ "$DRY" = 1 ]; then run env "${ev[@]}" "$py" -m pytest -q -p no:cacheprovider "${gargs[@]}"; continue; fi
    out=$(cd "$tree" && env "${ev[@]}" "$py" -m pytest -q -p no:cacheprovider "${gargs[@]}" 2>&1); rc=$?
    printf '%s\n' "$out" >> "$log"
    np=$((np + $(printf '%s' "$out" | grep -Eo '[0-9]+ passed' | tail -n1 | cut -d' ' -f1 || echo 0)))
    nf=$((nf + $(printf '%s' "$out" | grep -Eo '[0-9]+ (failed|errors?)' | awk '{s+=$1} END{print s+0}')))
    ns=$((ns + $(printf '%s' "$out" | grep -Eo '[0-9]+ skipped' | tail -n1 | cut -d' ' -f1 || echo 0)))
    [ $rc -eq 0 ] || ok=0
  done
  if [ "$DRY" = 0 ]; then
    if [ $ok = 1 ] && [ $nf -eq 0 ]; then say "  ok   pytest: $np passed, $ns skipped${subset:+ (subset: $subset)}"
    else say "  FAIL pytest: $np passed, $nf failed or errors, $ns skipped (see $log)"; ok=0; fi
    st jset "$rec" parts.pytest "{\"passed\": $np, \"failed\": $nf, \"skipped\": $ns, \"subset\": \"$subset\"}"
    st jset "$rec" subset "$([ -n "$subset" ] && echo true || echo false)"
  fi
  # 3. selftest.py (Lemma P, compiled kernels against the reference, one candidate)
  if (cd "$tree" && run "$py" scripts/completeness/selftest.py) >>"$log" 2>&1; then
    say "  ok   scripts/completeness/selftest.py"; [ "$DRY" = 0 ] && st jset "$rec" parts.selftest_py '"passed"'
  else say "  FAIL scripts/completeness/selftest.py (see $log)"; ok=0; [ "$DRY" = 0 ] && st jset "$rec" parts.selftest_py '"failed"'; fi
  # 4. peak RSS per certify worker: route test, one certify shard of 8 with one worker, at the default budget
  if [ "$DRY" = 1 ]; then
    run "$py" -I "$tree/scripts/full/run_state.py" rss -- "$py" scripts/completeness/run_shard.py certify --route test --shard 0 --of 8 --workers 1 --budget 300000
  else
    rm -rf "$tree/results/completeness"
    out=$(cd "$tree" && "$py" -I scripts/full/run_state.py rss -- "$py" scripts/completeness/run_shard.py gen --route test --shard 0 --of 1 2>>"$log" | tail -n 1)
    kib1=$(printf '%s' "$out" | "$py" -c 'import json,sys; print(json.load(sys.stdin)["peak_rss_kib"])' 2>/dev/null || echo 0)
    out=$(cd "$tree" && "$py" -I scripts/full/run_state.py rss -- "$py" scripts/completeness/run_shard.py certify --route test --shard 0 --of 8 --workers 1 --budget 300000 2>>"$log" | tail -n 1)
    kib2=$(printf '%s' "$out" | "$py" -c 'import json,sys; d=json.load(sys.stdin); print(d["peak_rss_kib"] if d["exit"]==0 else 0)' 2>/dev/null || echo 0)
    if [ "$kib2" -gt 0 ]; then
      say "  ok   peak RSS: gen $((kib1 / 1024)) MiB, one certify worker $((kib2 / 1024)) MiB (route test, shard 0 of 8)"
      st jset "$rec" peak_rss_kib "$(( kib1 > kib2 ? kib1 : kib2 ))"
    else say "  FAIL peak RSS measurement (route test certify failed; see $log)"; ok=0; fi
  fi
  # 5. Zenodo sandbox round trip: milestone M3
  if [ -e "$TOKDIR/zenodo-sandbox.token" ]; then say "  --   Zenodo sandbox round trip: not in this version (F8 milestone M3)"
  else say "  --   Zenodo sandbox round trip: no sandbox token, not run"; fi
  [ "$DRY" = 1 ] && return 0
  st jset "$rec" parts.zenodo_sandbox '"not run (M3)"'
  st jset "$rec" versions "$("$py" -I "$tree/scripts/full/run_state.py" versions)"
  st jset "$rec" host "\"$(hostname)\""
  st jset "$rec" finished_utc "\"$(utc)\""
  if [ $ok = 1 ]; then st jset "$rec" passed true; say "selftest passed${subset:+ (pytest subset only)} for $commit"; return 0
  else say "selftest FAILED for $commit (log $log)"; return 1; fi
}

# ------------------------------------------------------------------------------------------------- survival
survival_auto() {
  if [ "$(os_name)" != Linux ]; then echo "none (macOS launchd is milestone M4)"; return; fi
  local linger; linger=$("$LOGINCTL" show-user "$(id -un)" -p Linger --value 2>/dev/null || echo unknown)
  if [ "$linger" = yes ] && "$SYSTEMCTL" --user show-environment >/dev/null 2>&1; then echo systemd
  elif have tmux; then echo tmux
  else echo "none (no lingering and no tmux)"; fi
}
unit_text() {
  cat <<EOF
# written by scripts/full/run.sh start (run $RUN_ID); removed by nothing: run.sh stop disables it
[Unit]
Description=kisscrit completeness run $RUN_ID (scripts/full/run.sh)
StartLimitIntervalSec=3600
StartLimitBurst=5

[Service]
Type=simple
ExecStart=/bin/bash $TREE/scripts/full/run.sh _run
WorkingDirectory=$RUN_HOME
Environment=KISSCRIT_RUN_HOME=$RUN_HOME
StandardOutput=append:$RUN_HOME/logs/run.log
StandardError=append:$RUN_HOME/logs/run.log
Restart=on-failure
RestartSec=120
Nice=5
KillMode=control-group
TimeoutStopSec=90

[Install]
WantedBy=default.target
EOF
}
# cron_usable: crontab exists and answers -l with a table or "no crontab for ..." (not "not allowed")
cron_usable() { local out; have "$CRONTAB" || return 1; out=$("$CRONTAB" -l 2>&1) && return 0; printf '%s' "$out" | grep -qi 'no crontab'; }
cron_line() { printf 'KISSCRIT_RUN_HOME=%q /bin/bash %q start --boot >> %q 2>&1 %s' "$RUN_HOME" "$TREE/scripts/full/run.sh" "$RUN_HOME/logs/boot.log" "$CRON_TAG"; }
launch() {
  local mode=$1 inner
  inner="KISSCRIT_RUN_HOME=$(printf %q "$RUN_HOME") exec /bin/bash $(printf %q "$TREE/scripts/full/run.sh") _run >> $(printf %q "$RUN_HOME/logs/run.log") 2>&1"
  case $mode in
    systemd)
      unit_text | write_file "$UNIT_FILE"
      run "$SYSTEMCTL" --user daemon-reload
      run "$SYSTEMCTL" --user enable "$UNIT"
      run "$SYSTEMCTL" --user restart "$UNIT"
      say "started as the systemd user service $UNIT (lingering is on: it survives logout and restarts after a reboot)" ;;
    tmux)
      tmux -L "$TMUX_SOCK" has-session -t "$TMUX_SESSION" 2>/dev/null && run tmux -L "$TMUX_SOCK" kill-session -t "$TMUX_SESSION"
      run tmux -L "$TMUX_SOCK" new-session -d -s "$TMUX_SESSION" "$inner"
      if cron_usable; then
        if [ "$DRY" = 1 ]; then printf '[dry-run] crontab: replace the tagged line by: %s\n' "$(cron_line)"
        else { "$CRONTAB" -l 2>/dev/null | grep -vF "$CRON_TAG"; cron_line; printf '\n'; } | "$CRONTAB" -; fi
        say "started in the detached tmux session $TMUX_SESSION (socket $TMUX_SOCK); an @reboot crontab line restarts it after a reboot"
      else
        say "started in the detached tmux session $TMUX_SESSION (socket $TMUX_SOCK); WARNING: crontab is not usable, so the run does NOT restart after a reboot (run.sh start resumes it by hand)"
      fi ;;
    foreground)
      say "running in the foreground (log $RUN_HOME/logs/run.log)"
      if [ "$DRY" = 1 ]; then run /bin/bash "$TREE/scripts/full/run.sh" _run
      else KISSCRIT_RUN_HOME="$RUN_HOME" /bin/bash "$TREE/scripts/full/run.sh" _run 2>&1 | tee -a "$RUN_HOME/logs/run.log"; return "${PIPESTATUS[0]}"; fi ;;
    *) die "no way to keep the run alive: $mode (docs/setup-guide.md A4)" ;;
  esac
}

# ------------------------------------------------------------------------------------------------- start
cmd_start() {
  local route="" ref="" run_id="" fresh=0 shards=40 budget=300000 workers=0 wps=2 checker=$DEFAULT_CHECKER only=""
  local need=50 survival=auto new=0 boot=0 allow_partial=0 poll=30 every=3600 prepare=0 commit outs n f
  local -a pass=("$@")
  while [ $# -gt 0 ]; do
    case $1 in
      --route) route=$2; shift ;; --commit) ref=$2; shift ;; --run-id) run_id=$2; shift ;; --fresh) fresh=1 ;;
      --shards) shards=$2; shift ;; --budget) budget=$2; shift ;; --workers) workers=$2; shift ;;
      --workers-per-shard) wps=$2; shift ;; --checker) checker=$2; shift ;; --only-shards) only=$2; shift ;;
      --need-disk-gb) need=$2; shift ;; --survival) survival=$2; shift ;; --new) new=1 ;; --boot) boot=1 ;;
      --allow-partial-selftest) allow_partial=1 ;; --poll) poll=$2; shift ;; --progress-every) every=$2; shift ;;
      --prepare-only) prepare=1 ;;
      *) die "start: unknown option $1" ;;
    esac; shift
  done
  [ "$(os_name)" = Linux ] || die "this version supports Linux only; macOS (launchd) is milestone M4"
  run mkdir -p "$RUN_HOME/logs" "$RUN_HOME/runs"
  # ---- resume the current run
  if [ -z "$route" ]; then
    [ -z "$ref" ] || die "start: --commit needs --route (or give neither to resume the current run)"
    load_current || { [ $boot = 1 ] && exit 0; die "no current run: start a run with --route and --commit"; }
    if [ $boot = 1 ] && [ -f "$(stopped_flag)" ]; then say "$(utc) boot: run $RUN_ID was stopped by hand; not restarted"; exit 0; fi
    [ -f "$RDIR/state/done" ] && { say "run $RUN_ID is finished ($(cat "$RDIR/state/done")); start a new one with --route --commit --new"; exit 0; }
    if [ -z "${KC_INTREE:-}" ] && [ "$DRY" = 0 ] && [ "$SELF" != "$TREE/scripts/full/run.sh" ]; then
      KC_INTREE=1 exec bash "$TREE/scripts/full/run.sh" start "${pass[@]}"
    fi
    lock_free || { say "run $RUN_ID is already running (bash run.sh status)"; exit 0; }
    run rm -f "$(stopped_flag)"
    [ $boot = 1 ] && say "$(utc) boot: resuming run $RUN_ID"
    launch "$KC_SURVIVAL"; return
  fi
  # ---- a new run
  case $route in test|cond) ;; full) ;; *) die "start: --route must be test, cond or full" ;; esac
  [ -n "$ref" ] || die "start: a new run needs --commit <tag or full hash> (the pinned commit)"
  [[ $shards =~ ^[1-9][0-9]*$ && $budget =~ ^[1-9][0-9]*$ && $workers =~ ^[0-9]+$ && $wps =~ ^[1-9][0-9]*$ ]] || die "start: --shards, --budget, --workers-per-shard must be positive integers"
  if [ -n "$only" ]; then
    for n in ${only//,/ }; do
      if ! [[ $n =~ ^[0-9]+$ ]] || [ "$n" -ge "$shards" ]; then die "start: --only-shards: $n is not a shard of 0..$((shards-1))"; fi
    done
  fi
  if load_current && [ ! -f "$RDIR/state/done" ] && [ "${KC_INTREE:-}" != 1 ]; then
    [ $new = 1 ] || die "run $RUN_ID is not finished: 'bash run.sh start' resumes it; add --new to begin another (it is kept in $RDIR)"
    lock_free || die "run $RUN_ID is running: bash run.sh stop first"
  fi
  load_source
  commit=$(resolve_commit "$ref") || exit 1
  if [ "$route" = full ]; then
    if [ "$SRC_KIND" = archive ]; then f=$([ -f "$SRC_DIR/results/full/F6-READY" ] && echo y)
    else f=$(git -C "$SRC_DIR" cat-file -e "$commit:results/full/F6-READY" 2>/dev/null && echo y); fi
    [ "$f" = y ] || die "route full needs checklist F6 (gen for full), which is not done: commit $commit has no results/full/F6-READY. Use --route cond or test."
  fi
  if [ "$SRC_KIND" = gh ] && code_dirty "$SRC_DIR"; then die "the clone $SRC_DIR has uncommitted code changes; the run needs a clean pinned commit"; fi
  [ "$(free_gb "$RUN_HOME")" -ge "$need" ] || die "only $(free_gb "$RUN_HOME") GB free under $RUN_HOME (want $need; --need-disk-gb)"
  [ -n "$run_id" ] || run_id="$(date -u +%Y%m%dT%H%MZ)-${commit:0:7}"
  [[ $run_id =~ ^[A-Za-z0-9._-]+$ ]] || die "start: --run-id may contain only letters, digits, . _ -"
  set_run "$run_id"; TREE="$RDIR/tree"; local branch="run/$route-$run_id"
  # selftest of this commit
  local srec="$RUN_HOME/selftest/$commit.json"
  if [ "$DRY" = 0 ]; then
    [ "$(st jget "$srec" passed)" = true ] || die "no passed self-test for $commit: run  bash run.sh selftest --commit $commit  first"
    if [ "$(st jget "$srec" subset)" = true ] && [ $allow_partial = 0 ]; then die "the self-test of $commit ran a pytest subset only; rerun it in full (or --allow-partial-selftest, recorded)"; fi
  fi
  # the run tree: a worktree on the run branch (gh), or a git repository made from the archive tree
  if [ -z "${KC_INTREE:-}" ]; then
    [ -e "$TREE" ] && die "$TREE exists already (choose another --run-id)"
    run mkdir -p "$RDIR/state" "$RDIR/logs" "$RDIR/quarantine"
    if [ "$SRC_KIND" = gh ]; then
      run git -C "$SRC_DIR" worktree add -q -b "$branch" "$TREE" "$commit" || die "could not create the worktree"
    else
      make_tree "$commit" "$TREE"
      run git -C "$TREE" init -q
      run git -C "$TREE" checkout -q -b "$branch"
      run git -C "$TREE" add -A
      gitc commit -q -m "archive of $commit (run $run_id; tree sent with git archive)"
    fi
    if [ "$DRY" = 0 ]; then
      [ -f "$TREE/scripts/full/run.sh" ] || die "commit $commit has no scripts/full/run.sh: pin a commit at or after F8"
      KC_INTREE=1 exec bash "$TREE/scripts/full/run.sh" start "${pass[@]}" --run-id "$run_id"
    fi
  fi
  # ---- from here on: the run tree's own copy of run.sh
  outs=$(st outputs "$TREE" "$route" 2>/dev/null || true)
  if [ -n "$outs" ]; then
    if [ $fresh = 0 ]; then
      say "The pinned tree already holds results for route $route ($(printf '%s\n' "$outs" | wc -l) files, e.g. $(printf '%s\n' "$outs" | head -n 1))."
      say "run_shard.py would treat those candidates as done.  Add --fresh to remove them on the run branch $branch only (one logged commit)."
      if [ "$SRC_KIND" = gh ]; then run git -C "$SRC_DIR" worktree remove --force "$TREE"; run git -C "$SRC_DIR" branch -q -D "$branch"; fi
      run rm -rf "$RDIR"
      exit 2
    fi
    printf '%s\n' "$outs" | write_file "$RDIR/logs/fresh-removed.txt"
    if [ "$DRY" = 0 ]; then
      while IFS= read -r f; do git -C "$TREE" rm -q --cached --ignore-unmatch -- "$f" >/dev/null; rm -f "$TREE/$f"; done <<< "$outs"
    else printf '[dry-run] git rm and rm the %d files listed above\n' "$(printf '%s\n' "$outs" | wc -l)"; fi
    gitc commit -q --allow-empty -m "run $run_id: --fresh: remove $(printf '%s\n' "$outs" | wc -l) existing $route result files from the run branch" -m "$(printf '%s\n' "$outs" | head -n 200)"
  fi
  [ "$survival" = auto ] && survival=$(survival_auto)
  case $survival in systemd|tmux|foreground) ;; *) die "no way to keep the run alive: $survival (docs/setup-guide.md A4)" ;; esac
  local py; py=$(run_python) || die "no Python for the run: run setup first"
  local gshards; case $route in test) gshards=1 ;; cond) gshards=4 ;; full) gshards=40 ;; esac
  local ccommit cver
  if [ "$SRC_KIND" = gh ]; then ccommit=$(git -C "$TREE" log -1 --format=%H "$commit" -- kisscrit/core "$checker" 2>/dev/null)
  else ccommit=$commit; fi
  [ -n "$ccommit" ] || ccommit=$commit
  if [ "$DRY" = 0 ]; then
    git -C "$TREE" ls-files --error-unmatch -- "$checker" >/dev/null 2>&1 || die "--checker $checker is not a tracked file of the pinned tree (only versioned scripts run)"
  fi
  cver=$(sed -n "s/^CHECKER_VERSION = '\([^']*\)'.*/\1/p" "$TREE/$checker" 2>/dev/null | head -n 1)
  local only_json="[]"; [ -n "$only" ] && only_json="[${only}]"
  # the configuration (local) and the run record (results/full/run-record.json, committed)
  if [ "$DRY" = 0 ]; then
    rm -f "$CFG"
    st jset "$CFG" run_id "\"$run_id\""; st jset "$CFG" route "\"$route\""; st jset "$CFG" ref "\"$ref\""
    st jset "$CFG" commit "\"$commit\""; st jset "$CFG" source "\"$SRC_KIND\""; st jset "$CFG" branch "\"$branch\""
    st jset "$CFG" tree "\"$TREE\""; st jset "$CFG" shards "$shards"; st jset "$CFG" gen_shards "$gshards"
    st jset "$CFG" budget "$budget"; st jset "$CFG" workers "$workers"; st jset "$CFG" wps "$wps"
    st jset "$CFG" checker_cmd "\"$checker\""; st jset "$CFG" checker_default "$([ "$checker" = "$DEFAULT_CHECKER" ] && echo true || echo false)"
    st jset "$CFG" only_shards "$only_json"; st jset "$CFG" survival "\"$survival\""; st jset "$CFG" python "\"$py\""
    st jset "$CFG" poll "$poll"; st jset "$CFG" progress_every "$every"; st jset "$CFG" need_disk_gb "$need"; st jset "$CFG" push 0
    st jset "$CFG" checker_commit "\"$ccommit\""
    st jset "$CFG" verdict_label "\"checker verdict at commit $ccommit; not certified\""
    st jset "$CFG" created_utc "\"$(utc)\""
    printf '%s\n' "$run_id" > "$CUR_FILE"
    write_record "$srec" "$cver" "$allow_partial"
    commit_results "run $run_id: run record (route $route, commit ${commit:0:12}, checker $cver)" results/full
  else
    say "[dry-run] would write $CFG and $TREE/results/full/run-record.json and set $CUR_FILE to $run_id"
  fi
  eval "$(st env "$CFG" 2>/dev/null || true)"
  say "run $run_id: route $route, commit $commit, branch $branch, tree $TREE, survival $survival"
  if [ $prepare = 1 ]; then say "prepared, not started (--prepare-only):  bash run.sh start  launches it"; return 0; fi
  launch "$survival"
}
write_record() {  # write_record SELFTEST_RECORD CHECKER_NAME ALLOW_PARTIAL: the run record, from the config
  local r; r=$(REC); local py; py=$(run_python)
  "$HPY" -I - "$CFG" "$r" <<'PY'
import json, sys
c = json.load(open(sys.argv[1])); out = sys.argv[2]
rec = dict(kind='kisscrit run record (scripts/full/run.sh)', **{k: c[k] for k in sorted(c)})
import os; os.makedirs(os.path.dirname(out), exist_ok=True)
json.dump(rec, open(out, 'w'), indent=1, sort_keys=True)
PY
  # shellcheck disable=SC1091
  st jset "$r" host "{\"name\": \"$(hostname)\", \"nproc\": $(nproc_), \"mem_total_kib\": $(awk '/MemTotal/{print $2}' /proc/meminfo), \"os\": \"$( (. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-unknown}") )\", \"kernel\": \"$(uname -r)\", \"user\": \"$(id -un)\"}"
  st jset "$r" versions "$("$py" -I "$HELPER" versions)"
  st jset "$r" checker "{\"script\": \"$(st jget "$CFG" checker_cmd)\", \"name\": \"$2\", \"commit\": \"$(st jget "$CFG" checker_commit)\", \"default\": $(st jget "$CFG" checker_default), \"note\": \"the checker that decides; the C kernel front end (kisscrit/ckernel.py) is not used (audit B-8)\"}"
  st jset "$r" selftest "$(cat "$1")"
  st jset "$r" allow_partial_selftest "$([ "$3" = 1 ] && echo true || echo false)"
  mkdir -p "$TREE/results/full"
  if [ -z "${KISSCRIT_PYTHON:-}" ] && [ -x "$MF/bin/conda" ]; then
    "$MF/bin/conda" list -p "$ENV_DIR" --export > "$TREE/results/full/conda-list.txt" 2>/dev/null
    st jset "$r" package_list '"conda list --export (results/full/conda-list.txt)"'
  else
    "$py" -m pip freeze > "$TREE/results/full/conda-list.txt" 2>/dev/null
    st jset "$r" package_list '"pip freeze: KISSCRIT_PYTHON test mode, no conda env (results/full/conda-list.txt)"'
  fi
  st jset "$r" package_list_sha256 "\"$(st sha256 "$TREE/results/full/conda-list.txt" | cut -d' ' -f1)\""
}

# commit_results MESSAGE PATH...: commit result paths on the run branch (never code, never the certificates)
commit_results() {
  local msg=$1 p; shift; local -a add=()
  for p in "$@"; do
    case $p in results/completeness/certificates*) die "internal: certificates are never committed ($p)" ;; results/*) ;; *) die "internal: commit outside results/ refused ($p)" ;; esac
    [ -e "$TREE/$p" ] && add+=("$p")
  done
  [ ${#add[@]} -gt 0 ] || return 0
  run git -C "$TREE" add -f -- "${add[@]}"
  if [ "$DRY" = 0 ] && git -C "$TREE" diff --cached --quiet; then return 0; fi
  gitc commit -q -m "$msg" -m "Run tooling: scripts/full/run.sh (F8). Checker verdicts are labelled: $(st jget "$CFG" verdict_label)."
}

# ------------------------------------------------------------------------------------------------- the run itself
LOG() { printf '%s  %s\n' "$(date '+%F %T %Z')" "$*"; }
fail_phase() { LOG "FAILED: $*"; st jset "$STATE" last_failure "\"$(utc) $*\""; exit 1; }
phase_begin() {
  LOG "PHASE $1 start"; st jset "$STATE" phase "\"$1\""
  [ -n "$(st jget "$STATE" "phases.$1.first_start_epoch")" ] || st jset "$STATE" "phases.$1.first_start_epoch" "$(date +%s)"
  PHASE_T0=$(date +%s)
}
phase_end() {
  local prev; prev=$(st jget "$STATE" "phases.$1.wall_seconds"); prev=${prev:-0}
  st jset "$STATE" "phases.$1.wall_seconds" "$(( prev + $(date +%s) - PHASE_T0 ))"
  st jset "$STATE" "phases.$1.end_utc" "\"$(utc)\""
  st jset "$(REC)" "phases.$1" "$(st jget "$STATE" "phases.$1")"
  LOG "PHASE $1 done"
}
shard_list() { if [ -n "$KC_ONLY_SHARDS" ]; then printf '%s\n' "${KC_ONLY_SHARDS// /$'\n'}"; else seq 0 $((KC_SHARDS - 1)); fi; }
progress() { st progress --config "$CFG" --state "$STATE" --out "$TREE/results/full/progress.json"; }

cmd_run_inner() {
  load_current || die "_run: no current run"
  local PY=$KC_PYTHON
  mkdir -p "$RDIR/state" "$RDIR/logs" "$RDIR/quarantine" "$RUN_HOME/tmp" "$RUN_HOME/cache"
  exec 9>>"$RUN_HOME/run.lock"
  if ! flock -n 9; then LOG "another copy holds $RUN_HOME/run.lock (processes of an earlier start may still be running); waiting"; flock 9; fi
  [ -f "$(stopped_flag)" ] && { LOG "run $RUN_ID is stopped ($(cat "$(stopped_flag)")); bash run.sh start resumes it"; exit 0; }
  [ -f "$RDIR/state/done" ] && { LOG "run $RUN_ID is finished"; exit 0; }
  local PGID; PGID=$(ps -o pgid= $$ | tr -d ' ')
  printf '%s %s %s\n' "$$" "$PGID" "$(hostname)" > "$RDIR/state/pid"
  trap 'LOG "terminated by a signal"; [ "$PGID" = "$$" ] && kill -TERM -- -"$PGID" 2>/dev/null; kill $(jobs -p) 2>/dev/null; exit 143' TERM INT HUP
  ulimit -Sn "$(nofile_target)" 2>/dev/null || LOG "could not raise the open-files limit (soft $(ulimit -Sn))"
  export TMPDIR="$RUN_HOME/tmp" NUMBA_CACHE_DIR="$RUN_HOME/cache/numba" PYTHONPYCACHEPREFIX="$RUN_HOME/cache/pyc" MPLCONFIGDIR="$RUN_HOME/cache/mpl"
  cd "$TREE" || fail_phase "cannot cd to $TREE"
  code_dirty "$TREE" && fail_phase "the code in $TREE differs from commit $KC_COMMIT (git status); a running computation's scripts must not change"
  [ "$(git -C "$TREE" rev-parse --abbrev-ref HEAD)" = "$KC_BRANCH" ] || fail_phase "$TREE is not on the run branch $KC_BRANCH"
  [ "$(free_gb "$RUN_HOME")" -ge 1 ] || fail_phase "less than 1 GB free under $RUN_HOME"
  LOG "run $RUN_ID: route $KC_ROUTE, commit $KC_COMMIT, pid $$, open files $(ulimit -Sn), host $(hostname)"
  st jappend "$(REC)" resumes "{\"utc\": \"$(utc)\", \"pid\": $$, \"survival\": \"$KC_SURVIVAL\", \"open_files_soft\": \"$(ulimit -Sn)\"}"
  # workers: requested (default nproc), capped by memory and the measured peak RSS per worker
  local rss W P w bym
  rss=$(st jget "$RUN_HOME/selftest/$KC_COMMIT.json" peak_rss_kib); rss=${rss:-1048576}
  read -r W P w bym < <(st workers --nproc "$(nproc_)" --mem-avail-kib "$(mem_avail_kib)" --rss-kib "$rss" --wps "$KC_WPS" \
                          --requested "$KC_WORKERS" --shards "$(shard_list | wc -l)")
  st jset "$CFG" concurrent_shards "$P"; st jset "$CFG" workers_used "$W"
  st jset "$(REC)" resources "{\"nproc\": $(nproc_), \"mem_available_kib\": $(mem_avail_kib), \"peak_rss_kib_selftest\": $rss, \"workers\": $W, \"workers_by_memory\": $bym, \"concurrent_shards\": $P, \"workers_per_shard\": $w}"
  LOG "workers $W (memory allows $bym), $P shards at once, $w workers per shard, peak RSS per worker from the self-test $((rss / 1024)) MiB"
  # repair what a kill may have left (partial lines, a certificate without its result line)
  "$HPY" -I "$HERE/run_repair.py" "$TREE/results" "$KC_ROUTE" "$RDIR/quarantine" --log "$RDIR/logs/repair.jsonl" \
    | while read -r l; do LOG "$l"; done
  [ "${PIPESTATUS[0]}" -eq 0 ] || fail_phase "run_repair.py failed"

  # ---- gen
  phase_begin gen
  local k pids=() rc
  for k in $(seq 0 $((KC_GEN_SHARDS - 1))); do
    [ -f "$RDIR/state/gen-$k.done" ] && continue
    "$PY" scripts/completeness/run_shard.py gen --route "$KC_ROUTE" --shard "$k" --of "$KC_GEN_SHARDS" > "$RDIR/logs/gen-$k.log" 2>&1 &
    pids+=("$!:$k")
  done
  for p in "${pids[@]+"${pids[@]}"}"; do
    wait "${p%%:*}"; rc=$?; k=${p#*:}
    [ $rc -eq 0 ] || fail_phase "gen shard $k exited $rc (log $RDIR/logs/gen-$k.log)"
    touch "$RDIR/state/gen-$k.done"
  done
  commit_results "run $RUN_ID: gen $KC_ROUTE ($KC_GEN_SHARDS shards)" results/completeness/gen
  phase_end gen

  # ---- certify: shards concurrently, P at a time, w workers each; a commit per finished shard
  phase_begin certify
  local -A running=(); local todo=() last=0 failed=0 pid
  for k in $(shard_list); do [ -f "results/completeness/cert/done-$KC_ROUTE-$k-of-$KC_SHARDS.txt" ] || todo+=("$k"); done
  LOG "certify: ${#todo[@]} shards to do of $(shard_list | wc -l), budget $KC_BUDGET"
  while [ ${#todo[@]} -gt 0 ] || [ ${#running[@]} -gt 0 ]; do
    while [ ${#todo[@]} -gt 0 ] && [ ${#running[@]} -lt "$P" ]; do
      k=${todo[0]}; todo=("${todo[@]:1}")
      "$PY" scripts/completeness/run_shard.py certify --route "$KC_ROUTE" --shard "$k" --of "$KC_SHARDS" \
            --workers "$w" --budget "$KC_BUDGET" >> "$RDIR/logs/cert-$k.log" 2>&1 &
      running[$!]=$k; LOG "certify shard $k/$KC_SHARDS started (pid $!)"
    done
    sleep "$KC_POLL"
    for pid in "${!running[@]}"; do
      kill -0 "$pid" 2>/dev/null && continue
      wait "$pid"; rc=$?; k=${running[$pid]}; unset "running[$pid]"
      if [ $rc -eq 0 ] && [ -f "results/completeness/cert/done-$KC_ROUTE-$k-of-$KC_SHARDS.txt" ]; then
        progress >/dev/null
        commit_results "run $RUN_ID: certify $KC_ROUTE shard $k/$KC_SHARDS" "results/completeness/cert/cert-$KC_ROUTE-$k-of-$KC_SHARDS.jsonl" \
          "results/completeness/cert/done-$KC_ROUTE-$k-of-$KC_SHARDS.txt" results/full/progress.json
        LOG "certify shard $k/$KC_SHARDS done; $(progress)"
      else LOG "certify shard $k/$KC_SHARDS exited $rc (log $RDIR/logs/cert-$k.log)"; failed=$((failed + 1)); fi
    done
    if [ $(( $(date +%s) - last )) -ge "$KC_PROGRESS_EVERY" ]; then
      progress >/dev/null; commit_results "run $RUN_ID: progress" results/full/progress.json; last=$(date +%s)
    fi
  done
  [ $failed -eq 0 ] || fail_phase "$failed certify shards failed; rerunning start resumes them"
  phase_end certify
  "$PY" scripts/completeness/merge.py --route "$KC_ROUTE" --of "$KC_SHARDS" > "$RDIR/logs/merge.log" 2>&1 || fail_phase "merge.py failed ($RDIR/logs/merge.log)"
  progress >/dev/null
  report interim
  commit_results "run $RUN_ID: certify finished; merge summary and interim report" "results/completeness/summary-$KC_ROUTE.json" \
    results/full/progress.json results/full/report-interim.txt results/full/run-record.json

  # ---- check: the configured checker on each certificate file; a commit per file
  phase_begin check
  local cf name nv nt
  for k in $(shard_list); do
    cf="results/completeness/certificates/certificates-$KC_ROUTE-$k-of-$KC_SHARDS.jsonl"
    [ -f "$RDIR/state/check-$k.done" ] && continue
    [ -f "$cf" ] || fail_phase "missing $cf"
    LOG "check $cf with $KC_CHECKER_CMD ($W processes)"
    "$PY" "$KC_CHECKER_CMD" --route "$KC_ROUTE" --workers "$W" "$cf" >> "$RDIR/logs/check-$k.log" 2>&1 || fail_phase "checker exited non-zero on $cf (log $RDIR/logs/check-$k.log)"
    touch "$RDIR/state/check-$k.done"
    name=""; for f in results/completeness/checks/*-"$(basename "$cf")"; do [ -f "$f" ] && { name=$f; break; }; done
    nt=$(wc -l < "${name:-/dev/null}"); nv=$(grep -c '"verified": true' "${name:-/dev/null}")
    progress >/dev/null
    commit_results "run $RUN_ID: check $KC_ROUTE shard $k/$KC_SHARDS ($nv of $nt verified by $(basename "$KC_CHECKER_CMD"))" ${name:+"$name"} results/full/progress.json
  done
  phase_end check

  # ---- manifest (V5) and summary (V7), with the same checker script
  phase_begin manifest
  local -a files=(); for k in $(shard_list); do files+=("results/completeness/certificates/certificates-$KC_ROUTE-$k-of-$KC_SHARDS.jsonl"); done
  "$PY" "$KC_CHECKER_CMD" --route "$KC_ROUTE" --summary "${files[@]}" > results/full/check-summary.json 2>>"$RDIR/logs/manifest.log" || fail_phase "checker --summary failed"
  "$PY" "$KC_CHECKER_CMD" --route "$KC_ROUTE" --manifest > results/full/manifest.json 2>>"$RDIR/logs/manifest.log" || fail_phase "checker --manifest failed"
  phase_end manifest
  st jset "$STATE" phase '"done"'; progress >/dev/null
  report final
  commit_results "run $RUN_ID: manifest, check summary, final report ($(st jget results/full/manifest.json complete | sed 's/^/complete=/'))" \
    results/full/manifest.json results/full/check-summary.json results/full/report.txt results/full/progress.json results/full/run-record.json
  utc > "$RDIR/state/done"
  LOG "run $RUN_ID finished: $TREE/results/full/report.txt ($(st jget "$CFG" verdict_label))"
  [ "$KC_SURVIVAL" = systemd ] && "$SYSTEMCTL" --user disable "$UNIT" >/dev/null 2>&1
  return 0
}

report() {  # report interim|final: a plain-text report, local and in results/full
  local kind=$1 out f; [ "$kind" = interim ] && f=report-interim.txt || f=report.txt
  out="$TREE/results/full/$f"
  {
    echo "kisscrit run $RUN_ID: $kind report, written $(date '+%F %T %Z')"
    echo "verdict label: $(st jget "$CFG" verdict_label)"
    echo "route $KC_ROUTE; commit $KC_COMMIT (ref $KC_REF, source $KC_SOURCE); branch $KC_BRANCH; host $(hostname); $(nproc_) cores"
    echo "checker: $KC_CHECKER_CMD ($([ "$KC_CHECKER_DEFAULT" = 1 ] && echo "the default reference checker" || echo "NOT the default checker")), checker commit $(st jget "$CFG" checker_commit)"
    echo "settings: shards $KC_SHARDS${KC_ONLY_SHARDS:+ (only $KC_ONLY_SHARDS)}, gen shards $KC_GEN_SHARDS, budget $KC_BUDGET, workers $(st jget "$CFG" workers_used), shards at once $(st jget "$CFG" concurrent_shards), workers per shard $KC_WPS"
    echo "wall seconds: gen $(st jget "$STATE" phases.gen.wall_seconds), certify $(st jget "$STATE" phases.certify.wall_seconds), check $(st jget "$STATE" phases.check.wall_seconds)"
    echo "sizes (this route only): $(st sizes "$TREE" "$KC_ROUTE")"
    echo "merge.py:"; sed 's/^/  /' "$RDIR/logs/merge.log" 2>/dev/null
    echo "progress:"; sed 's/^/  /' results/full/progress.json
    if [ "$kind" = final ]; then
      echo "checker summary (V7):"; sed 's/^/  /' results/full/check-summary.json
      echo "manifest (V5):"; sed 's/^/  /' results/full/manifest.json
    else echo "the checker has not run yet: this report is written before the check phase, which can take hours"; fi
    echo "repairs before resume: $(wc -l < "$RDIR/logs/repair.jsonl" 2>/dev/null || echo 0) (details: $RDIR/logs/repair.jsonl)"
  } > "$out" 2>&1
  cp "$out" "$RDIR/$f"
  LOG "$kind report: $out"
}

# ------------------------------------------------------------------------------------------------- status, stop
# group_live PGID: the process group has a member that is not a zombie (kill -0 and pgrep count zombies, which a
# container's init may reap late)
# shellcheck disable=SC2009  # ps shows the state (Z); pgrep does not filter zombies everywhere
group_live() { [ -n "$1" ] && ps -o stat= -g "$1" 2>/dev/null | grep -qv '^Z'; }
alive() {  # alive: the run's processes exist (by the pid file's process group)
  local pid pg; [ -f "$RDIR/state/pid" ] || return 1; read -r pid pg _ < "$RDIR/state/pid"
  group_live "$pg" || { ps -o stat= -p "$pid" 2>/dev/null | grep -qv '^Z'; }
}
cmd_status() {
  local final=0; [ "${1:-}" = --final ] && final=1
  load_current || { say "no run yet (bash run.sh start --route ... --commit ...)"; return 0; }
  if [ -z "${KC_INTREE:-}" ] && [ "$DRY" = 0 ] && [ -f "$TREE/scripts/full/run.sh" ] && [ "$SELF" != "$TREE/scripts/full/run.sh" ]; then
    KC_INTREE=1 exec bash "$TREE/scripts/full/run.sh" status "$@"
  fi
  say "run $RUN_ID: route $KC_ROUTE, commit $KC_COMMIT, branch $KC_BRANCH"
  say "tree: $TREE"
  local surv="$KC_SURVIVAL"
  case $surv in
    systemd) surv="systemd user service $UNIT: $("$SYSTEMCTL" --user is-active "$UNIT" 2>/dev/null || echo inactive)" ;;
    tmux) surv="tmux session $TMUX_SESSION (socket $TMUX_SOCK): $(tmux -L "$TMUX_SOCK" has-session -t "$TMUX_SESSION" 2>/dev/null && echo running || echo absent); @reboot crontab line: $("$CRONTAB" -l 2>/dev/null | grep -qF "$CRON_TAG" && echo present || echo absent)" ;;
  esac
  say "survival: $surv"
  if alive; then say "processes: running (pid file $RDIR/state/pid)"; else say "processes: none"; fi
  if lock_free; then say "lock: free"; else say "lock: held"; fi
  [ -f "$(stopped_flag)" ] && say "stopped: $(cat "$(stopped_flag)")"
  [ -f "$RDIR/state/done" ] && say "finished: $(cat "$RDIR/state/done")"
  local lf; lf=$(st jget "$STATE" last_failure); [ -n "$lf" ] && say "last failure: $lf"
  say "phase: $(st jget "$STATE" phase)"
  st progress --config "$CFG" --state "$STATE" --out "$RUN_HOME/tmp/status-progress.json" | sed 's/^progress: /progress: /'
  say "verdict label: $(st jget "$CFG" verdict_label)"
  say "last log lines ($RUN_HOME/logs/run.log):"; tail -n 5 "$RUN_HOME/logs/run.log" 2>/dev/null | sed 's/^/  /'
  if [ $final = 1 ]; then
    if [ -f "$TREE/results/full/manifest.json" ]; then
      say "manifest: complete=$(st jget "$TREE/results/full/manifest.json" complete), not_verified=$(st jget "$TREE/results/full/manifest.json" not_verified), sets_equal_no_duplicates=$(st jget "$TREE/results/full/manifest.json" sets_equal_no_duplicates), terminal_unidentified=$(st jget "$TREE/results/full/manifest.json" terminal_unidentified)"
      say "  (from $TREE/results/full/manifest.json; $(st jget "$CFG" verdict_label))"
    else say "manifest: not written yet (the run has not reached the manifest phase)"; fi
  fi
}
cmd_stop() {
  local reason="stopped by hand (run.sh stop)" pid pg
  while [ $# -gt 0 ]; do case $1 in --reason) reason=$2; shift ;; *) die "stop: unknown option $1" ;; esac; shift; done
  load_current || die "no current run"
  run mkdir -p "$RDIR/state"
  [ "$DRY" = 0 ] && printf '%s %s\n' "$(utc)" "$reason" > "$(stopped_flag)"
  case $KC_SURVIVAL in
    systemd) run "$SYSTEMCTL" --user stop "$UNIT"; run "$SYSTEMCTL" --user disable "$UNIT" ;;
  esac
  if [ -f "$RDIR/state/pid" ]; then
    read -r pid pg _ < "$RDIR/state/pid"
    if group_live "$pg"; then
      run kill -TERM -- "-$pg"
      for _ in $(seq 1 30); do [ "$DRY" = 1 ] && break; group_live "$pg" || break; sleep 1; done
      if [ "$DRY" = 0 ] && group_live "$pg"; then run kill -KILL -- "-$pg"; fi
    fi
  fi
  [ "$KC_SURVIVAL" = tmux ] && { tmux -L "$TMUX_SOCK" has-session -t "$TMUX_SESSION" 2>/dev/null && run tmux -L "$TMUX_SOCK" kill-session -t "$TMUX_SESSION"; }
  if [ "$DRY" = 0 ]; then
    st jappend "$(REC)" stops "{\"utc\": \"$(utc)\", \"reason\": $(printf '%s' "$reason" | "$HPY" -c 'import json,sys; print(json.dumps(sys.stdin.read()))')}"
    lock_free && commit_results "run $RUN_ID: stopped ($reason)" results/full/run-record.json
  fi
  say "run $RUN_ID stopped: $reason.  bash run.sh start  resumes it (finished shards are kept)."
}

usage() { sed -n '2,/^set -uo/p' "$SELF" | sed '$d' | sed 's/^# \{0,1\}//'; }

# ------------------------------------------------------------------------------------------------- main
main() {
  local -a args=(); local a sub
  for a in "$@"; do [ "$a" = --dry-run ] && DRY=1 || args+=("$a"); done
  set -- "${args[@]+"${args[@]}"}"
  sub=${1:-help}; [ $# -gt 0 ] && shift
  case $sub in
    setup) cmd_setup "$@" ;;
    preflight) SRC_DIR=$RUN_HOME/repo; cmd_preflight "$@" ;;
    selftest) cmd_selftest "$@" ;;
    start) cmd_start "$@" ;;
    status) cmd_status "$@" ;;
    stop) cmd_stop "$@" ;;
    check-tokens) cmd_check_tokens ;;
    zenodo-sandbox-test) cmd_zenodo_sandbox_test ;;
    _run) [ "$DRY" = 1 ] && { say "[dry-run] _run is the long-running process; it is not dry-run"; exit 0; }; cmd_run_inner ;;
    help|-h|--help) usage ;;
    *) usage >&2; exit 2 ;;
  esac
}
main "$@"
