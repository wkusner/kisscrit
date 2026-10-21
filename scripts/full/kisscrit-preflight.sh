#!/usr/bin/env bash
# kisscrit-preflight.sh: check whether this Linux machine can run the kisscrit `full` completeness run
# without admin rights, and draft a request to the sysadmin for anything missing.
#
# 2026-10-16, version 1 (repo copy of the Project doc claude/kisscrit-preflight.sh; untrusted tooling, no part
# of the proof).  Read-only except --try-install, which installs a test toolchain in ~/kisscrit-preflight/
# (about 2-3 GB; delete it afterwards with: rm -rf ~/kisscrit-preflight).  `run.sh setup` reuses these checks.
#
# Usage:   bash kisscrit-preflight.sh                 (quick checks, about 1 minute)
#          bash kisscrit-preflight.sh --try-install   (also installs Miniforge + FLINT and compiles a test)
# Output:  printed, and saved to ~/kisscrit-preflight-report.txt.  It contains your user name and the host
#          name, nothing secret.  Never needs sudo.

set -u
usage() {
  cat <<'EOF'
Usage: bash kisscrit-preflight.sh [--try-install] [--need-disk-gb N] [--results-dir DIR]
  --try-install       also install Miniforge + FLINT in ~/kisscrit-preflight (about 2-3 GB)
  --need-disk-gb N    free space the run needs, in GB (default 50: an estimate, to be revised after milestone F5/F6)
  --results-dir DIR   where the run will keep its files (default ~/kisscrit-run); its filesystem is checked
EOF
}
TRY_INSTALL=0; NEED_DISK_GB=50; RESDIR="$HOME/kisscrit-run"
while [ $# -gt 0 ]; do
  case $1 in
    --try-install) TRY_INSTALL=1 ;;
    --need-disk-gb)
      [ $# -ge 2 ] || { echo "--need-disk-gb needs a whole number of GB" >&2; exit 2; }
      case $2 in ''|*[!0-9]*) echo "--need-disk-gb: not a whole number: $2" >&2; exit 2 ;; esac
      NEED_DISK_GB=$2; shift ;;
    --results-dir)
      [ $# -ge 2 ] || { echo "--results-dir needs a directory" >&2; exit 2; }
      RESDIR=$2; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done
REPORT="$HOME/kisscrit-preflight-report.txt"
exec > >(tee "$REPORT") 2>&1

PASS=(); WARN=(); FAIL=(); ASK=()
ok()   { printf '  [ ok ] %s\n' "$*"; PASS+=("$*"); }
warn() { printf '  [warn] %s\n' "$*"; WARN+=("$*"); }
bad()  { printf '  [FAIL] %s\n' "$*"; FAIL+=("$*"); }
ask()  { ASK+=("$*"); }
hdr()  { printf '\n== %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }

echo "kisscrit preflight, $(date -u +%Y-%m-%dT%H:%MZ), user $(id -un), host $(hostname -f 2>/dev/null || hostname)"

# ---------------------------------------------------------------------------------------------- system
hdr "System"
. /etc/os-release 2>/dev/null && echo "  OS: ${PRETTY_NAME:-unknown};  kernel $(uname -r);  arch $(uname -m);  $(getconf GNU_LIBC_VERSION 2>/dev/null || echo 'glibc ?')"
NCPU=$(nproc 2>/dev/null || echo 0); echo "  cores (nproc): $NCPU"
# a cgroup CPU limit on this user would cap the run below nproc
CG=$(awk -F: '$1=="0"{print $3}' /proc/self/cgroup 2>/dev/null)
if [ -n "$CG" ] && [ -r "/sys/fs/cgroup$CG/cpu.max" ]; then
  read -r Q P < "/sys/fs/cgroup$CG/cpu.max"
  if [ "$Q" = "max" ]; then ok "no cgroup CPU limit on this session"
  else LIM=$(( Q / P )); warn "cgroup CPU limit: about $LIM cores (of $NCPU)"; ask "cpu"; fi
fi
for d in "/sys/fs/cgroup/user.slice/user-$(id -u).slice"; do
  [ -r "$d/cpu.max" ] && { read -r Q P < "$d/cpu.max"; [ "$Q" != "max" ] && { warn "user slice CPU limit: about $((Q/P)) cores"; ask "cpu"; }; }
  [ -r "$d/memory.max" ] && { M=$(cat "$d/memory.max"); [ "$M" != "max" ] && warn "user slice memory limit: $((M/1024/1024/1024)) GB"; }
done
MEM_GB=$(awk '/MemTotal/{printf "%d", $2/1048576}' /proc/meminfo); AVAIL_GB=$(awk '/MemAvailable/{printf "%d", $2/1048576}' /proc/meminfo)
echo "  memory: ${MEM_GB} GB total, ${AVAIL_GB} GB available"
[ "$MEM_GB" -ge $(( NCPU / 2 )) ] && ok "memory at least 0.5 GB per core" || warn "less than 0.5 GB of memory per core; the run may use fewer workers than cores"

# ---------------------------------------------------------------------------------------------- disk
hdr "Disk"
# free space per filesystem; the results directory may not exist yet, so its nearest existing parent is measured
# (nothing is created)
nearest() { local p=$1; while [ ! -e "$p" ] && [ "$p" != "/" ] && [ "$p" != "." ]; do p=$(dirname "$p"); done; printf '%s' "$p"; }
dfline()  { df -Pk "$(nearest "$1")" 2>/dev/null | awk 'NR==2'; }
free_gb() { dfline "$1" | awk '{printf "%d", $4/1048576}'; }
fs_desc() { dfline "$1" | awk '{m=$6; for(i=7;i<=NF;i++) m=m" "$i; printf "mount %s (%s)", m, $1}'; }
HOME_FREE=$(free_gb "$HOME"); TMP_FREE=$(free_gb /tmp); RES_FREE=$(free_gb "$RESDIR")
HOME_FREE=${HOME_FREE:-0}; TMP_FREE=${TMP_FREE:-0}; RES_FREE=${RES_FREE:-0}
echo "  free space; the run needs ${NEED_DISK_GB} GB (an estimate, to be revised after milestone F5/F6):"
echo "    home:        ${HOME_FREE} GB free, $HOME, $(fs_desc "$HOME")"
echo "    /tmp:        ${TMP_FREE} GB free, $(fs_desc /tmp), type $(findmnt -n -o FSTYPE /tmp 2>/dev/null || echo '?')  (the run keeps its temporary files under ~/kisscrit-run, not /tmp)"
echo "    results dir: ${RES_FREE} GB free, $RESDIR, $(fs_desc "$RESDIR")"
if [ "$RES_FREE" -ge "$NEED_DISK_GB" ]; then ok "at least $NEED_DISK_GB GB free where the results go ($RESDIR)"
else bad "only ${RES_FREE} GB free where the results go ($RESDIR; want $NEED_DISK_GB+)"; ask "disk"; fi
if have quota; then Q=$(quota -s 2>/dev/null | tail -n +3); [ -n "$Q" ] && echo "  quota:" && echo "$Q" | sed 's/^/    /'; fi
TF="$HOME/.kisscrit-preflight-write-test"; if echo x > "$TF" 2>/dev/null; then rm -f "$TF"; ok "home is writable"; else bad "cannot write to home"; fi
# exec from home (noexec mounts would break the toolchain)
TS="$HOME/.kisscrit-preflight-exec-test.sh"; printf '#!/bin/sh\nexit 0\n' > "$TS"; chmod +x "$TS"
if "$TS" 2>/dev/null; then ok "programs can run from home"; else bad "home is mounted noexec"; ask "noexec"; fi; rm -f "$TS"

# ---------------------------------------------------------------------------------------------- limits
hdr "Limits"
echo "  max user processes: $(ulimit -u);  open files: $(ulimit -n);  virtual memory: $(ulimit -v)"
[ "$(ulimit -u)" = "unlimited" ] || [ "$(ulimit -u)" -ge $(( NCPU * 8 )) ] && ok "process limit is enough" || warn "process limit is low for $NCPU workers"

# ---------------------------------------------------------------------------------------------- commands
hdr "Commands (the run installs git, gh, Python, the C compiler and FLINT itself; these must come from the system)"
for c in bash tar gzip; do have "$c" && ok "$c" || bad "$c missing"; done
if have curl; then ok "curl"; elif have wget; then ok "wget (curl missing; wget will do)"; else bad "neither curl nor wget"; ask "pkg:curl"; fi
have tmux && ok "tmux" || { warn "tmux missing"; ask "pkg:tmux"; }
for c in git gcc make python3; do have "$c" && echo "  (optional, present) $c" || echo "  (optional, absent)  $c"; done

# ---------------------------------------------------------------------------------------------- keeping the run alive
hdr "Running after logout and after a reboot"
LINGER=$(loginctl show-user "$(id -un)" -p Linger --value 2>/dev/null || echo unknown)
echo "  linger: $LINGER"
KUP=$( (cat /etc/systemd/logind.conf /etc/systemd/logind.conf.d/*.conf /usr/lib/systemd/logind.conf.d/*.conf 2>/dev/null) \
       | grep -E '^\s*KillUserProcesses' | tail -1 | cut -d= -f2 | tr -d ' ')
echo "  KillUserProcesses: ${KUP:-not set (Ubuntu default: no)}"
USERMGR=no; systemctl --user show-environment >/dev/null 2>&1 && USERMGR=yes
echo "  systemd user manager reachable: $USERMGR"
CRON_OK=no
if have crontab; then
  OUT=$(crontab -l 2>&1); RC=$?
  if [ $RC -eq 0 ] || echo "$OUT" | grep -qi "no crontab"; then CRON_OK=yes; fi
  echo "$OUT" | grep -qi "not allowed" && CRON_OK=denied
fi
CRON_RUN=no; pgrep -x cron >/dev/null 2>&1 && CRON_RUN=yes
echo "  crontab usable: $CRON_OK;  cron daemon running: $CRON_RUN"
if [ "$LINGER" = "yes" ] && [ "$USERMGR" = "yes" ]; then
  ok "lingering is on: the run can be a systemd user service (best option)"
elif [ "$CRON_OK" = "yes" ] && [ "$CRON_RUN" = "yes" ] && have tmux && [ "${KUP:-no}" != "yes" ]; then
  ok "tmux + cron @reboot will keep the run alive (lingering would be better)"; ask "linger-optional"
else
  bad "no way yet to keep the run alive after logout and reboots"; ask "linger-required"
fi
# a direct test: does a detached process survive this check?  (only meaningful if you log out and back in)
echo "  tip: to test survival by hand: run 'nohup sleep 600 &', log out, log back in, and run 'pgrep -u $(id -un) sleep'."

# ---------------------------------------------------------------------------------------------- network
hdr "Outbound HTTPS"
for v in https_proxy HTTPS_PROXY http_proxy no_proxy; do [ -n "${!v:-}" ] && echo "  proxy setting $v=${!v}"; done
# probe URL WHAT [range]: 2xx and 3xx are ok; 4xx and 5xx are a warning with the code (the host answered, but not
# with success: a proxy or firewall often answers 403 or 407 for a blocked host); no response at all is a FAIL and
# goes into the sysadmin request.  "range" asks for the first byte only (Range: bytes=0-0), so a download URL is
# not downloaded.  Only the host name of the final URL is printed: the final URL of a release download is a signed
# link, and it must not reach the report.
probe() {
  local url=$1 what=$2 mode=${3:-get} out last code err="" rc=0 fhost="" uhost="" via="" hint=""
  if have curl; then
    local -a o=(-sS -o /dev/null -w '\n%{http_code} %{url_effective}' --max-time 20 -L)
    [ "$mode" = range ] && o+=(-r 0-0)
    out=$(curl "${o[@]}" "$url" 2>&1); rc=$?
    last=${out##*$'\n'}; err=${out%$'\n'*}; err=${err//$'\n'/ }
    code=${last%% *}; fhost=${last#* }; fhost=${fhost#*://}; fhost=${fhost%%[/?#]*}
  else
    local -a o=(-q -S -T 20 -O /dev/null)
    [ "$mode" = range ] && o+=(--header='Range: bytes=0-0')
    out=$(wget "${o[@]}" "$url" 2>&1); rc=$?
    code=$(printf '%s\n' "$out" | awk '/^ *HTTP\//{c=$2} END{print c}')
    fhost=$(printf '%s\n' "$out" | awk '/^ *[Ll]ocation:/{l=$2} END{print l}'); fhost=${fhost#*://}; fhost=${fhost%%[/?#]*}
    err=$(printf '%s\n' "$out" | tail -n 1)
  fi
  uhost=${url#*://}; uhost=${uhost%%[/?#]*}
  [ -n "$fhost" ] && [ "$fhost" != "$uhost" ] && via=", final host $fhost"
  case $code in
    2??|3??) ok "$what ($url): HTTP $code$via" ;;
    4??|5??)
      case $code in 403|407) hint=" (a proxy or firewall may be blocking this host)" ;; esac
      warn "$what ($url): HTTP $code$via$hint" ;;
    *) bad "$what ($url): no response (exit $rc${err:+: ${err:0:120}})"; ask "net:$url" ;;
  esac
}
MF_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
# github.com itself is probed through the git smart-HTTP endpoint of a public repository, the path that git clone,
# git fetch and git push use: a bare "https://github.com/" returned 400 in the cloud container (a response, so it
# was counted as reachable), and a repository page returned 403 there; this endpoint returned 200.  It needs no
# login, no git and no private repository.
probe "https://github.com/conda-forge/miniforge/info/refs?service=git-upload-pack" "GitHub over HTTPS (the git protocol the run uses)"
probe "https://api.github.com/"                              "GitHub API (gh login)"
probe "$MF_URL"                                              "Miniforge installer (GitHub release; first byte only)" range
probe "https://conda.anaconda.org/conda-forge/noarch/current_repodata.json" "conda-forge packages"
probe "https://zenodo.org/api/records?size=1"                "Zenodo (certificate archive)"
probe "https://sandbox.zenodo.org/api/records?size=1"        "Zenodo sandbox (tests)"
probe "https://pypi.org/simple/"                             "PyPI (fallback for Python packages)"

# ---------------------------------------------------------------------------------------------- optional install test
if [ $TRY_INSTALL -eq 1 ]; then
  hdr "Install test (Miniforge + FLINT + C compiler, in ~/kisscrit-preflight)"
  P="$HOME/kisscrit-preflight"; mkdir -p "$P"; T0=$(date +%s)
  if [ ! -x "$P/mf/bin/conda" ]; then
    curl -fsSL -o "$P/mf.sh" "$MF_URL" \
      && bash "$P/mf.sh" -b -p "$P/mf" >"$P/install.log" 2>&1 && ok "Miniforge installed" || bad "Miniforge install failed (see $P/install.log)"
  fi
  if [ -x "$P/mf/bin/conda" ]; then
    "$P/mf/bin/conda" create -y -q -p "$P/env" -c conda-forge python=3.12 numpy scipy numba python-flint libflint c-compiler make git gh \
      >>"$P/install.log" 2>&1 && ok "environment created (python-flint, libflint, compiler, numba, git, gh)" || bad "conda environment failed (see $P/install.log)"
  fi
  E="$P/env"
  if [ -f "$E/include/flint/arb.h" ]; then ok "FLINT headers present (libflint ships arb.h)"; else bad "FLINT headers missing from libflint"; fi
  if [ -x "$E/bin/python" ]; then
    "$E/bin/python" -c "import flint, numba, scipy; print('  python-flint', flint.__version__, ' numba', numba.__version__)" 2>&1
  fi
  CC=$(ls "$E"/bin/*-gcc 2>/dev/null | head -1); [ -z "$CC" ] && CC="$E/bin/cc"
  cat > "$P/t.c" <<'EOF'
#include <stdio.h>
#include <fenv.h>
#include <flint/arb.h>
int main(void){ arb_t x,s; arb_init(x); arb_init(s); arb_set_d(x,0.731); arb_sin(s,x,64);
  printf("  arb_sin(0.731) = "); arb_printd(s,20); printf("\n");
  return fegetround()==FE_TONEAREST ? 0 : 1; }
EOF
  if [ -x "$CC" ] && "$CC" -O2 -fno-fast-math -ffp-contract=off -I"$E/include" -L"$E/lib" -Wl,-rpath,"$E/lib" \
       -o "$P/t" "$P/t.c" -lflint -lm >>"$P/install.log" 2>&1 && "$P/t"; then
    ok "compiled and ran a FLINT/Arb test program with the conda compiler (no system gcc needed)"
  else bad "could not compile or run the FLINT test (see $P/install.log)"; fi
  echo "  install test took $(( $(date +%s) - T0 )) s;  size: $(du -sh "$P" 2>/dev/null | cut -f1)  (delete with: rm -rf $P)"
fi

# ---------------------------------------------------------------------------------------------- summary
hdr "Summary"
echo "  ${#PASS[@]} ok, ${#WARN[@]} warnings, ${#FAIL[@]} failures"
for f in "${FAIL[@]+"${FAIL[@]}"}"; do echo "  FAIL: $f"; done
for w in "${WARN[@]+"${WARN[@]}"}"; do echo "  warn: $w"; done

# ---------------------------------------------------------------------------------------------- request to sysadmin
declare -A SEEN; LINES=()
add() { [ -z "${SEEN[$1]:-}" ] && { SEEN[$1]=1; LINES+=("$2"); }; }
NET=()
for a in "${ASK[@]+"${ASK[@]}"}"; do
  case $a in
    net:*) NET+=("${a#net:}") ;;
    pkg:tmux) add tmux "- Install tmux:  sudo apt install tmux" ;;
    pkg:cron) add cron "- Install and enable cron:  sudo apt install cron && sudo systemctl enable --now cron" ;;
    pkg:curl) add curl "- Install curl:  sudo apt install curl" ;;
    cron-allow) add cronallow "- Allow my account to use crontab (it is denied now; /etc/cron.allow or /etc/cron.deny)." ;;
    linger-required)
      ALT=""
      have tmux || ALT="sudo apt install tmux; "
      have crontab || ALT="${ALT}sudo apt install cron && sudo systemctl enable --now cron; "
      [ "$CRON_OK" = "denied" ] && ALT="${ALT}allow my account to use crontab (/etc/cron.allow); "
      [ "$CRON_RUN" = "no" ] && have crontab && ALT="${ALT}sudo systemctl enable --now cron; "
      add linger "- Enable lingering for my account, so a long computation keeps running after I log out and restarts after a reboot:  sudo loginctl enable-linger $(id -un)${ALT:+
    (If you prefer not to enable lingering, this alternative also works: ${ALT%; }.)}" ;;
    linger-optional) add linger "- (Optional but preferred) enable lingering for my account, so the computation can run as a systemd user service:  sudo loginctl enable-linger $(id -un)" ;;
    disk) add disk "- About ${NEED_DISK_GB} GB of disk space I can write to (estimate, to be revised after milestone F5/F6; ${RES_FREE} GB free now where the results would go, ${RESDIR}), e.g. a directory on a data volume." ;;
    cpu) add cpu "- Lift or explain the CPU limit on my account (cgroup): the computation is meant to use all $NCPU cores for several days." ;;
    noexec) add noexec "- My home directory appears to be mounted noexec; I need a location where I can install and run programs." ;;
  esac
done
[ ${#NET[@]} -gt 0 ] && LINES=("- Allow outbound HTTPS (port 443) from this machine to: $(printf '%s ' "${NET[@]}" | sed -E 's#https://([^/]+)[^ ]*#\1#g')" "${LINES[@]+"${LINES[@]}"}")

hdr "Draft request to the sysadmin"
if [ ${#LINES[@]} -eq 0 ]; then
  echo "  Nothing needed: this machine can run kisscrit without admin help."
else
  cat <<EOF
  Subject: $(hostname -s): small setup request for a long-running computation

  Hi,

  I'll be running a multi-day research computation on $(hostname -s) (all cores, installed entirely in my home
  directory; no root needed for the software). A preflight check found a few things I can't set up myself:

$(printf '  %s\n' "${LINES[@]}")

  Nothing else is needed. The software (Python, a C compiler, FLINT) installs under my home directory.
  Thanks!
  $(id -un)
EOF
fi
echo; echo "Report saved to $REPORT"
