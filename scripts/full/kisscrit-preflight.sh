#!/usr/bin/env bash
# kisscrit-preflight.sh: check whether this Linux machine can run the kisscrit `full` completeness run
# without admin rights, and draft a request to the sysadmin for anything missing.
#
# 2026-10-21, version 2 (revises version 1 of 2026-10-16; control P-preflight-v2; untrusted tooling, no part of
# the proof).  Not yet run on the target server (Ubuntu 26.04).  Read-only by default: it writes only a temporary
# file or two in $HOME (a write test and an exec test, removed at once) and the report.  --try-install also
# installs a test toolchain in ~/kisscrit-preflight/ (about 3 GB; delete it afterwards with:
# rm -rf ~/kisscrit-preflight) and builds and tests the repo's C kernel there.  `run.sh setup` reuses these checks.
#
# Usage:   bash kisscrit-preflight.sh                 (quick checks, about 1 minute)
#          bash kisscrit-preflight.sh --try-install --repo DIR
#                                                     (also installs Miniforge + FLINT, builds csrc/ with that
#                                                      compiler and the repo's flags, runs make test, scans for
#                                                      fused multiply-add, runs the C kernel's pytest files)
#          Options: --need-disk-gb N, --results-dir DIR, --repo DIR (see --help).
# Output:  printed, and saved to ~/kisscrit-preflight-report.txt.  It contains your user name and the host
#          name, nothing secret.  Never needs sudo.

set -u
usage() {
  cat <<'EOF'
Usage: bash kisscrit-preflight.sh [--try-install] [--need-disk-gb N] [--results-dir DIR] [--repo DIR]
  --try-install       also install Miniforge + FLINT in ~/kisscrit-preflight (about 2-3 GB), then build the repo's
                      csrc/ with that compiler and the repo's flags, run its C tests, look for fused multiply-add
                      instructions in the built library, and run the C kernel's pytest files
  --repo DIR          a copy of the kisscrit repo (csrc/, kisscrit/, tests/, scripts/completeness/) for --try-install; default: the repo this
                      script sits in (scripts/full/), if it is run from one
  --need-disk-gb N    free space the run needs, in GB (default 50: an estimate, to be revised after milestone F5/F6)
  --results-dir DIR   where the run will keep its files (default ~/kisscrit-run); its filesystem is checked
EOF
}
TRY_INSTALL=0; NEED_DISK_GB=50; RESDIR="$HOME/kisscrit-run"; REPO=""
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
    --repo)
      [ $# -ge 2 ] || { echo "--repo needs a directory" >&2; exit 2; }
      REPO=$2; shift ;;
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
SLICE="/sys/fs/cgroup/user.slice/user-$(id -u).slice"
[ -r "$SLICE/cpu.max" ] && { read -r Q P < "$SLICE/cpu.max"; [ "$Q" != "max" ] && { warn "user slice CPU limit: about $((Q/P)) cores"; ask "cpu"; }; }
[ -r "$SLICE/memory.max" ] && { M=$(cat "$SLICE/memory.max"); [ "$M" != "max" ] && warn "user slice memory limit: $((M/1024/1024/1024)) GB"; }
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

# ---------------------------------------------------------------------------------------------- build helpers
# fma_check LIB LABEL [ENVDIR]: look for fused multiply-add instructions in a built library; the repo's flag
# -ffp-contract=off must leave none.  x86-64: vfmadd*, vfmsub*, vfnmadd*, vfnmsub*.  aarch64: fmadd, fmsub, fnmadd,
# fnmsub, and also the vector forms fmla and fmls (a superset of the four named in the task).  Any other
# architecture, a missing objdump, or no disassembly: "not checked", as a warning.
fma_check() {
  local lib=$1 label=$2 e=${3:-} arch od pat mn total n
  arch=$(uname -m); od=$(command -v objdump 2>/dev/null || true)
  [ -z "$od" ] && [ -n "$e" ] && od=$(find "$e/bin" -name '*objdump' -type f 2>/dev/null | head -n 1)
  case $arch in
    x86_64|amd64)  pat='^vfn?m(add|sub)' ;;
    aarch64|arm64) pat='^(fmadd|fmsub|fnmadd|fnmsub|fmla|fmls)' ;;
    *) warn "$label: fused multiply-add check not checked (architecture $arch is not covered)"; return 0 ;;
  esac
  [ -n "$od" ] || { warn "$label: fused multiply-add check not checked (objdump missing; architecture $arch)"; return 0; }
  mn=$("$od" -d "$lib" 2>/dev/null | awk -F'\t' 'NF>=3{split($3,a," "); print a[1]}')
  total=$(printf '%s\n' "$mn" | grep -c .)
  [ "$total" -gt 0 ] || { warn "$label: fused multiply-add check not checked (objdump printed no instructions for $lib)"; return 0; }
  n=$(printf '%s\n' "$mn" | grep -Ec "$pat")
  if [ "$n" -eq 0 ]; then ok "$label: no fused multiply-add instruction in $(basename "$lib") (arch $arch, $total instructions scanned)"
  else bad "$label: $n fused multiply-add instructions in $(basename "$lib") (arch $arch); -ffp-contract=off was not honoured"; fi
}

# build_and_test LABEL CC FLINT_PREFIX PYTHON SRC [ENVDIR]: build SRC/csrc with CC and the repo's own flags, run
# make test, look for fused multiply-add instructions, then run the four bit-for-bit pytest files with PYTHON
# (empty: skip them).  FLINT_PREFIX is empty when FLINT is on the compiler's default path.  The result of the
# pytest step is "ok" only when every test passed: a skipped test is a warning, never an ok.
build_and_test() {
  local label=$1 cc=$2 fp=$3 py=$4 src=$5 e=${6:-} blog="$5/build.log" plog="$5/pytest.log" t0 rc w np nf ne ns
  local -a mk ev
  mk=(env -u CFLAGS -u CPPFLAGS -u LDFLAGS -u LIBS "${MAKEBIN:-make}" -C "$src/csrc" "CC=$cc")
  [ -n "$fp" ] && mk+=("FLINT_PREFIX=$fp")
  : > "$blog"
  if "${mk[@]}" all >>"$blog" 2>&1; then ok "$label: csrc built with $("$cc" --version 2>&1 | head -n 1)"
  else bad "$label: csrc build failed (see $blog)"; tail -n 8 "$blog" | sed 's/^/    /'; return 1; fi
  if grep -q -- ' -O2 ' "$blog" && grep -q -- '-fno-fast-math' "$blog" && grep -q -- '-ffp-contract=off' "$blog"
  then ok "$label: build used the repo flags (-O2 -fno-fast-math -ffp-contract=off)"
  else bad "$label: the build log does not show -O2 -fno-fast-math -ffp-contract=off (see $blog)"; fi
  w=$(grep -c 'warning:' "$blog"); [ "$w" -eq 0 ] || warn "$label: the C build printed $w warnings (see $blog)"
  if "${mk[@]}" test >>"$blog" 2>&1 && grep -q 'all containment tests passed' "$blog"
  then ok "$label: make -C csrc test: all containment tests passed"
  else bad "$label: make -C csrc test failed (see $blog)"; tail -n 8 "$blog" | sed 's/^/    /'; fi
  fma_check "$src/build/libkisscheck.so" "$label" "$e"
  if [ -z "$py" ]; then warn "$label: the pytest files were not run (no Python with pytest)"; return 0; fi
  # test_cmodel, test_creplay and test_ckernel read the cond certificates; without them they skip or fail for a
  # reason that has nothing to do with this machine, so say that instead of running them
  if ! compgen -G "$src/results/completeness/certificates/certificates-cond-*-of-40.jsonl" >/dev/null; then
    warn "$label: the pytest files were not run: no results/completeness/certificates/certificates-cond-*-of-40.jsonl in the repo copy (see docs/setup-guide.md A4)"
    return 0
  fi
  # tests that compile their own copies of the C files look for cc/gcc, FLINT's headers and FLINT's library
  mkdir -p "$TMPD/bin"
  printf '#!/bin/sh\nexec "%s" "$@"\n' "$cc" > "$TMPD/bin/cc"; cp "$TMPD/bin/cc" "$TMPD/bin/gcc"; chmod +x "$TMPD/bin/cc" "$TMPD/bin/gcc"
  ev=("PATH=$TMPD/bin:$(dirname "$py"):$PATH" "CC=$cc")
  [ -n "$fp" ] && ev+=("CPATH=$fp/include${CPATH:+:$CPATH}" "LIBRARY_PATH=$fp/lib${LIBRARY_PATH:+:$LIBRARY_PATH}" "LD_LIBRARY_PATH=$fp/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}")
  echo "  running tests/test_carith.py test_cmodel.py test_creplay.py test_ckernel.py (several minutes) ..."
  t0=$(date +%s)
  ( cd "$src" && env "${ev[@]}" "$py" -m pytest -q -p no:cacheprovider -rsfE \
      tests/test_carith.py tests/test_cmodel.py tests/test_creplay.py tests/test_ckernel.py ) >"$plog" 2>&1; rc=$?
  np=$(grep -Eo '[0-9]+ passed' "$plog" | tail -n 1 | cut -d' ' -f1);  np=${np:-0}
  nf=$(grep -Eo '[0-9]+ failed' "$plog" | tail -n 1 | cut -d' ' -f1);  nf=${nf:-0}
  ne=$(grep -Eo '[0-9]+ errors?' "$plog" | tail -n 1 | cut -d' ' -f1); ne=${ne:-0}
  ns=$(grep -Eo '[0-9]+ skipped' "$plog" | tail -n 1 | cut -d' ' -f1); ns=${ns:-0}
  echo "  pytest: exit $rc, $np passed, $nf failed, $ne errors, $ns skipped, $(( $(date +%s) - t0 )) s (log: $plog)"
  if [ "$rc" -eq 0 ] && [ "$np" -gt 0 ] && [ "$ns" -eq 0 ]; then ok "$label: the four pytest files: all $np tests passed"
  elif [ "$rc" -eq 0 ] && [ "$np" -gt 0 ]; then
    warn "$label: the four pytest files: $np passed but $ns skipped, so not all pass (log: $plog)"
    grep '^SKIPPED' "$plog" | sed -E 's/^SKIPPED \[[0-9]+\] //' | sort | uniq -c | head -n 6 | sed 's/^/    skipped: /'
  else
    bad "$label: the four pytest files did not pass (exit $rc; $np passed, $nf failed, $ne errors, $ns skipped; log: $plog)"
    grep -E '^(FAILED|ERROR) ' "$plog" | sort -u | head -n 8 | sed 's/^/    /'
    tail -n 2 "$plog" | sed 's/^/    /'
  fi
}

# ---------------------------------------------------------------------------------------------- optional install test
if [ $TRY_INSTALL -eq 1 ]; then
  hdr "Install test (Miniforge + FLINT + C compiler, in ~/kisscrit-preflight)"
  P="$HOME/kisscrit-preflight"; mkdir -p "$P"; T0=$(date +%s)
  TMPD=$(mktemp -d "${TMPDIR:-/tmp}/kisscrit-preflight.XXXXXX") || TMPD="$P/tmp"
  mkdir -p "$TMPD"; trap 'rm -rf "$TMPD"' EXIT
  if [ ! -x "$P/mf/bin/conda" ]; then
    curl -fsSL -o "$P/mf.sh" "$MF_URL" \
      && bash "$P/mf.sh" -b -p "$P/mf" >"$P/install.log" 2>&1 && ok "Miniforge installed" || bad "Miniforge install failed (see $P/install.log)"
  fi
  if [ -x "$P/mf/bin/conda" ]; then
    # pytest, networkx, mpmath and sympy are what the repo's own tests import besides numpy, scipy, numba and python-flint
    "$P/mf/bin/conda" create -y -q -p "$P/env" -c conda-forge python=3.12 numpy scipy numba python-flint libflint c-compiler make git gh \
      pytest networkx mpmath sympy >>"$P/install.log" 2>&1 \
      && ok "environment created (python-flint, libflint, compiler, numba, git, gh, pytest)" || bad "conda environment failed (see $P/install.log)"
  fi
  E="$P/env"
  if [ -f "$E/include/flint/arb.h" ]; then ok "FLINT headers present (libflint ships arb.h)"; else bad "FLINT headers missing from libflint"; fi
  CC=$(ls "$E"/bin/*-gcc 2>/dev/null | head -1); [ -z "$CC" ] && CC="$E/bin/cc"
  cat > "$P/t.c" <<'CEOF'
#include <stdio.h>
#include <fenv.h>
#include <flint/arb.h>
int main(void){ arb_t x,s; arb_init(x); arb_init(s); arb_set_d(x,0.731); arb_sin(s,x,64);
  printf("  FLINT header %s, library %s\n", FLINT_VERSION, flint_version);
  printf("  arb_sin(0.731) = "); arb_printd(s,20); printf("\n");
  return fegetround()==FE_TONEAREST ? 0 : 1; }
CEOF
  FLINT_RUN=0
  if [ -x "$CC" ] && "$CC" -O2 -fno-fast-math -ffp-contract=off -I"$E/include" -L"$E/lib" -Wl,-rpath,"$E/lib" \
       -o "$P/t" "$P/t.c" -lflint -lm >>"$P/install.log" 2>&1 && "$P/t" >"$P/t.out" 2>&1; then
    ok "compiled and ran a FLINT/Arb test program with the conda compiler (no system gcc needed)"; FLINT_RUN=1; sed -n 2p "$P/t.out"
  else bad "could not compile or run the FLINT test (see $P/install.log)"; fi

  # --- the real thing: the repo's csrc/ built with this compiler and the repo's flags, and its tests
  echo "  build record (this machine, for the report):"
  echo "    machine $(uname -m);  OS ${PRETTY_NAME:-unknown};  $(getconf GNU_LIBC_VERSION 2>/dev/null || echo 'glibc ?');  cores $NCPU;  memory ${MEM_GB} GB"
  [ -x "$CC" ] && echo "    compiler (the environment's): $("$CC" --version 2>&1 | head -n 1)"
  SYSCC=$(command -v cc 2>/dev/null || command -v gcc 2>/dev/null || true)
  [ -n "$SYSCC" ] && echo "    system compiler (not used): $SYSCC, $("$SYSCC" --version 2>&1 | head -n 1)"
  [ "$FLINT_RUN" -eq 1 ] && sed -n 1p "$P/t.out"
  if [ -x "$E/bin/python" ]; then
    "$E/bin/python" -I - <<'PY' 2>&1
import sys
v = [('Python', sys.version.split()[0])]
for mod, name in (('numpy', 'numpy'), ('numba', 'numba'), ('flint', 'python-flint'), ('scipy', 'scipy'), ('pytest', 'pytest')):
    try: v.append((name, __import__(mod).__version__))
    except Exception as ex: v.append((name, 'import failed (%s)' % type(ex).__name__))
print('    ' + ';  '.join('%s %s' % kv for kv in v))
PY
  fi
  if [ -z "$REPO" ]; then
    d=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)
    [ -f "$d/csrc/Makefile" ] && [ -f "$d/kisscrit/__init__.py" ] && REPO=$d
  fi
  if [ -z "$REPO" ] || [ ! -f "$REPO/csrc/Makefile" ] || [ ! -f "$REPO/kisscrit/__init__.py" ] || [ ! -d "$REPO/tests" ] \
       || [ ! -f "$REPO/scripts/completeness/diff_kernel.py" ]; then
    warn "C build and tests not run: no kisscrit repo copy found (give --repo DIR with csrc/, kisscrit/, tests/, scripts/completeness/; see docs/setup-guide.md A4)"
  elif [ ! -x "$CC" ] || [ ! -f "$E/include/flint/arb.h" ]; then
    warn "C build and tests not run: the environment has no compiler or no FLINT headers (see above)"
  else
    SRC="$P/src"; rm -rf "$SRC"; mkdir -p "$SRC"
    echo "  repo copy: $REPO -> $SRC (csrc, kisscrit, tests, scripts/completeness, pyproject.toml; the certificates are linked, not copied)"
    tar --exclude=__pycache__ -C "$REPO" -cf - csrc kisscrit tests scripts/completeness pyproject.toml | tar -C "$SRC" -xf -
    if [ -d "$REPO/results/completeness/certificates" ]; then
      mkdir -p "$SRC/results/completeness"; ln -s "$REPO/results/completeness/certificates" "$SRC/results/completeness/certificates"
    fi
    MAKEBIN=$(command -v make 2>/dev/null || true); [ -z "$MAKEBIN" ] && MAKEBIN="$E/bin/make"
    PYT=""; [ -x "$E/bin/python" ] && "$E/bin/python" -c 'import pytest' 2>/dev/null && PYT="$E/bin/python"
    build_and_test "conda compiler" "$CC" "$E" "$PYT" "$SRC" "$E"
  fi
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
    pkg:curl) add curl "- Install curl:  sudo apt install curl" ;;
    # cron needs no key of its own: a missing cron, a denied crontab and a stopped cron daemon matter only when
    # lingering is not available, and the linger-required text below already lists all three as its alternative.
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
