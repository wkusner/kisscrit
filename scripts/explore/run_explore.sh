#!/usr/bin/env bash
# run_explore.sh: overnight numerical exploration of critical configurations of N points on S^(D-1) on many cores.
# Status: untrusted (exploration tooling; no part of the proof).  Everything it finds is numerical (double
# precision) and unverified; no proof may cite it.
#
# What it does: it runs scripts/explore/explore_worker.py (the existing random_ascent and surgery of
# kisscrit/search.py, unchanged) on CORES workers in batches until the wall-clock budget ends; after each batch it
# merges the worker census files (explore_merge.py); at the end it writes the refine inputs (explore_export.py).
# Docs: docs/explore-run.md.
#
# Subcommands:
#   start   [flags]  launch (or resume) the run in a detached tmux session; --foreground runs it in this shell
#   status           elapsed time, workers alive, batches, types at the last merge, last Euler line
#   stop             end the run cleanly (workers dump their census on SIGTERM), then merge and export
#   merge            merge the worker files found so far into result/merged-census.json
#   export           write refine-input.json, refine-queue.tsv and summary.txt from the merged census
# Flags (each has a default, set in the variables below; a flag overrides its variable):
#   --n N --d D --cores C --nice K --hours H --runs R --rounds K --near M --batch-minutes B --surgery-types T
#   --seed S --run-home DIR --conda-env NAME --tree DIR --foreground
#
# Priority rule: the exploration never preempts the completeness runs.  Workers run at nice 19, and CORES may not
# exceed nproc.  While another server job runs, lower CORES: `stop`, then `start --cores K` resumes the same run
# with K workers (only --cores and --nice may change on a resume).
#
# Everything lives under RUN_HOME/explore (default $HOME/kisscrit-run, as in scripts/full/run.sh; never /tmp).
# result/ holds the only files meant to leave the machine (manifest, merged census, refine inputs, summary); they
# contain no host, user or path strings.  A second `start` with the same RUN_HOME resumes the run.

set -uo pipefail
export LC_ALL=C

# ----------------------------------------------------------------------------------------------- parameters
N=12                    # number of points
D=3                     # ambient dimension; points lie on S^(D-1), so 3 means the sphere S^2
CORES=16                # workers run in parallel; the default is capped at nproc, an explicit --cores above nproc is refused
NICE=19                 # nice level of every worker (19 is the lowest priority)
HOURS=10                # wall-clock budget of the whole run, in hours (decimals allowed)
RUNS=60                 # random ascents per worker and batch
ROUNDS=2                # surgery rounds per surgery call
NEAR=8                  # nearest non-contacts that surgery tries (the value of scripts/ascent_census.py)
BATCH_MINUTES=30        # wall-clock length of one batch; workers dump and stop at its end
SURGERY_TYPES=60        # at most this many not-yet-used types per worker are given to surgery in a batch
SEED=20261027           # base seed; worker j is seeded from (SEED, j)
RUN_HOME="${KISSCRIT_RUN_HOME:-$HOME/kisscrit-run}"   # as in scripts/full/run.sh
CONDA_ENV=kisscrit-run  # conda environment (under KISSCRIT_MINIFORGE, default $HOME/miniforge3) that runs the workers
TREE_ARG=""             # run the code in this directory as it is; empty: a snapshot of the repository's HEAD
FOREGROUND=0            # 1: run the loop in this shell instead of a tmux session
TMUX_SOCK="${KISSCRIT_TMUX_SOCKET:-kisscrit}"         # tmux socket name, as in scripts/full/run.sh
TMUX_SESSION=kisscrit-explore

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
HERE="$(dirname "$SELF")"
REPO="$(cd "$HERE/../.." && pwd)"
declare -A GIVEN=()     # flags given on this command line (a resume refuses changed ones)

die()  { printf 'run_explore.sh: %s\n' "$*" >&2; exit 1; }
say()  { printf '%s\n' "$*"; }
utc()  { date -u '+%Y-%m-%dT%H:%M:%SZ'; }

# ----------------------------------------------------------------------------------------------- arguments
parse_flags() {
  while [ $# -gt 0 ]; do
    local flag=$1
    case $flag in
      --foreground) FOREGROUND=1; shift; continue ;;
      --n|--d|--cores|--nice|--hours|--runs|--rounds|--near|--batch-minutes|--surgery-types|--seed|--run-home|--conda-env|--tree) ;;
      *) die "unknown flag $flag (see the header of this script)" ;;
    esac
    [ $# -ge 2 ] || die "$flag needs a value"
    local v=$2; shift 2
    case $flag in
      --n) N=$v ;; --d) D=$v ;; --cores) CORES=$v ;; --nice) NICE=$v ;; --hours) HOURS=$v ;; --runs) RUNS=$v ;;
      --rounds) ROUNDS=$v ;; --near) NEAR=$v ;; --batch-minutes) BATCH_MINUTES=$v ;; --surgery-types) SURGERY_TYPES=$v ;;
      --seed) SEED=$v ;; --run-home) RUN_HOME=$v ;; --conda-env) CONDA_ENV=$v ;; --tree) TREE_ARG=$v ;;
    esac
    GIVEN[$flag]=$v
  done
}

is_int()  { [[ ${1:-} =~ ^[0-9]+$ ]]; }
is_pos()  { [[ ${1:-} =~ ^[0-9]+(\.[0-9]+)?$ ]] && awk -v x="$1" 'BEGIN{exit !(x > 0)}'; }

# validate: refuse parameters that cannot work, before anything is created
validate() {
  # nproc reads OMP_NUM_THREADS and OMP_THREAD_LIMIT, which a shell or a conda environment may set to 1
  local np; np=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
  is_int "$N" && [ "$N" -ge 3 ] || die "--n must be an integer >= 3 (got '$N')"
  is_int "$D" && [ "$D" -ge 2 ] || die "--d must be an integer >= 2 (got '$D')"
  is_int "$CORES" && [ "$CORES" -ge 1 ] || die "--cores must be an integer >= 1 (got '$CORES')"
  if [ -z "${GIVEN[--cores]:-}" ] && [ "$CORES" -gt "$np" ]; then   # the default is capped, an explicit value is not
    say "note: this machine has $np cores, so the default of $CORES workers becomes $np"; CORES=$np
  fi
  [ "$CORES" -le "$np" ] || die "--cores $CORES is more than nproc ($np): refused, so the exploration cannot starve other jobs"
  is_int "$NICE" && [ "$NICE" -le 19 ] || die "--nice must be an integer from 0 to 19 (got '$NICE')"
  is_pos "$HOURS" || die "--hours must be a positive number (got '$HOURS')"
  is_int "$RUNS" || die "--runs must be an integer >= 0 (got '$RUNS')"
  is_int "$ROUNDS" || die "--rounds must be an integer >= 0 (got '$ROUNDS')"
  is_int "$NEAR" && [ "$NEAR" -ge 1 ] || die "--near must be an integer >= 1 (got '$NEAR')"
  is_pos "$BATCH_MINUTES" || die "--batch-minutes must be a positive number (got '$BATCH_MINUTES')"
  is_int "$SURGERY_TYPES" || die "--surgery-types must be an integer >= 0 (got '$SURGERY_TYPES')"
  is_int "$SEED" || die "--seed must be an integer >= 0 (got '$SEED')"
}

# ----------------------------------------------------------------------------------------------- layout
set_dirs() {
  EXP="$RUN_HOME/explore"; STATE="$EXP/state"; BATCHES="$EXP/batches"; LOGS="$EXP/logs"; RESULT="$EXP/result"
  PARAMS="$STATE/params.env"
}
# the settings file is written by write_params below, so its location is not constant for shellcheck
# shellcheck disable=SC1090
load_params() { [ -f "$PARAMS" ] || die "no exploration run in this run-home (nothing to $1)"; set -a; . "$PARAMS"; set +a; }
lock_free() { ( exec 9>>"$STATE/explore.lock"; flock -n 9 ) 2>/dev/null; }
# alive_workers: how many of the worker processes of the current batch are still running
alive_workers() {
  local n=0 p
  if [ -f "$STATE/pids" ]; then while read -r p; do kill -0 "$p" 2>/dev/null && n=$((n + 1)); done < "$STATE/pids"; fi
  echo "$n"
}

pick_python() {
  local env_python="${KISSCRIT_MINIFORGE:-$HOME/miniforge3}/envs/$CONDA_ENV/bin/python"
  if [ -n "${KISSCRIT_PYTHON:-}" ]; then PY=$KISSCRIT_PYTHON
  elif [ -x "$env_python" ]; then PY=$env_python
  else
    PY=$(command -v python3 || true)
    [ -n "$PY" ] || die "no python found (conda environment '$CONDA_ENV' is missing and python3 is not on PATH)"
    say "note: conda environment '$CONDA_ENV' not found; using $(basename "$PY") from PATH"
  fi
}

# write_params: the settings of the run, read by every later subcommand and by the loop
write_params() {
  {
    for v in N D CORES NICE HOURS RUNS ROUNDS NEAR BATCH_MINUTES SURGERY_TYPES SEED RUN_ID COMMIT TREE PY START_EPOCH DEADLINE_EPOCH; do
      printf '%s=%q\n' "$v" "${!v}"
    done
  } > "$PARAMS.tmp" && mv -f "$PARAMS.tmp" "$PARAMS"
}

# make_tree: pin the code.  Default: a snapshot of the repository's HEAD (the commit is what the manifest records);
# --tree DIR runs DIR as it is.  Uncommitted changes are not in a snapshot.
make_tree() {
  if [ -n "$TREE_ARG" ]; then
    TREE=$(cd "$TREE_ARG" && pwd) || die "--tree $TREE_ARG is not a directory"
    COMMIT=$(git -C "$TREE" rev-parse HEAD 2>/dev/null || echo unknown)
  elif git -C "$REPO" rev-parse HEAD >/dev/null 2>&1; then
    COMMIT=$(git -C "$REPO" rev-parse HEAD); TREE="$EXP/tree"; mkdir -p "$TREE"
    git -C "$REPO" archive HEAD kisscrit scripts/explore | tar -x -C "$TREE" || die "git archive failed"
    git -C "$REPO" archive HEAD data/n12_above1.json 'results/census_n*.json' 2>/dev/null | tar -x -C "$TREE" 2>/dev/null || true
    if [ -n "$(git -C "$REPO" status --porcelain -- kisscrit scripts/explore)" ]; then
      say "note: uncommitted changes in kisscrit/ or scripts/explore/ are not in the snapshot of commit ${COMMIT:0:12}"
    fi
  else
    TREE=$REPO; COMMIT=$(tr -d ' \n' < "$REPO/COMMIT" 2>/dev/null || echo unknown)
  fi
}

write_manifest() {
  "$PY" - "$RESULT/manifest.json" <<'EOF'
import json, os, platform, sys
import numpy, scipy, networkx
e = os.environ
def num(k): return float(e[k]) if '.' in e[k] else int(e[k])
m = dict(tool='scripts/explore/run_explore.sh', run_id=e['RUN_ID'], commit=e['COMMIT'],
         label='numerical (double precision); nothing in this run is verified; no value is a proof input',
         parameters=dict(N=num('N'), d=num('D'), cores=num('CORES'), nice=num('NICE'), hours=num('HOURS'), runs=num('RUNS'),
                         rounds=num('ROUNDS'), near=num('NEAR'), batch_minutes=num('BATCH_MINUTES'),
                         surgery_types=num('SURGERY_TYPES'), seed=num('SEED')),
         started_utc=e['STARTED_UTC'], deadline_utc=e['DEADLINE_UTC'],
         versions=dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
                       networkx=networkx.__version__, system=platform.system()),
         restarts=[])
tmp = sys.argv[1] + '.tmp'
json.dump(m, open(tmp, 'w'), indent=1); os.replace(tmp, sys.argv[1])
EOF
}

# note_restart: record a resume (and any new cores or nice) in the manifest
note_restart() {
  "$PY" - "$RESULT/manifest.json" "$(utc)" "$CORES" "$NICE" <<'EOF'
import json, os, sys
p = sys.argv[1]; m = json.load(open(p))
m['restarts'].append(dict(utc=sys.argv[2], cores=int(sys.argv[3]), nice=int(sys.argv[4])))
json.dump(m, open(p + '.tmp', 'w'), indent=1); os.replace(p + '.tmp', p)
EOF
}

# ----------------------------------------------------------------------------------------------- merge, export
do_merge() {
  "$PY" "$TREE/scripts/explore/explore_merge.py" --out "$RESULT/merged-census.json" --lower-dir "$TREE/results" "$BATCHES" \
        > "$RESULT/merge-report.txt" 2>> "$EXP/loop.log"
}
do_export() {
  local known=()
  [ "$N" = 12 ] && [ -f "$TREE/data/n12_above1.json" ] && known+=(--known "$TREE/data/n12_above1.json")
  [ -f "$TREE/results/census_n$N.json" ] && known+=(--known "$TREE/results/census_n$N.json")
  [ -f "$RESULT/merged-census.json" ] || { echo "nothing to export: no merged census yet (run merge first)" >&2; return 1; }
  "$PY" "$TREE/scripts/explore/explore_export.py" "$RESULT/merged-census.json" --out-dir "$RESULT" ${known[@]+"${known[@]}"}
}

# ----------------------------------------------------------------------------------------------- the loop
log() { printf '%s %s\n' "$(utc)" "$*" >> "$EXP/loop.log"; }

run_batch() {   # run_batch BATCH SECONDS: CORES workers, each with its slice of the unused merged types
  local b=$1 secs=$2 i idx start=() bdir; bdir="$BATCHES/b$(printf %04d "$b")"; mkdir -p "$bdir"
  [ -f "$RESULT/merged-census.json" ] && start=(--start "$RESULT/merged-census.json")
  local next; next=$(cat "$STATE/next_index" 2>/dev/null || echo 0); echo $((next + CORES)) > "$STATE/next_index"
  : > "$STATE/pids"
  for ((i = 0; i < CORES; i++)); do
    idx=$((next + i))
    nice -n "$NICE" "$PY" "$TREE/scripts/explore/explore_worker.py" "$N" "$D" --out "$bdir/w$(printf %06d "$idx").json" \
         --seed "$SEED" --index "$idx" --runs "$RUNS" --rounds "$ROUNDS" --near "$NEAR" --budget "$secs" \
         ${start[@]+"${start[@]}"} --slice-index "$i" --slice-count "$CORES" --surgery-types "$SURGERY_TYPES" \
         --commit "$COMMIT" --run-id "$RUN_ID" --batch "$b" >> "$LOGS/b$(printf %04d "$b")-w$(printf %03d "$i").log" 2>&1 &
    echo $! >> "$STATE/pids"
  done
  wait
}

loop() {
  load_params run; exec 9>>"$STATE/explore.lock"; flock -n 9 || die "another copy of this run is active (see status)"
  rm -f "$STATE/stop" "$STATE/finished"; log "loop started: commit ${COMMIT:0:12}, $CORES workers, nice $NICE"
  local reason=deadline b secs t0 left
  while :; do
    [ -e "$STATE/stop" ] && { reason=stopped; break; }
    left=$((DEADLINE_EPOCH - $(date +%s))); [ "$left" -gt 10 ] || break
    b=$(find "$BATCHES" -maxdepth 1 -name 'b[0-9]*' 2>/dev/null | wc -l)
    secs=$(awk -v m="$BATCH_MINUTES" -v l="$left" 'BEGIN{s = m * 60; if (l < s) s = l; printf "%d", s}')
    log "batch $b: $CORES workers for $secs s"; t0=$(date +%s)
    run_batch "$b" "$secs"; do_merge
    log "batch $b done: $(grep -m1 '^merge:' "$RESULT/merge-report.txt")"
    # workers that end within 10 s of a batch of 30 s or more, with no stop requested, are failing (see the logs)
    if [ ! -e "$STATE/stop" ] && [ "$secs" -ge 30 ] && [ $(($(date +%s) - t0)) -lt 10 ]; then
      log "workers exit at once; see the logs"; reason=failed; break; fi
  done
  do_merge || log "final merge failed: no readable worker census"
  do_export >> "$EXP/loop.log" 2>&1 || log "final export skipped"
  echo "$reason" > "$STATE/finished"; log "loop finished: $reason"
  [ "$reason" != failed ] || exit 3
}

# ----------------------------------------------------------------------------------------------- subcommands
# resume_flags: on a resume the run's own parameters win; a flag that differs from them is refused, except
# --cores and --nice, which may be lowered or raised (docs/explore-run.md, "priority rule")
resume_flags() {
  local f var want
  for f in --n:N --d:D --hours:HOURS --runs:RUNS --rounds:ROUNDS --near:NEAR --batch-minutes:BATCH_MINUTES \
           --surgery-types:SURGERY_TYPES --seed:SEED; do
    want=${GIVEN[${f%%:*}]:-}; var=${f##*:}
    [ -z "$want" ] || [ "$want" = "${!var}" ] || die "this run-home holds a run with ${f%%:*} ${!var}; only --cores and --nice may change on a resume"
  done
  [ -z "${GIVEN[--tree]:-}" ] || die "--tree is fixed when the run starts"
}

cmd_start() {
  parse_flags "$@"; set_dirs
  if [ -f "$PARAMS" ]; then
    local cores_new=${GIVEN[--cores]:-} nice_new=${GIVEN[--nice]:-}
    load_params start; resume_flags
    [ -z "$cores_new" ] || CORES=$cores_new; [ -z "$nice_new" ] || NICE=$nice_new
    validate
    [ "$(cat "$STATE/finished" 2>/dev/null)" != deadline ] || die "this run reached its deadline; use another --run-home for a new run"
    lock_free || die "the run is active (see status)"
    write_params; note_restart; say "resuming run $RUN_ID with $CORES workers"
  else
    validate; pick_python; mkdir -p "$STATE" "$BATCHES" "$LOGS" "$RESULT"
    RUN_ID="explore-n$N-d$D-$(date -u +%Y%m%dT%H%M%SZ)"; make_tree
    START_EPOCH=$(date +%s); DEADLINE_EPOCH=$((START_EPOCH + $(awk -v h="$HOURS" 'BEGIN{printf "%d", h * 3600}')))
    STARTED_UTC=$(date -u -d "@$START_EPOCH" '+%Y-%m-%dT%H:%M:%SZ'); DEADLINE_UTC=$(date -u -d "@$DEADLINE_EPOCH" '+%Y-%m-%dT%H:%M:%SZ')
    export N D CORES NICE HOURS RUNS ROUNDS NEAR BATCH_MINUTES SURGERY_TYPES SEED RUN_ID COMMIT STARTED_UTC DEADLINE_UTC
    write_params; write_manifest; say "started run $RUN_ID at commit ${COMMIT:0:12}"
  fi
  if [ "$FOREGROUND" = 1 ]; then loop; return; fi
  command -v tmux >/dev/null 2>&1 || die "tmux is not installed; use --foreground inside screen or nohup"
  tmux -L "$TMUX_SOCK" has-session -t "$TMUX_SESSION" 2>/dev/null && die "tmux session $TMUX_SESSION already exists"
  tmux -L "$TMUX_SOCK" new-session -d -s "$TMUX_SESSION" "bash '$SELF' _loop --run-home '$RUN_HOME'"
  say "running in the detached tmux session $TMUX_SESSION (socket $TMUX_SOCK); watch with: bash $(basename "$SELF") status"
}

cmd_status() {
  parse_flags "$@"; set_dirs
  [ -f "$PARAMS" ] || { say "no exploration run in this run-home"; return 0; }
  load_params status
  local now elapsed state; now=$(date +%s); elapsed=$((now - START_EPOCH))
  if [ -f "$STATE/finished" ]; then state="finished ($(cat "$STATE/finished"))"
  elif lock_free; then state="not running (resume with start)"
  else state=running; fi
  say "run $RUN_ID (commit ${COMMIT:0:12}): $state"
  local batches; batches=$(find "$BATCHES" -maxdepth 1 -name 'b[0-9]*' 2>/dev/null | wc -l)
  say "elapsed $((elapsed / 3600))h$(((elapsed % 3600) / 60))m of $HOURS h; workers alive $(alive_workers) of $CORES; batches started $batches"
  if [ -f "$RESULT/merge-report.txt" ]; then
    say "at the last merge ($(date -u -r "$RESULT/merge-report.txt" '+%Y-%m-%dT%H:%MZ')):"
    grep -E '^(merge:|Euler check)' "$RESULT/merge-report.txt" | sed 's/^/  /'
  else say "no merge yet (the first batch is still running)"; fi
}

cmd_stop() {
  parse_flags "$@"; set_dirs; load_params stop
  touch "$STATE/stop"
  if [ -f "$STATE/pids" ]; then while read -r p; do kill -TERM "$p" 2>/dev/null || true; done < "$STATE/pids"; fi
  local i; for i in $(seq 1 900); do lock_free && break; sleep 1; done
  lock_free || die "the run did not stop within 15 minutes (see status)"
  tmux -L "$TMUX_SOCK" kill-session -t "$TMUX_SESSION" 2>/dev/null || true
  [ -f "$STATE/finished" ] || { do_merge && do_export; echo stopped > "$STATE/finished"; }
  say "stopped; merged census and refine inputs are in result/"
}

main() {
  local cmd=${1:-help}; shift || true
  case $cmd in
    start) cmd_start "$@" ;;
    status) cmd_status "$@" ;;
    stop) cmd_stop "$@" ;;
    merge) parse_flags "$@"; set_dirs; load_params merge; do_merge; cat "$RESULT/merge-report.txt" ;;
    export) parse_flags "$@"; set_dirs; load_params export; do_export || exit 1 ;;
    _loop) parse_flags "$@"; set_dirs; loop ;;
    help|-h|--help) sed -n '2,/^$/p' "$SELF" | sed 's/^# \{0,1\}//' ;;
    *) die "unknown subcommand '$cmd' (start, status, stop, merge, export, help)" ;;
  esac
}
main "$@"
