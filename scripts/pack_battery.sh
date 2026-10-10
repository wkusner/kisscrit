#!/bin/bash
# usage: bash scripts/pack_battery.sh  (logs to results/pack_battery.log)
cd "$(dirname "$0")/.."
run(){ timeout 3000 python3 scripts/pack_explore.py "$@" 2>&1 | head -40; }
for n in 2 3 4 5; do run sq $n 40 1500; done
for n in 2 3 4; do run tor $n 40 1500; done
for n in 2 3 4; do run hex $n 40 1500; done
for n in 2 3 4 5; do run disk $n 40 1500; done
for n in 2 3 4; do run cube $n 40 1500; done
for n in 2 3 4; do run ball $n 40 1500; done
for n in 6 7; do run sq $n 60 2500; done
for n in 5; do run tor $n 60 2500; run hex $n 60 2500; done
for n in 6 7; do run disk $n 60 2500; done
