#!/bin/bash
cd "$(dirname "$0")/.."
python3 scripts/pack_maxima.py tor3 4 5 30
for n in 2 3; do python3 scripts/pack_explore.py cube $n 40 1500; done
for n in 2 3; do python3 scripts/pack_explore.py tor3 $n 40 1500; done
for n in 2 3; do python3 scripts/pack_explore.py ball $n 40 1500; done
for w in sq tor hex disk; do python3 scripts/pack_explore.py $w 4 60 6000; done
python3 scripts/pack_euler.py
