#!/bin/bash
cd "$(dirname "$0")/.."
while pgrep -f "[p]ack_maxima_all" >/dev/null; do sleep 30; done
for n in 2 3; do python3 scripts/pack_explore.py cube $n 40 1500; done
for n in 2 3; do python3 scripts/pack_explore.py tor3 $n 40 1500; done
for n in 2 3; do python3 scripts/pack_explore.py ball $n 40 1500; done
for w in sq tor hex disk; do python3 scripts/pack_explore.py $w 4 60 6000; done
python3 scripts/pack_euler.py
