#!/bin/bash
cd "$(dirname "$0")/.."
python3 scripts/pack_maxima.py sq 7 10 60
python3 scripts/pack_maxima.py disk 2 10 60
python3 scripts/pack_maxima.py tor 2 8 60
python3 scripts/pack_maxima.py hex 2 8 60
python3 scripts/pack_maxima.py cube 2 8 40
python3 scripts/pack_maxima.py ball 2 8 40
python3 scripts/pack_maxima.py tor3 2 5 40
