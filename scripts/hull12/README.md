Scripts behind docs/hull.md (N = 12). They expect tri12.json (= data/triangulations_n12.json) in the
working directory and, for flipwalk.py, seeded.json (results/) and ts_a.json from survey output.
survey.py: theta*(T) for all 7595 triangulations (resumable, k % W == w split across W workers).
flipwalk.py: from realizations with theta > 60, flip each hull edge and re-optimize seeded at the parent.
sample_hulls.py: ascend random configurations to just above 60 degrees and tally hull types.
