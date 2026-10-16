"""Completeness of the census of critical configurations of 12 points above r = 1 (work in progress).

Status: untrusted (exploration or tooling; no part of the proof).

Pipeline (see docs/completeness-runbook.md):
  candidates.py  enumerate candidate stressed graphs (plane graphs with 3-, 4-, 5-gon faces) from triangulations
  canon.py       canonical codes of plane graphs (deduplication across triangulations and machines)
  model.py       face-angle model, reference implementation (numpy intervals)
  fast.py        compiled kernels (numba): corner enclosures, contraction, Gordan/Motzkin exclusion, B&B
  runner.py      one candidate, all pentagon charts
"""
