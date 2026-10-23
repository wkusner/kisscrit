# attic: retired prototypes (not maintained)

This folder holds early prototypes of the completeness search that the current code no longer uses. They were moved here, with their git history, in the CLEAN-A session (2026-10-22, checklist B1); nothing imports them.

| file | what it was |
|---|---|
| `kisscrit/completeness/ad.py` | forward-mode automatic differentiation over intervals or floats |
| `kisscrit/completeness/bnb.py` | an early branch and bound over the face-parameter box of a candidate stressed graph |
| `kisscrit/completeness/contract.py` | a linear-relaxation contractor for the vertex angle-sum equations |
| `scripts/completeness/lpfilter.py` | an LP relaxation test for a candidate face structure |
| `scripts/completeness/rattler_check.py` | an early rattler check |

These files are **not maintained**, are not tested, are not part of the trusted base and are not part of the code tour (`docs/code-tour.md`). Nothing in `STATUS.md` or the paper rests on them. The search in use is `kisscrit/completeness/` (`fast.py`, `runner.py`) and the checker is `kisscrit/core/`.
