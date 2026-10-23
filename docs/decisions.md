# Decisions

One line each, dated. Newest last. Only decisions Wöden stated or accepted; proposals go in `docs/open-questions.md`.

Documents named below and not in this repository (the transfer prompt v1.2, the checker specification v0.2, the pipeline specification v0.3, the control board with its decision ids D-0nn) are internal working documents, available on request; they are not copied into the repository.

- 2026-10-15 — Run the `full` route (all 7,595 triangulations). Hypothesis H is removed by it; `cond` becomes a pilot and a regression set.
- 2026-10-15 — Drop the hull lemma (checklist C1, C2) from the proof path; keep it in the paper as a conjecture with the survey, marked unused.
- 2026-10-15 — Rebuild the checker as a C kernel on FLINT/Arb, with every error bound from a certified source (spec v0.3); the Python reference checker, patched with V1–V5, stays as the oracle.
- 2026-10-15 — Report progress to GitHub during the run (`results/full/progress.json`, hourly, rates only, no ETA). Archive certificates on Zenodo as a draft the script never publishes; only Wöden publishes.
- 2026-10-15 — The range 1 < r < r_min is handled by compactness plus second-order isolation at the balanced 60° configurations (paper, Appendix B), with a band run outside the certified balls.
- 2026-10-15 — Every manual step lives in `docs/setup-guide.md`, updated in the same commit as the code that changes it.
- 2026-10-16 — The `cond` candidate count is 45,260 (distinct codes); 90,836 was the number of gen lines including duplicates and is not to be used.
- 2026-10-16 — The run machine is a dedicated Linux server (Ubuntu 26.04, direct SSH, no scheduler, 32 or 64 cores, no root for Wöden; anything needing root goes to his sysadmin). The iMac (macOS 12, 4 cores) is the fallback. The run script targets Linux first: no-root install (Miniforge, `libflint`, conda-forge `c-compiler`, git, gh); Zenodo tokens from `~/.config/kisscrit/zenodo.token` and `zenodo-sandbox.token` (mode 600, looser refused); survival by `systemd --user` when lingering is on, else `tmux` plus a `@reboot` crontab line; workers default to `nproc`, capped by the memory measured in the self-test.
- 2026-10-16 — Ubuntu 26.04 specifics (from release coverage, not yet tested on the server): `/tmp` is a tmpfs capped at half of RAM, so scratch files go under `~/kisscrit-run`; the default `sha256sum` is the slower uutils rewrite, so hashes are computed with Python `hashlib`; the system Python (3.14) is not used.
- 2026-10-16 — Preflight: Wöden runs `kisscrit-preflight.sh` on the server and pastes its report; `run.sh setup` reuses its checks (the script lives at `scripts/full/kisscrit-preflight.sh`).
- 2026-10-16 — Code standards, "readable by a referee" (transfer prompt v1.2): every file carries a header (purpose, trusted or untrusted, lemmas and appendix sections implemented); every trusted function states its mathematical claim; every inequality and rounding direction cites its justification; plain code with paper-notation names; `docs/code-tour.md` maps paper statements to code; a fresh agent reviews readability at F1, F4 and F6 before the item is ticked; a test enforces the headers. Applies from F1 onward, including the reference checker.
- 2026-10-20 — F4 performance gate: accept the spec v0.3 §4.4 target (checker CPU time ≤ search CPU time on the same candidates); measured ratio 0.24 on 242 fresh `cond` certificates, about 0.83 without the two numba-compile candidates (2-core Xeon VM). F4 ticked.
- 2026-10-21 — The `cond` certificates stay in git as the pilot record. The Zenodo rule of 2026-10-15 (certificates on Zenodo, not in git) applies to the `full` run. (Wöden, 15:37 CT; board D-020.)
- 2026-10-20 — The `cond` certificates are not regenerated wholesale for the F4 differential: compare boxes and leaves on the existing certificates plus a fresh sample of at least 200 from the current tag. (Wöden, 10:57 CT; board D-016.)
- 2026-10-20 — F4 scope change: `kc_hull` in `csrc/kc_arith.c` follows Python `min`/`max` (keeps the first argument on ties and NaN); own commit, regression tests from the four failing inputs. (Wöden, 12:33 CT; board D-017.)
- 2026-10-21 — C4 is split into C4a (Krawczyk balls for regular terminal boxes) and C4b (an isolation argument for the singular terminal boxes: 458 in the `cond` run, 456 DOD and 2 M, reported by DES-F2, plus the lattice rule D5). (Wöden, 21:35 CT; board D-021.)
- 2026-10-21 — The paper source moves into the repo (`paper/`); `01_Paper` keeps archived drafts. Not started. (Wöden, 21:35 CT; board D-022.)
- 2026-10-21 — Fable 5.1 reviews the V3-FIX and BALL-KILL diffs after DES-F2 (done by AUD-2, 2026-10-23). (Wöden, 21:35 CT; board D-023.)
- 2026-10-21 — The repo bib is canonical; the other bib is merged into it and the old one goes to `_Archive`; both are listed with their differences before the merge. Not started. (Wöden, 21:35 CT; board D-024.)
- 2026-10-22 — C4 design: the per-leaf trichotomy (lemma L23) replaces the cluster ball of D1; the BALL-KILL ball data stay as data. (Wöden, "Accept all", 09:15 CT; board D-025.)
- 2026-10-22 — Exact M: store 80-digit intervals now (the residual form of L24); exact re-verification later, optional. (Board D-026.)
- 2026-10-22 — A new critical configuration whose certificate the checker has accepted is recorded and the run continues; it is flagged in `progress.json` and the manifest and never reads as "complete". (Board D-027.)
- 2026-10-22 — Accept the √N loss of the Procrustes step for now; the quantitative lemma only if the numbers need it. (Board D-028.)
- 2026-10-22 — DOM-1 before any further run (done: V1b, f3e8d0d). (Board D-029.)
- 2026-10-22 — The manifest refuses ("not complete") when two check generations for a code disagree; no newest-wins. Not implemented (DOM-3). (Board D-030.)
- 2026-10-22 — Replace checklist C3 by R1–R4; run R1 first (53 candidates at 3M/30M budget, about 38 CPU-hours), which needs a `--codes` option in `run_shard.py` (C3-1). (Board D-031.)
- 2026-10-22 — H5 (an explicit list of unresolved codes) is accepted as the fallback paper statement only. (Board D-032.)
- 2026-10-23 — All fixes proposed by the AUD-2 audit are approved (sync document items S-01 to S-20 and S8, review section 2), as docs commits and the listed data and docstring corrections. (Wöden, 01:44 CT; board D-033.)
- 2026-10-23 — Be cautious about the word "certified": it is reserved for a certificate accepted by the checker; interval verification is called "verified in interval arithmetic"; existing file names and JSON keys are not renamed; nothing new says "certified". (Wöden, 01:44 CT; board D-034.)
- 2026-10-23 — Publishing the repository is acceptable after a double check of the fixes and a last scan; the visibility switch, the Zenodo link and the release tag are Wöden's. (Wöden, 01:44 CT; board D-035.)
