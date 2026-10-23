# Pilot snapshot, route `cond`, 2026-10-22

Text files copied from the server pilot run (commit 322f8bc of the pilot script's work tree, 64 cores, started 2026-10-22 08:17 CDT, done 12:47 CDT) as published by the pilot's `publish-v1` script on the branch `pilot/cond-2026-10-22` (last snapshot commit 9268c00, which forked from f3e8d0d; the branch was removed from the public remote and is kept as a git bundle outside the repository).

**Edits for publication (2026-10-23):** in the `.txt` and `.log` files the server's user home path is replaced by `<run-home>` and the server's host name by `<server>`. Nothing else was changed; `summary-cond.json`, `checked-lines.txt` are byte-identical to the branch.

**Status of the numbers:** computed on the server; not reproduced in this repository (the pilot's certificate files are not here). STATUS.md states the counts used; the explanation of the 467/468 bookkeeping and of the three candidates that moved relative to the iMac run is AUD-2's comparison (reported, not reproduced). This is a pilot of the conditional route `cond`; it proves nothing about completeness.

**On the word "certified" in these logs:** the files are verbatim output of the pilot script. There, `certified` (as in "certified 45260" and the key `certified` of `summary-cond.json`) is the script's historical name for "records produced by the certify phase", not a certificate that the checker has accepted. The checker run is reported in `PILOT-REPORT.txt`: 45,204 of 45,260 verified, 56 not verified (structural), 0 failed boxes (computed on the server, not reproduced here; see `STATUS.md`). Nothing here establishes completeness.
