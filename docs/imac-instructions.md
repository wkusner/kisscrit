# Running the completeness computation on the iMac

How to run the N = 12 completeness computation on the late-2015 iMac (4 cores, 3.2 GHz Intel Core i5), from Terminal only. This replaces GitHub Desktop.

2026-10-11.

## What to download

Download both files from the chat into **~/Downloads**:

- `mac_run.sh`
- `kisscrit-update.bundle`

## The one line

Open Terminal (Applications → Utilities → Terminal) and run:

```
bash ~/Downloads/mac_run.sh
```

You will be asked for input only once. When the GitHub login starts:

1. Terminal shows an 8-character code.
2. Press Return, and your browser opens.
3. Paste the code there and approve.

After that, leave it running. The script keeps the Mac awake while it works.

## What the script does

Every step is skipped if it is already done.

1. **Installs a private Python (Miniforge) in `~/miniforge3`.** No admin password, no Homebrew, and nothing changes outside that folder.
2. **Creates a Python environment called `kisscrit`.** It contains numpy, scipy, networkx, mpmath, numba, git and the GitHub command-line tool `gh`. This takes about 1.5 GB of disk.
3. **Logs in to GitHub** in the browser, once.
4. **Gets the code.**
   - If your account has no repository called `kisscrit`, it creates one (**private**) from the bundle and pushes it.
   - Otherwise it clones it.
   - To use a different repository, put `REPO=critconfig` in front of the line (or whatever its name is). The code may also sit in a subfolder of that repository.
5. **Runs the self-test**, which takes 1–3 minutes.
6. **Phase A: generates candidate graphs.** It uses the four icosahedral-type triangulations, in parallel, then commits and pushes the result.
7. **Phase B: certifies the candidates in 40 shards.** It claims one shard at a time, uses all 4 cores, and commits and pushes after each shard.
8. **Prints the summary** and pushes it.

## Stopping, resuming, and checking progress

- **To stop:** press Control-C, or close the window.
- **To resume:** run the same line again. Finished work is skipped, and a shard that was interrupted is finished before a new one is claimed.
- **To check progress:** logs are in `~/kisscrit-run/logs/`. Every finished shard appears on GitHub as a commit.

## Expected time

These are rough estimates; phase A gives the real count.

- **Phase A:** about 2–3 hours. Its longest part is the icosahedron itself, which has about 1.2 million tilings.
- **Phase B:** probably between half a day and a few days on the iMac. It depends on how many distinct candidates phase A finds.
  - In tests, most candidates take well under a second.
  - About 1–2% hit the budget, at roughly 3–10 minutes each.

## Adding another machine later

Run the same line on the other Mac. It claims different shards, and git keeps the results apart. Claude cloud sessions can also run the manual commands in `docs/completeness-runbook.md`.

## What the run proves

A clean finish proves only part of completeness. Completeness also needs:

- **The hull lemma.** The run covers only the four icosahedral-type hulls, so we also have to prove that no other hull triangulation occurs above r = 1.
- **Candidates that hit the budget.** These need the stronger certificate.
- **The range 1 < r < 1.0001.** The run doesn't cover it.

The run tells us exactly how big the second item is.

## If something goes wrong

- **"push failed":** results stay committed on the iMac and are pushed on the next run.
- **"Could not claim a shard":** usually the network. Rerun.
- **To start from scratch:** delete `~/kisscrit-run`. If you also want a fresh Python, delete `~/miniforge3`.

The script has been tested end to end on Linux, with a local git remote standing in for GitHub. The Mac-specific parts have not been run on a Mac yet: the Miniforge install, `caffeinate`, and the browser login.
