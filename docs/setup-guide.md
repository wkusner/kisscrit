# Setup guide: the kisscrit `full` run (Linux server over SSH, GitHub progress, Zenodo)

Purpose: everything Wöden does by hand for the restarted completeness run, in order: the exact command, what success looks like, and what to do if it fails.

This is the maintainer's runbook, written for one person and one server; it is kept in the repository as the record of the manual steps. "The Project guide" below is an internal document, available on request.

Date: 2026-10-16. Version 2.3 (repo copy of the Project guide v2.1; this copy wins when the two differ, and the code changes it in the same commit as any code that changes a manual step).
- v2.3 (2026-10-22, F8 M1): Part B written for `scripts/full/run.sh` version 1 (setup, selftest, start, status, stop; no pushes or Zenodo yet); C2, C4 and D1 use its commands.
- v2.1: the server runs Ubuntu 26.04 and Wöden has no admin rights. Step A4 uses `scripts/full/kisscrit-preflight.sh`, which tests the server and drafts the request to the sysadmin.
- v2: the run machine is a dedicated Linux server reached by SSH. It has no scheduler, Wöden is probably its only user, and it has 32 or 64 cores. The iMac becomes the fallback. v1 (2026-10-15) assumed the iMac.
- Written before the run script exists. Steps marked **(confirm)** depend on a website's current layout, which was not checked.

## The plan in one paragraph

A Claude Code session rebuilds the pipeline: a fast checker in C on FLINT/Arb (every error bound from Arb), faster candidate generation, and certificates with witnesses. It then hands you one script, `run.sh`. You copy that script to the server, run `setup` once, then `start`. The run survives logging out and reboots. While it runs:
- every finished piece is committed to GitHub;
- a progress file there is updated at least hourly;
- the certificate files are uploaded to a private Zenodo **draft**.

At the end you publish the Zenodo record (this gives a DOI and cannot be undone), make the repository public, and cite both in the paper.

| When | Part |
|---|---|
| Now, while you set up the server (no hurry) | A: accounts, tokens, server checks |
| When the code session says `run.sh` is ready | B: install, test, start |
| During the run | C: watching it, and fixing common problems |
| At the end | D: publish |
| Any time (optional) | E: the iMac's old run |

---

## Part A: now

### A1. Zenodo account

1. Go to https://zenodo.org. Sign up, or log in with GitHub or ORCID. GitHub makes step D3 easier.
2. Link your ORCID in your Zenodo profile if you have one.

**Success:** "Upload" or "New upload" appears in the top bar.

### A2. Zenodo sandbox account (for tests)

The sandbox is a separate test copy of Zenodo. Its DOIs are fake (prefix 10.5072), and it can be wiped at any time. Every test upload goes there first.

1. Register **separately** at https://sandbox.zenodo.org. Your zenodo.org login does not carry over.

### A3. Two access tokens (in your laptop's browser)

**(confirm the menu names)**

1. On each site, open your account menu (top right), then **Applications**, then **Personal access tokens**, then **New token**.
2. Name it `kisscrit-run`, or `kisscrit-run-sandbox` on the sandbox.
3. Scopes: tick **deposit:write** and **deposit:actions**.
   - The script uses deposit:write to create the draft and upload.
   - It uses deposit:actions only to edit or discard drafts.
   - It never publishes. You do that by hand in D2.
4. Copy the token. It is shown only once. Keep it in your clipboard for A5.

**Don't paste tokens into a chat, an email, a file in the repo, or a shared folder.**

### A4. Server checks (over SSH): run the preflight script

The server runs Ubuntu 26.04 and you have no admin rights. So the first step is to find out what the sysadmin must do. `scripts/full/kisscrit-preflight.sh` (also in the Meta inbox and the Project) does that:
- it tests everything the run needs, without sudo;
- it saves a report;
- it ends with a ready-to-send email to the sysadmin that lists only what is missing;
- with `--try-install`, it also builds the repo's C kernel with the server's compiler and runs its tests (below).

This version was tested only in a cloud container (Ubuntu 24.04), not on the server: the first run on the server is yours.

From your Mac:

```
scp ~/Meta_Inbox/kisscrit-preflight.sh you@server:~/
# for --try-install only (needed while the repository is private; once it is public a plain `git clone` works): from your laptop, inside your clone of kisscrit; about 360 MB, mostly the certificates
git archive HEAD csrc kisscrit tests scripts/completeness pyproject.toml results/completeness/certificates | ssh you@server 'mkdir -p ~/kisscrit-src && tar -x -C ~/kisscrit-src'
ssh you@server
bash kisscrit-preflight.sh                  # quick checks, about a minute
bash kisscrit-preflight.sh --try-install --repo ~/kisscrit-src    # optional: also installs FLINT and a C compiler into ~/kisscrit-preflight, builds and tests the C kernel
```

Options (all optional):
- `--need-disk-gb N`: the free space the run needs (default 50; an estimate, to be revised after milestone F5/F6). The report lists free space for your home, `/tmp` and the filesystem of the results directory, and the sysadmin request uses N.
- `--results-dir DIR`: where the run will keep its files (default `~/kisscrit-run`). The disk requirement is checked on its filesystem; nothing is created there.
- `--repo DIR`: the copy of the repository for `--try-install` (above). Without it, the script uses the repository it sits in, if it is run from one; otherwise the build and tests are reported as not run.

`--try-install` needs GitHub and conda-forge to be reachable. It answers the open questions: whether the toolchain installs without root, whether conda-forge's `libflint` ships its headers, and whether the C kernel builds and passes its tests with the server's compiler. After the environment is created it:
- builds `csrc/` with the environment's compiler and the repo's flags (`-O2 -fno-fast-math -ffp-contract=off`) and runs `make -C csrc test`;
- scans the built library for fused multiply-add instructions (there must be none; on a machine other than x86-64 or aarch64, or without `objdump`, it says "not checked");
- runs `tests/test_carith.py`, `test_cmodel.py`, `test_creplay.py` and `test_ckernel.py`. The result is "ok" only if every test passes; skipped tests, or missing certificates, give a warning, never an "ok";
- prints a build record for the report: machine, OS and glibc, compiler, FLINT, Python, numpy, numba and python-flint versions, cores and memory.

In the cloud container it took 271 s (the four pytest files 163 s, 81 tests) and left 3.1 GB in `~/kisscrit-preflight`; the server may differ.

Delete its test install afterwards with `rm -rf ~/kisscrit-preflight`.

The report is saved to `~/kisscrit-preflight-report.txt`.
- Send the "Draft request to the sysadmin" section to your sysadmin.
- Paste the whole report into the code chat. It contains your user name and the host name, nothing secret.

The rest of this section explains what the script checks, in case you want to do it by hand. Log in (`ssh you@server`) and run:

```
nproc                 # cores: 32 or 64 expected
free -g               # memory in GB
df -h ~ /tmp          # free disk in your home directory and in /tmp
which git tmux crontab curl
curl -s -o /dev/null -w '%{http_code}\n' "https://github.com/conda-forge/miniforge/info/refs?service=git-upload-pack"
curl -s -o /dev/null -w '%{http_code}\n' "https://zenodo.org/api/records?size=1"
curl -s -o /dev/null -w '%{http_code}\n' "https://sandbox.zenodo.org/api/records?size=1"
loginctl show-user "$USER" -p Linger
```

What you want to see:
- **Disk:** at least 50 GB free in `~` (or where `--results-dir` points). Certificates are written locally before they upload. This figure is an estimate; the code session revises it after milestones F5/F6.
- **`which`:** prints a path for each of git, tmux, crontab and curl.
  - Missing ones go in the sysadmin request (`sudo apt install tmux cron curl`).
  - git and gh also come with the run environment, so tmux and cron are the ones that matter.
- **The three `curl` lines:** each prints an HTTP status code. The script classifies them the same way: 2xx and 3xx are fine; 4xx and 5xx mean the host answered but not with success (a proxy or firewall often answers 403 or 407 for a blocked host), and the script warns with the code; `000`, an error or a hang means the server cannot reach that site, and the script reports a failure.
  - The first line goes to github.com through the git protocol's own address, because the bare `https://github.com` returned 400 in the cloud container, where it was wrongly counted as reachable.
  - If one hangs or errors, ask IT to allow outbound HTTPS to github.com, zenodo.org and sandbox.zenodo.org.
  - This is the one thing that would block the run.
- **`Linger=yes`:** lets the run continue as a background service after you log out.
  - If it says `no`: you have no root, so this goes in the sysadmin request (`sudo loginctl enable-linger <you>`).
  - The alternative is a `tmux` session plus a restart-on-reboot cron line, which needs tmux and cron installed and cron allowed for your account. That works just as well; the sysadmin can choose.

Send the code chat the output of these commands, with nothing secret in it.

### A5. Store the tokens on the server

```
mkdir -p ~/.config/kisscrit && chmod 700 ~/.config/kisscrit
cat > ~/.config/kisscrit/zenodo.token
```

Paste the **real** token, press Return, then Control-D. Then the sandbox token:

```
cat > ~/.config/kisscrit/zenodo-sandbox.token
```

Paste it, Return, Control-D. Then lock down both files:

```
chmod 600 ~/.config/kisscrit/*.token
wc -c ~/.config/kisscrit/*.token      # each about 60; prints sizes, not the tokens
```

The script refuses to run if these files are readable by anyone but you.

As the only user, you're the only one who can read them. Root on the machine always can, which is acceptable for revocable tokens with only these scopes.

To replace a token, repeat the `cat >` line. To revoke one, use the same Applications page on Zenodo.

### A6. Optional: an SSH key from your laptop

If you don't already log in with a key:

```
ssh-copy-id you@server
```

Run this on your laptop, once. After that, `ssh you@server` doesn't ask for a password.

---

## Part B: when the code session says `run.sh` is ready

`scripts/full/run.sh` version 1 (F8 milestone M1, 2026-10-22). What it does now: setup, preflight, self-test, start, status and stop on routes `test` and `cond`, with resume, one run at a time, and survival by systemd or tmux. What it does not do yet: push to GitHub (milestone M2) and Zenodo (milestone M3). Until then, results are committed only locally, on the run branch `run/<route>-<run-id>` in `~/kisscrit-run/runs/<run-id>/tree`. Route `full` is refused until `results/full/F6-READY` exists (checklist F6). macOS is milestone M4. Tested in a cloud container only (Ubuntu 24.04, systemd and cron faked); the first server run is yours.

### B1. Copy the scripts to the server

From your laptop, inside your clone of kisscrit:

```
scp scripts/full/run.sh scripts/full/kisscrit-preflight.sh you@server:~/
```

`run.sh` needs nothing else for `setup`. Every later step runs the copy of `run.sh` inside the pinned tree.

The pinned commit must be pushed to GitHub and must contain `scripts/full/run.sh`. Pin it as a tag (recommended, e.g. `git tag run-cond-1 <hash> && git push origin run-cond-1`) or as the full 40-character hash. Branch names and short hashes are refused.

**Without GitHub on the server** (as for the pilot): send the tree with a COMMIT file, then use `--archive ~/kisscrit-src` in B2. No pushes are possible in this mode.

```
git archive HEAD | ssh you@server 'mkdir -p ~/kisscrit-src && tar -x -C ~/kisscrit-src' && git rev-parse HEAD | ssh you@server 'cat > ~/kisscrit-src/COMMIT'
```

### B2. Setup (once)

```
ssh you@server
bash ~/run.sh setup
```

It does these steps, all under your home directory, without root:
- runs the preflight (`kisscrit-preflight.sh`) and stops if anything fails: an unreachable host, or too little disk;
- installs Miniforge into `~/miniforge3`, downloading with curl or, if curl is missing, wget;
- creates the conda environment `kisscrit-run` (python 3.12, numpy, scipy, numba, python-flint, libflint, c-compiler, make, networkx, mpmath, sympy, pytest, git, gh), logged to `~/kisscrit-run/logs/conda.log`;
- logs in to GitHub. `gh` prints a one-time code. On your laptop, open https://github.com/login/device, type the code and approve. Nothing has to open a browser on the server;
- clones `wkusner/kisscrit` into `~/kisscrit-run/repo`.

With `--archive ~/kisscrit-src` it uses the archive tree instead and skips the GitHub login.

**Success:** the last lines give the cores, the memory, the survival method (`systemd` when lingering is on, otherwise `tmux`) and `next: bash run.sh selftest --commit ...`.

### B3. Self-test the pinned commit

```
bash ~/run.sh selftest --commit run-cond-1
```

In a scratch copy of the commit under `~/kisscrit-run/selftest/<hash>/`, it does these steps:
- builds `csrc/` with the environment's compiler and runs `make -C csrc test`;
- runs the whole pytest suite. `tests/test_ball_kill.py` runs in its own process, because it changes a global Arb precision that other tests depend on;
- runs `scripts/completeness/selftest.py`;
- measures the peak memory of one certify worker on route test. `start` uses this to cap the number of workers.

The Zenodo sandbox round trip comes with milestone M3. In the cloud container the self-test took about 2 minutes with a pytest subset; with the whole suite allow 10 minutes.

**Success:** `selftest passed for <hash>`. The record is `~/kisscrit-run/selftest/<hash>.json`.

For a manual `pytest -q` on the server, set `KISSCRIT_REQUIRE_ALL=1`: `tests/conftest.py` then turns every skipped test (C kernel not built, libflint, certificates or gen files absent) into a failure, except the two allowed skips in its `ALLOWED_SKIP_MARKERS`: the opt-in `DOM1_REGEN` test and the shellcheck lint of `run.sh` (reason "shellcheck not installed"; the server's conda environment has no shellcheck, and `bash -n` still runs). `run.sh selftest` sets it for its pytest step, so on the server a missing tool (tmux, flock, libflint) makes the self-test fail instead of skipping; a missing shellcheck only skips the lint (CLEAN-B, CLEAN-B2).

### B4. Start

The pilot-equivalent run on `cond`:

```
bash ~/run.sh start --route cond --commit run-cond-1 --fresh
```

- `--fresh` is needed because the repository already holds `cond` results; without them `run_shard.py` would skip every candidate. It removes them on the run branch only, in one commit, and the removed file names go to `~/kisscrit-run/runs/<run-id>/logs/fresh-removed.txt`. `main` is never touched.
- Defaults: 40 shards, budget 300,000, workers `nproc` capped by memory, 2 workers per shard. Override with `--shards`, `--budget`, `--workers`, `--workers-per-shard`.
- The reference checker `scripts/completeness/check_certificates.py` decides. `--checker SCRIPT` names another tracked script of the pinned tree, and the run record names whichever was used. The C kernel front end is not used (audit finding B-8).
- Every result is labelled "checker verdict at commit X; not certified".

Then log out. The run continues.
- **Resume:** after a stop, a kill or a crash, `bash ~/run.sh start` (no options) continues. Finished shards are skipped.
- **Damaged lines:** before each resume, `run_repair.py` moves any damaged or unpaired result lines to `~/kisscrit-run/runs/<run-id>/quarantine/`, logged in `logs/repair.jsonl`.
- **A new run:** while a run is unfinished, a new one needs `--new`.

---

## Part C: during the run

### C1. On GitHub (from anywhere, including your phone)

**Not yet:** `run.sh` version 1 commits only locally on the run branch, so nothing below appears on GitHub until milestone M2. Until then, use C2.

- **`results/full/progress.json`**, updated at least hourly. It shows the phase, shards done out of total, candidates by status, boxes, CPU time and the measured rate. There is no time estimate on purpose.
- **Commits.** One per finished shard.
- **`results/full/zenodo-manifest.json`**. One entry per certificate file uploaded to Zenodo, with its checksum.

### C2. On the server

```
bash ~/run.sh status                     # run, phase, shards, survival (systemd or tmux), last log lines
tail -f ~/kisscrit-run/logs/run.log      # Control-C leaves tail, not the run
tmux -L kisscrit attach -t kisscrit-run  # only in tmux mode; detach with Control-B then D
bash ~/run.sh stop --reason "why"        # stops cleanly; the reason goes into the run record
```

### C3. The Zenodo draft

Look in your uploads on zenodo.org **(confirm the page name)**. The draft's files grow as shards finish.

**Don't press Publish.**

### C4. Common problems

| Symptom | Likely cause | What to do |
|---|---|---|
| No new commits for hours, and `status` is stuck on one shard | a long shard, or a hang | Compare the log's last lines an hour apart. If nothing changes: `run.sh stop`, then `run.sh start` (no options). The run resumes. |
| "push failed" in the log | expired GitHub login, or the network | `gh auth status`; if needed, `gh auth login`. The script retries. |
| "Zenodo 401" | token expired or revoked | New token (A3), then replace the file (A5). The script retries. |
| "No space left on device" | disk full | Free space. Delete uploaded certificate files only if `run.sh status` lists them as uploaded and verified. |
| The server rebooted | maintenance | Nothing to do: systemd or the @reboot cron line restarts the run (not after `run.sh stop`). Check `run.sh status`. |

To stop for good: `bash ~/run.sh stop`.

---

## Part D: at the end

### D1. Check that the run is complete

```
bash ~/run.sh status --final
```

**Success means all of these hold:**
- every candidate is verified by the new checker, or is terminal and identified;
- the gen, certificate and check sets agree;
- every certificate file is in the Zenodo manifest.

Anything else means the proof is not finished. Take the output to the code chat.

### D2. Publish the Zenodo record (irreversible)

1. Open the draft on zenodo.org.
2. Fill in the metadata. The code session gives you the text.
   - Title.
   - Creators, with your ORCID.
   - Description.
   - Upload type: Dataset.
   - License: CC BY 4.0 is usual for data; your call.
   - Related identifiers: the GitHub repository, and later the paper.
3. Press **Publish**. This gives the DOI.
   - Files can be corrected only within 30 days.
   - After that, changes need a new version, with its own DOI under the same concept DOI.

### D3. Archive the code

The repository becomes public here (checklist item B4).

1. Add `CITATION.cff` and `.zenodo.json`. The code session prepares both.
2. On Zenodo, open account menu → GitHub and switch the repository on **(confirm the menu location)**. Whether Zenodo can archive a private repository was not confirmed; make it public first.
3. Make a GitHub release, for example `v1.0`. Zenodo archives it and gives the code its own DOI.

### D4. Cite both DOIs in the paper

The data DOI (D2) and the code DOI (D3) go into the completeness section and the certification appendix.

---

## Part E: the iMac's old run (optional, any time)

The iMac may still be checking the old `cond` certificates with the old Python checker. Those results are informational and play no part in the proof. You can let it finish, or stop it.

```
tail -n 30 ~/kisscrit-run/logs/check.log       # where it is
launchctl bootout gui/$(id -u)/org.kisscrit.run  # stop it
launchctl list | grep kisscrit                   # prints nothing once stopped
```

Keep the `~/kisscrit-run` folder on the iMac. If the server falls through, the iMac can run `run.sh` too. For that it needs Apple's Command Line Tools (`xcode-select --install`), and the tokens go in its Keychain:

```
security add-generic-password -a "$USER" -s kisscrit-zenodo -w
security add-generic-password -a "$USER" -s kisscrit-zenodo-sandbox -w
```

Each line prompts you to paste a token.

---

## Reference

- **Zenodo** (help.zenodo.org and developers.zenodo.org, checked 2026-10-15):
  - limits: 100 files and 50 GB per record, raisable to 200 GB through a storage-quota request;
  - API scopes: deposit:write and deposit:actions;
  - publishing cannot be undone.
- **GitHub release assets** (docs.github.com, checked 2026-10-15): each file under 2 GiB, up to 1,000 per release, no total limit. This is the fallback if one Zenodo record is not enough.
- **Ubuntu 26.04 notes** (from release coverage, not yet tested on the server): `/tmp` is a tmpfs capped at half of RAM, so the run keeps scratch files under `~/kisscrit-run`; the default `sha256sum` is a slower rewrite, so the scripts hash with Python; the system Python is not used.
- **Where things live:**
  - code and progress: GitHub `wkusner/kisscrit`;
  - certificates: the Zenodo draft, then the published record;
  - tokens: `~/.config/kisscrit/` on the server, or the iMac Keychain;
  - logs: `~/kisscrit-run/logs/` on the run machine.

## Maintenance

Version 2.3, 2026-10-22. Revise:
- when `run.sh` exists (exact names and flags);
- after the F5 pilot (the disk figure);
- after the A4 server checks (anything the server lacks);
- when run.sh milestones M2 (GitHub pushes), M3 (Zenodo) and M4 (macOS) land: C1, C3, D1 and Part E.
