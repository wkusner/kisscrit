# Setup guide: the kisscrit `full` run (accounts, iMac, progress, Zenodo)

Purpose: everything Wöden does by hand for the restarted completeness run, in the order he does it, with the exact command, what success looks like, and what to do if it fails.

Date: 2026-10-16. Version 1.1 (repo copy of the Project draft v1 of 2026-10-15). This copy wins when the two differ. The code changes it in the same commit as any code that changes a manual step.
- Written before the new run script exists. Steps marked **(confirm)** depend on a website's current layout, which was not checked.

## The plan in one paragraph

A Claude Code session rebuilds the pipeline: a fast certified checker in C on FLINT/Arb, faster candidate generation, and certificates with witnesses. It then hands you one script. That script runs the `full` route on a machine you choose (the iMac, or a Linux server if the iMac is too slow). While it runs:
- every finished piece is committed to GitHub;
- a progress file there is updated at least hourly;
- the certificate files are uploaded to a private Zenodo **draft**.

At the end you publish the Zenodo record (this gives a DOI, and cannot be undone), make the repository public, and cite both in the paper.

| When | Part |
|---|---|
| Now (about 30 minutes, no hurry) | A: accounts and tokens, iMac preparation |
| When the code session says the run script is ready | B: stop the old run, install, test, start |
| During the run | C: watching it, and fixing common problems |
| At the end | D: publish |

---

## Part A: now

### A1. Zenodo account

1. Go to https://zenodo.org and sign up, or log in with GitHub or ORCID. Using GitHub makes the code-archiving step in D3 easier.
2. Link your ORCID in your Zenodo profile if you have one. It goes into the record's creator field.

**Success:** you can see "Upload" or "New upload" in the top bar.

### A2. Zenodo sandbox account (for tests)

The sandbox is a separate test copy of Zenodo. Its DOIs are fake (prefix 10.5072) and it can be wiped at any time. Every test upload goes there first.

1. Go to https://sandbox.zenodo.org and register **separately**. Your zenodo.org login does not carry over.

### A3. Two access tokens

Make one token on each site. **(confirm the menu names)**

1. Click your account menu (top right), then **Applications**, then **Personal access tokens**, then **New token**.
2. Name: `kisscrit-run` (sandbox: `kisscrit-run-sandbox`).
3. Scopes: tick **deposit:write** and **deposit:actions**.
   - deposit:write lets the script create the draft and upload files.
   - deposit:actions lets it edit or discard drafts. The script will never publish; you do that by hand in D2.
4. Copy the token. It is shown only once.

**Don't paste tokens into a chat, an email or a file in the repo.** They go only into the iMac's Keychain (A4).

### A4. Store the tokens in the iMac's Keychain

In Terminal on the iMac, run one line per token. Each line asks for a "password": paste the token and press Return.

```
security add-generic-password -a "$USER" -s kisscrit-zenodo -w
security add-generic-password -a "$USER" -s kisscrit-zenodo-sandbox -w
```

To check them without printing them (each line should print a number around 60):

```
security find-generic-password -s kisscrit-zenodo -w | wc -c
security find-generic-password -s kisscrit-zenodo-sandbox -w | wc -c
```

If a token is wrong or expired, replace it by adding `-U` before `-w`.

The first time the run script reads a token, macOS may ask whether "security" may use the Keychain item. Choose **Always Allow**.

*If you run on Linux instead:* put each token in `~/.config/kisscrit/zenodo.token` (or `zenodo-sandbox.token`), then run `chmod 600` on that file.

### A5. Prepare the iMac

1. **Apple's Command Line Tools** (the C compiler). They are needed for the new checker.

   ```
   xcode-select --install
   ```

   A window opens; click Install. If it says they are already installed, that's fine.

   **Success:** `clang --version` prints a version.

2. **Disk space.**

   ```
   df -h ~
   ```

   Keep at least 30 GB free. The certificate files are written locally before they are uploaded. The code session will tell you the measured size per candidate after a pilot (milestone F5); revise this number then.

3. **Stay logged in.** The run is a LaunchAgent (as now). It runs only while your user is logged in, and it survives Terminal closing and reboots once you log in again.
   - Either turn on automatic login: System Preferences → Users & Groups → Login Options. **(confirm on macOS 12)**
   - Or log in after any reboot.

4. **Leave the current run alone for now.** It is checking the old `cond` certificates with the old Python checker. Its results are informational only. Part B stops it.

### A6. Which machine runs `full`: decide after the pilot

The iMac has 4 cores. How many candidates the `full` route has is not known until the new generator runs (milestone F6). The code session will report that count with measured rates, and then you choose:

- **The iMac.** Free, and set up already. Possibly weeks of wall-clock time.
- **A departmental or college Linux server.** Ask your IT or computing contact whether there is a machine with many cores you can use for a few days. The script will run on Linux.
- **A rented cloud machine**, paid by the hour.

The script is the same in each case, and the GitHub progress reporting and Zenodo uploads work identically.

---

## Part B: when the code session says the run script is ready

The code session will give you the script's exact name and download location. Below it is called `run.sh`.

### B1. Look at the old run once

```
tail -n 30 ~/kisscrit-run/logs/run.log
tail -n 30 ~/kisscrit-run/logs/check.log
```

Note where it is. If the summary in `check.log` says it finished, the old checker's results have been pushed.

### B2. Stop the old run

```
launchctl bootout gui/$(id -u)/org.kisscrit.run
```

**Success:** `launchctl list | grep kisscrit` prints nothing.

Do **not** delete `~/kisscrit-run`. It holds the old run's logs and its clone of the repository.

### B3. Install and self-test

```
bash ~/Downloads/run.sh setup
```

This does the following:
- installs or updates the private Python environment, now including FLINT and the compiler tools;
- checks the GitHub login, which may open your browser once, as before;
- compiles the checker and runs the self-tests;
- does a **sandbox Zenodo round trip**: it creates a draft on sandbox.zenodo.org, uploads a small file, checks it, and deletes it.

**Success:** the last line says the self-test passed, including Zenodo sandbox.

**If it fails:**
- *"clang not found"* or a compiler error: redo A5.1.
- *"Zenodo 401"* or *"403"*: the sandbox token is wrong or lacks a scope. Redo A3 and A4 for the sandbox.
- Any other error: copy the last 30 lines of the output into the code chat.

### B4. Start the run

```
bash ~/Downloads/run.sh start --route full
```

After a short check this hands the run to launchd, as the old script did. You can close Terminal.

---

## Part C: during the run

### C1. Watch it on GitHub (from anywhere, including your phone)

- **`results/full/progress.json`**, updated at least hourly. It shows the phase, shards done out of the total, candidates by status, boxes, CPU time and the measured rate. It has no time estimate on purpose, since earlier estimates were wrong.
- **Commits.** One per finished shard, from the iMac.
- **`results/full/zenodo-manifest.json`**. One entry per certificate file uploaded to Zenodo, with its checksum.

### C2. Watch it on the iMac

```
bash ~/Downloads/run.sh status
tail -f ~/kisscrit-run/logs/run.log
```

Press Control-C to leave `tail`; that doesn't stop the run.

### C3. Watch the Zenodo draft

On zenodo.org, go to your uploads. **(confirm the page name)** The run's draft is listed there, and its files grow as shards finish.

**Don't press Publish.**

### C4. Common problems

| Symptom | Likely cause | What to do |
|---|---|---|
| No new commits for hours, and `status` shows the same shard | a long shard, or a hang | Compare the log's last lines over an hour. If nothing changes, run `run.sh stop`, then `run.sh start --route full`. The run resumes. |
| Log shows "push failed" | expired GitHub login, or no network | `gh auth status`, then `gh auth login` if needed. The script retries on the next shard. |
| Log shows "Zenodo 401" | token expired or revoked | Make a new token (A3), then store it with `-U` (A4). The script retries. |
| "No space left on device" | disk full | Free space. Files already uploaded to Zenodo can be removed locally only if this guide says that is safe. |
| The Mac rebooted | power, or updates | Log in. The run restarts by itself. |

To stop for good: `bash ~/Downloads/run.sh stop`.

---

## Part D: at the end

### D1. Check that the run is complete

```
bash ~/Downloads/run.sh status --final
```

**Success:** every candidate is verified by the new checker, or is terminal and identified. The gen, certificate and check sets agree. Every certificate file is in the Zenodo manifest.

Anything else means the proof is not finished: bring the output to the code chat.

### D2. Publish the Zenodo record (irreversible)

1. Open the draft on zenodo.org.
2. Fill in the metadata. The code session will give you the text.
   - Title.
   - Creators, with your ORCID.
   - Description.
   - Upload type: Dataset.
   - License: CC BY 4.0 is the usual choice for data. Your call.
   - Related identifiers: the GitHub repository and, later, the paper.
3. Press **Publish**. You get a DOI.
   - You can correct files only within 30 days of publishing.
   - After that, changes need a new version, which gets its own DOI under the same concept DOI.

### D3. Archive the code

This is checklist item B4. The repository becomes public at this point.
1. Add `CITATION.cff` and `.zenodo.json`. The code session prepares both.
2. On Zenodo, go to account menu → GitHub, and switch the repository on. **(confirm the menu location)** Whether Zenodo can archive a private repository was not confirmed; plan on making it public first.
3. Make a GitHub release, for example `v1.0`. Zenodo archives that release and gives the code its own DOI.

### D4. Cite both DOIs in the paper

The paper cites the data DOI (D2) and the code DOI (D3). Both go into the completeness section and the certification appendix.

---

## Reference

- **Zenodo limits** (help.zenodo.org and developers.zenodo.org, checked 2026-10-15):
  - 100 files and 50 GB per record, raisable to 200 GB through a storage-quota request;
  - the API needs the scopes deposit:write (create and upload) and deposit:actions (publish, edit, discard);
  - publishing cannot be undone.
- **GitHub release assets:** each file under 2 GiB, up to 1,000 per release, with no total limit (docs.github.com, checked 2026-10-15). This is a fallback store if Zenodo's per-record size is not enough.
- **Where things live:**
  - code and progress: GitHub `wkusner/kisscrit`;
  - certificates: the Zenodo draft, then the published record;
  - tokens: the iMac Keychain;
  - logs: `~/kisscrit-run/logs/`.

## Maintenance

Version 1.1, 2026-10-16. Revise:
- when the run script exists (exact names and flags);
- after the F5 pilot (disk-space figure);
- after the machine choice (A6).
