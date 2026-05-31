#!/usr/bin/env bash
# Create the PRIVATE GitHub repo and push dssatutils + the v0.1.0 tag.
# Run this ONCE, from the dssatutils repo root, after you have a way to auth.
#
# Two options — pick ONE.
set -euo pipefail
cd "$(dirname "$0")/.."   # -> dssatutils repo root

OWNER=alwinhopf
REPO=dssatutils

# ---- OPTION A: GitHub CLI (recommended) --------------------------------------
# Requires: brew install gh && gh auth login
if command -v gh >/dev/null 2>&1; then
  echo "Using gh to create private repo $OWNER/$REPO ..."
  gh repo create "$OWNER/$REPO" --private --source=. --remote=origin --push
  git push origin v0.1.0
  echo "Done (gh). Repo is private; v0.1.0 pushed."
  exit 0
fi

# ---- OPTION B: manual (no gh) ------------------------------------------------
# 1. In the GitHub UI create an EMPTY private repo: github.com/$OWNER/$REPO
#    (no README/license/gitignore — this repo already has them).
# 2. Then run:
echo "gh not found — using manual remote. Make sure the empty private repo exists."
git remote remove origin 2>/dev/null || true
git remote add origin "https://github.com/$OWNER/$REPO.git"
git branch -M main
git push -u origin main
git push origin v0.1.0
echo "Done (manual). If push prompts for credentials, use a PAT with 'repo' scope."
