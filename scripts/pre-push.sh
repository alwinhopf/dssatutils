#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# scripts/pre-push.sh — Fast local pre-push test gate for dssatutils
# Replicates the fast-unit PR CI lane locally.
# -----------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=== [dssatutils pre-push] Running offline Python test suite ==="
if command -v pytest >/dev/null 2>&1; then
  pytest tests -v -m "not live" --junit-xml=test-results.xml
elif command -v python3 >/dev/null 2>&1; then
  python3 -m pytest tests -v -m "not live" --junit-xml=test-results.xml
else
  echo "::error:: Neither pytest nor python3 found in PATH" >&2
  exit 1
fi

echo ""
echo "=== [dssatutils pre-push] Checking R offline test suite ==="
RSCRIPT_BIN=""
if [ -n "${RSCRIPT:-}" ] && [ -x "$RSCRIPT" ]; then
  RSCRIPT_BIN="$RSCRIPT"
elif command -v Rscript >/dev/null 2>&1; then
  RSCRIPT_BIN="$(command -v Rscript)"
elif [ -x "/usr/local/bin/Rscript" ]; then
  RSCRIPT_BIN="/usr/local/bin/Rscript"
elif [ -x "/opt/homebrew/bin/Rscript" ]; then
  RSCRIPT_BIN="/opt/homebrew/bin/Rscript"
fi

if [ -n "$RSCRIPT_BIN" ]; then
  echo "Using Rscript: $RSCRIPT_BIN"
  "$RSCRIPT_BIN" --vanilla -e "
    if (requireNamespace('pkgload', quietly = TRUE) && requireNamespace('testthat', quietly = TRUE)) {
      pkgload::load_all('.', quiet = TRUE)
      testthat::test_dir('tests/testthat', reporter = 'summary', stop_on_failure = TRUE)
    } else {
      cat('Skipping R tests: pkgload or testthat not installed in R environment.\n')
    }
  "
else
  echo "Rscript not found; skipping R test suite."
fi

echo ""
echo "=== [dssatutils pre-push] All fast unit tests passed! ==="
