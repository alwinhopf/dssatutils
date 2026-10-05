#!/usr/bin/env bash
# Exact offline Python gate; native R is a separate CI parity lane.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="$(command -v python || command -v python3)"
"$PYTHON_BIN" -m pytest tests -v -m "not live" --junit-xml=test-results.xml 2>&1 | tee fast-tests.log
