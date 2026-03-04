#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
UNIT_TEST_DIR="$ROOT_DIR/src/unit_tests"

test_files=()
while IFS= read -r file; do
  test_files+=("$file")
done < <(find "$UNIT_TEST_DIR" -maxdepth 1 -type f -name 'test_*.f90' | sort)

if [ "${#test_files[@]}" -eq 0 ]; then
  echo "No unit-test files found under src/unit_tests; skipping assertion check."
  exit 0
fi

missing_assertions=()
for file in "${test_files[@]}"; do
  if ! rg -qi "call[[:space:]]+assert_[a-z0-9_]+[[:space:]]*\\(" "$file"; then
    missing_assertions+=("$file")
  fi
done

if [ "${#missing_assertions[@]}" -gt 0 ]; then
  echo "FAIL: assertion-free unit tests are not allowed."
  echo "Add at least one assertion call (e.g., call assert_true/assert_equal_*) in:"
  printf '  - %s\n' "${missing_assertions[@]}"
  exit 1
fi

echo "Assertion presence check passed for ${#test_files[@]} unit-test file(s)."
