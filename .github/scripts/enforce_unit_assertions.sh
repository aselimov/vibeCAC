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
smoke_only_logic_heavy=()
missing_negative_cases=()
for file in "${test_files[@]}"; do
  if ! rg -qi "call[[:space:]]+assert_[a-z0-9_]+[[:space:]]*\\(" "$file"; then
    missing_assertions+=("$file")
  fi

  case "$(basename "$file")" in
    test_compute.f90|test_potential.f90|test_integration.f90|test_box.f90|test_temp.f90|test_neighbors.f90|test_thermo.f90|test_berendsen.f90|test_math.f90)
      total_asserts="$(rg -io "call[[:space:]]+assert_[a-z0-9_]+[[:space:]]*\\(" "$file" | wc -l | tr -d ' ')"
      trivial_asserts="$(rg -io "call[[:space:]]+assert_true[[:space:]]*\\([[:space:]]*\\.true\\.|call[[:space:]]+assert_false[[:space:]]*\\([[:space:]]*\\.false\\." "$file" || true)"
      trivial_asserts="$(printf '%s\n' "$trivial_asserts" | wc -l | tr -d ' ')"
      if [ "$total_asserts" -gt 0 ] && [ "$total_asserts" -eq "$trivial_asserts" ]; then
        smoke_only_logic_heavy+=("$file")
      fi
      ;;
  esac

  case "$(basename "$file")" in
    test_input_validation.f90|test_time.f90|test_atom_types.f90)
      if ! rg -qi "call[[:space:]]+assert_[a-z0-9_]+[[:space:]]*\\([^\\n]*'(?:[^']*)(invalid|reject|rejected|fail|error|missing|non[- ]?numeric|zero|negative)" "$file"; then
        missing_negative_cases+=("$file")
      fi
      ;;
  esac
done

if [ "${#missing_assertions[@]}" -gt 0 ]; then
  echo "FAIL: assertion-free unit tests are not allowed."
  echo "Add at least one assertion call (e.g., call assert_true/assert_equal_*) in:"
  printf '  - %s\n' "${missing_assertions[@]}"
  exit 1
fi

if [ "${#smoke_only_logic_heavy[@]}" -gt 0 ]; then
  echo "FAIL: pure smoke tests are not allowed for logic-heavy routines."
  echo "Replace tautological assertions with behavioral checks (expected vs actual state/output) in:"
  printf '  - %s\n' "${smoke_only_logic_heavy[@]}"
  exit 1
fi

if [ "${#missing_negative_cases[@]}" -gt 0 ]; then
  echo "FAIL: validator/parser unit tests must include at least one explicit negative-case assertion."
  echo "Add a failing-input assertion (e.g., invalid/rejected/missing/error path) in:"
  printf '  - %s\n' "${missing_negative_cases[@]}"
  exit 1
fi

echo "Unit-test assertion checks passed for ${#test_files[@]} file(s)."
