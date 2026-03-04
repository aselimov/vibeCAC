#!/usr/bin/env bash
set -euo pipefail

if [[ "${GITHUB_EVENT_NAME:-}" != "pull_request" ]]; then
  echo "coverage delta: skipping (not a pull_request event)"
  exit 0
fi

BASE_REF="${GITHUB_BASE_REF:-}"
if [[ -z "$BASE_REF" ]]; then
  echo "coverage delta: missing GITHUB_BASE_REF"
  exit 0
fi

git fetch --no-tags --depth=200 origin "$BASE_REF"
MERGE_BASE="$(git merge-base "origin/$BASE_REF" HEAD)"

map_module_to_unit_test() {
  case "$1" in
    math) echo "math" ;;
    str) echo "str" ;;
    time) echo "time" ;;
    box) echo "box" ;;
    dump_logic) echo "dump_logic" ;;
    input_validation) echo "input_validation" ;;
    atom_types) echo "atom_types" ;;
    *) echo "" ;;
  esac
}

line_coverage_for() {
  local repo_root="$1"
  local module="$2"
  local target="$3"
  local gcov_bin

  gcov_bin="$(command -v gcov-15 2>/dev/null || command -v gcov)"

  (
    cd "$repo_root/src"
    make clean >/dev/null
    make unit-test UNIT_TEST="$target" UNIT_TEST_FFLAGS="-Wall -O0 -g -fcheck=all -fbacktrace --coverage" >/dev/null
    local out pct
    out="$($gcov_bin -b -c "$module.f90" -o obj || true)"
    pct="$(printf '%s\n' "$out" | awk '/Lines executed:/{gsub("executed:","",$2); gsub("%","",$2); print $2; exit}')"
    if [[ -z "$pct" ]]; then
      echo "NA"
    else
      echo "$pct"
    fi
  )
}

CHANGED_MODULES="$(git diff --name-only "$MERGE_BASE"...HEAD -- 'src/*.f90' | sed -E 's#^src/##; s#\.f90$##' | sort -u)"
if [[ -z "$CHANGED_MODULES" ]]; then
  echo "No changed src/*.f90 modules in this PR; skipping coverage delta report."
  exit 0
fi

WORKTREE_DIR="$(mktemp -d)"
trap 'git worktree remove --force "$WORKTREE_DIR" >/dev/null 2>&1 || true' EXIT

git worktree add --detach "$WORKTREE_DIR" "$MERGE_BASE" >/dev/null

RESULT_ROWS=""
TRACKED=0
for module in $CHANGED_MODULES; do
  target="$(map_module_to_unit_test "$module")"
  if [[ -z "$target" ]]; then
    continue
  fi

  TRACKED=$((TRACKED + 1))
  base_cov="$(line_coverage_for "$WORKTREE_DIR" "$module" "$target")"
  head_cov="$(line_coverage_for "$(pwd)" "$module" "$target")"

  if [[ "$base_cov" == "NA" || "$head_cov" == "NA" ]]; then
    delta="NA"
  else
    delta="$(awk -v h="$head_cov" -v b="$base_cov" 'BEGIN{printf "%.2f", h-b}')"
  fi

  RESULT_ROWS+="| ${module} | ${base_cov}% | ${head_cov}% | ${delta} |"$'\n'
done

if [[ $TRACKED -eq 0 ]]; then
  echo "Changed modules do not yet have unit-test coverage targets; no per-module deltas to report."
  exit 0
fi

SUMMARY_FILE="${GITHUB_STEP_SUMMARY:-}"
if [[ -n "$SUMMARY_FILE" ]]; then
  {
    echo "### Coverage delta for touched modules"
    echo
    echo "Compared merge-base \\`${MERGE_BASE:0:12}\\` -> \\`HEAD\\`"
    echo
    echo "| module | base line coverage | PR line coverage | delta (pp) |"
    echo "| --- | ---: | ---: | ---: |"
    printf "%s" "$RESULT_ROWS"
  } >> "$SUMMARY_FILE"
fi

echo "Coverage delta report generated for $TRACKED touched module(s)."
