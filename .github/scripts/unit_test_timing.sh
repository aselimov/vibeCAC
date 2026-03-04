#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -gt 0 ]; then
  TARGETS=("$@")
else
  TARGETS=(math str)
fi

ARTIFACT_DIR="src/artifacts"
RAW_REPORT="${ARTIFACT_DIR}/unit-test-timing.tsv"
MD_REPORT="${ARTIFACT_DIR}/unit-test-timing.md"

mkdir -p "${ARTIFACT_DIR}"
: > "${RAW_REPORT}"

(
  cd src
  make clean
)

warmup_target="${TARGETS[0]}"
echo "Preparing shared build artifacts with warm-up target: ${warmup_target}"
(
  cd src
  make unit-test UNIT_TEST="${warmup_target}"
)

for target in "${TARGETS[@]}"; do
  echo "Running unit test target: ${target}"
  start_epoch="$(date +%s)"
  (
    cd src
    make unit-test UNIT_TEST="${target}"
  )
  end_epoch="$(date +%s)"
  duration_sec="$((end_epoch - start_epoch))"
  printf "%s\t%s\n" "${target}" "${duration_sec}" >> "${RAW_REPORT}"
done

sorted_durations="$(awk '{print $2}' "${RAW_REPORT}" | sort -n)"
count="$(printf "%s\n" "${sorted_durations}" | awk 'NF{n++} END{print n+0}')"

if [ "${count}" -eq 0 ]; then
  echo "No timing data collected." >&2
  exit 1
fi

if [ $((count % 2)) -eq 1 ]; then
  middle_index="$((count / 2 + 1))"
  median="$(printf "%s\n" "${sorted_durations}" | awk -v idx="${middle_index}" 'NR==idx { print; exit }')"
else
  lower_index="$((count / 2))"
  upper_index="$((lower_index + 1))"
  lower="$(printf "%s\n" "${sorted_durations}" | awk -v idx="${lower_index}" 'NR==idx { print; exit }')"
  upper="$(printf "%s\n" "${sorted_durations}" | awk -v idx="${upper_index}" 'NR==idx { print; exit }')"
  median="$(awk -v a="${lower}" -v b="${upper}" 'BEGIN { printf "%.3f", (a + b) / 2.0 }')"
fi

outlier_threshold="$(awk -v m="${median}" 'BEGIN { t = 1.5 * m; if (t < 1.0) t = 1.0; printf "%.3f", t }')"

{
  echo "# Unit Test Timing Report"
  echo
  echo "| Target | Duration (s) | Relative to median | Outlier |"
  echo "| --- | ---: | ---: | --- |"
  while IFS=$'\t' read -r target duration; do
    ratio="$(awk -v d="${duration}" -v m="${median}" 'BEGIN { if (m == 0) printf "n/a"; else printf "%.2fx", d / m }')"
    is_outlier="no"
    if awk -v d="${duration}" -v t="${outlier_threshold}" 'BEGIN { exit !(d > t) }'; then
      is_outlier="yes"
    fi
    printf "| %s | %s | %s | %s |\n" "${target}" "${duration}" "${ratio}" "${is_outlier}"
  done < "${RAW_REPORT}"
  echo
  printf "Median: %s s; outlier threshold: > %s s (1.5x median, floor 1s).\n" "${median}" "${outlier_threshold}"
} > "${MD_REPORT}"

cat "${MD_REPORT}"

if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
  {
    echo "## Unit Test Timing"
    cat "${MD_REPORT}"
  } >> "${GITHUB_STEP_SUMMARY}"
fi
