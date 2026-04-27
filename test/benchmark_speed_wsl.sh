#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="$ROOT/test/full_flow/run"
SPEED_DIR="$ROOT/test/speed"
REPEATS="${1:-5}"

if [[ -f /home/shilulu/anaconda3/etc/profile.d/conda.sh ]]; then
  set +u
  # shellcheck disable=SC1091
  source /home/shilulu/anaconda3/etc/profile.d/conda.sh
  conda activate ldsc_py2
  set -u
fi

cmake -S "$ROOT" -B "$ROOT/build" >/dev/null
cmake --build "$ROOT/build" -j2 >/dev/null
python3 "$ROOT/test/generate_full_parity.py"

rm -rf "$SPEED_DIR"
mkdir -p "$SPEED_DIR/logs"
RESULTS="$SPEED_DIR/results.tsv"
SUMMARY="$SPEED_DIR/summary.tsv"
printf 'scope\tside\titeration\tcommand_index\tseconds\tcommand\n' > "$RESULTS"

elapsed_seconds() {
  local start_ns="$1"
  local end_ns="$2"
  awk -v s="$start_ns" -v e="$end_ns" 'BEGIN { printf "%.6f", (e - s) / 1000000000.0 }'
}

run_side() {
  local side="$1"
  local iter="$2"
  local command_file="$ROOT/test/full_flow/commands_${side}.sh"
  local side_dir="$RUN_DIR/$side"
  local log_file="$SPEED_DIR/logs/${side}_${iter}.log"

  rm -rf "$side_dir/out"
  mkdir -p "$side_dir/out"

  local total_start total_end
  total_start="$(date +%s%N)"
  local idx=0
  (
    cd "$side_dir"
    while IFS= read -r cmd; do
      [[ -z "$cmd" ]] && continue
      idx=$((idx + 1))
      local start end seconds
      start="$(date +%s%N)"
      eval "$cmd" >> "$log_file" 2>&1
      end="$(date +%s%N)"
      seconds="$(elapsed_seconds "$start" "$end")"
      printf 'command\t%s\t%s\t%s\t%s\t%s\n' "$side" "$iter" "$idx" "$seconds" "$cmd" >> "$RESULTS"
    done < "$command_file"
  )
  total_end="$(date +%s%N)"
  printf 'total\t%s\t%s\tall\t%s\tfull_flow\n' "$side" "$iter" "$(elapsed_seconds "$total_start" "$total_end")" >> "$RESULTS"
}

for iter in $(seq 1 "$REPEATS"); do
  run_side py "$iter"
  run_side cpp "$iter"
done

python3 - "$RESULTS" "$SUMMARY" <<'PY'
import csv
import statistics
import sys
from collections import defaultdict

results_path, summary_path = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(results_path, encoding="utf-8"), delimiter="\t"))
groups = defaultdict(list)
for row in rows:
    key = (row["scope"], row["command_index"], row["side"])
    groups[key].append(float(row["seconds"]))

summary = []
keys = sorted({(row["scope"], row["command_index"]) for row in rows}, key=lambda x: (x[0] != "total", x[1]))
for scope, command_index in keys:
    py = groups.get((scope, command_index, "py"), [])
    cpp = groups.get((scope, command_index, "cpp"), [])
    if not py or not cpp:
        continue
    py_med = statistics.median(py)
    cpp_med = statistics.median(cpp)
    speedup = py_med / cpp_med if cpp_med else float("inf")
    reduction = (1.0 - cpp_med / py_med) * 100.0 if py_med else 0.0
    command = next(row["command"] for row in rows if row["scope"] == scope and row["command_index"] == command_index and row["side"] == "py")
    summary.append({
        "scope": scope,
        "command_index": command_index,
        "py_median_s": f"{py_med:.6f}",
        "cpp_median_s": f"{cpp_med:.6f}",
        "speedup_x": f"{speedup:.3f}",
        "time_reduction_pct": f"{reduction:.1f}",
        "command": command,
    })

with open(summary_path, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=["scope", "command_index", "py_median_s", "cpp_median_s", "speedup_x", "time_reduction_pct", "command"],
    )
    writer.writeheader()
    writer.writerows(summary)

for row in summary:
    print("\t".join(row[k] for k in ["scope", "command_index", "py_median_s", "cpp_median_s", "speedup_x", "time_reduction_pct", "command"]))
PY
