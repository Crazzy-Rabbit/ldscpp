#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="$ROOT/test/wiki_flow/run"

if [[ -f /home/shilulu/anaconda3/etc/profile.d/conda.sh ]]; then
  set +u
  # shellcheck disable=SC1091
  source /home/shilulu/anaconda3/etc/profile.d/conda.sh
  conda activate ldsc_py2
  set -u
fi

cmake -S "$ROOT" -B "$ROOT/build"
cmake --build "$ROOT/build" -j2

python3 "$ROOT/test/generate_wiki_example.py"

run_commands() {
  local side="$1"
  local command_file="$ROOT/test/wiki_flow/commands_${side}.sh"
  (
    cd "$RUN_DIR/$side"
    while IFS= read -r cmd; do
      [[ -z "$cmd" ]] && continue
      printf '[%s] %s\n' "$side" "$cmd" | tee -a commands.executed.txt
      eval "$cmd"
    done < "$command_file"
  )
}

run_commands py
run_commands cpp

python3 "$ROOT/test/compare_wiki_flow.py" --run-dir "$RUN_DIR"
