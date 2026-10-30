#!/usr/bin/env bash
# Run the agentops launcher with the repo's virtualenv when there is one (Windows or POSIX
# layout), else the first python on PATH. Used by .claude/settings.json hooks and status line.
root="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
for py in "$root/.venv/Scripts/python.exe" "$root/.venv/bin/python" python3 python; do
  if command -v "$py" >/dev/null 2>&1; then
    exec "$py" "$root/tools/agentops_cli.py" --config "$root/agentops.toml" "$@"
  fi
done
echo "agentops: no python found" >&2
exit 0
