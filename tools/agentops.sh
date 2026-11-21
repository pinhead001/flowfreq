#!/usr/bin/env bash
# Run the shared agentops package (D:/a/agentops, pip install -e into .venv) on config/agents.toml
# with the repo's virtualenv when there is one (Windows or POSIX layout), else the first python on
# PATH. Used by the hooks and status line in .claude/settings.local.json; run.bat is the same for
# a person. Exits 0 when no python is found, so a missing install never blocks a session.
root="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
for py in "$root/.venv/Scripts/python.exe" "$root/.venv/bin/python" python3 python; do
  if command -v "$py" >/dev/null 2>&1; then
    exec "$py" -m agentops.cli --config "$root/config/agents.toml" "$@"
  fi
done
echo "agentops: no python found" >&2
exit 0
