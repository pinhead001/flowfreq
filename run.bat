@echo off
rem flowfreq entry points. Usage: run.bat <command> [args]
rem   gov status|next|plan|calibrate --agy-remaining PCT   governor (agentops, D:\a\agentops)
rem   cards add|list|run|move|reset ...                    agy task-card queue (config\cards\README.md)
rem   ledger [--tail N]                                     delegate call ledger
rem   statusline | hook <event>                             used by .claude\settings.json
rem   accept <pytest paths...>                              lint + mypy + pytest in THIS checkout (worktree-safe)
rem   check                                                 make check with the venv on PATH (what CI runs)
rem agentops is run from its own venv (uv run --project D:\a\agentops); flowfreq stays a pip project.
setlocal
set "HERE=%~dp0"
set "HERE=%HERE:~0,-1%"
set "AGENTOPS=D:\a\agentops"
set "AGENTOPS_CONFIG=%HERE%\config\agents.toml"
rem The main checkout's interpreter: a card worktree has no venv and no editable install of its own.
set "PY=%HERE%\.venv\Scripts\python.exe"
if not exist "%PY%" set "PY=D:\a\flowfreq\.venv\Scripts\python.exe"
if "%~1"=="" goto :help
if /i "%~1"=="gov" goto :agentops
if /i "%~1"=="cards" goto :agentops
if /i "%~1"=="ledger" goto :agentops
if /i "%~1"=="statusline" goto :agentops
if /i "%~1"=="hook" goto :agentops
if /i "%~1"=="accept" goto :accept
if /i "%~1"=="check" goto :check
goto :help

:agentops
uv run --project "%AGENTOPS%" agentops --config "%AGENTOPS_CONFIG%" %*
exit /b %ERRORLEVEL%

:accept
rem PYTHONPATH=%HERE% makes `import flowfreq` resolve to this checkout, not the main checkout's
rem editable install; PYTHONSAFEPATH matches the `pytest` console script CI runs; PYTHONUTF8 keeps
rem isort from skipping the non-ASCII files (CLAUDE.md, Build & Development Commands).
shift
set "PYTHONPATH=%HERE%"
set "PYTHONSAFEPATH=1"
set "PYTHONUTF8=1"
"%PY%" -c "import flowfreq, os, sys; p = os.path.normcase(os.path.abspath(flowfreq.__file__)); print(p); sys.exit(0 if p.startswith(os.path.normcase(os.path.abspath(os.environ['HERE']))) else 3)" || exit /b 3
"%PY%" -m black --check flowfreq tests || exit /b 1
"%PY%" -m isort --check-only flowfreq tests || exit /b 1
"%PY%" -m mypy flowfreq/ || exit /b 1
"%PY%" -m pytest -q -x %1 %2 %3 %4 %5 %6 %7 %8 %9
exit /b %ERRORLEVEL%

:check
for %%I in ("%PY%") do set "PATH=%%~dpI;%PATH%"
set "PYTHONPATH=%HERE%"
make check
exit /b %ERRORLEVEL%

:help
findstr /b /c:"rem " "%~f0"
exit /b 0
