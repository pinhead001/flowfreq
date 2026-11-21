@echo off
rem flowfreq agentops entry points: the shared package D:\a\agentops, pip-installed into .venv
rem (pip install -e D:\a\agentops; docs\AGENTOPS.md). Usage: run.bat <command> [args]
rem   gov status|next|plan|check|calibrate --agy-remaining PCT|profile [name|--clear]
rem   cards new ID --repo flowfreq --title T --out .agentops/cards | add FILE.. | list [--state S]
rem   cards run [--workers N] [--card ID] | show ID | diff ID [--stat]
rem   cards accept ID | reject ID --reason R | reset ID | move ID STATE
rem   usage [--json] | ledger [--tail N] [--card ID] [--kind agy|quota|review]
rem   statusline | hook <event>     used by .claude/settings.local.json
setlocal
set "HERE=%~dp0"
set "HERE=%HERE:~0,-1%"
set "AGENTOPS_CONFIG=%HERE%\config\agents.toml"
if "%~1"=="" goto :help
for %%V in (gov cards usage ledger statusline hook) do if /i "%~1"=="%%V" goto :agentops
goto :help
:agentops
"%HERE%\.venv\Scripts\python.exe" -m agentops.cli --config "%AGENTOPS_CONFIG%" %*
exit /b %ERRORLEVEL%
:help
findstr /b /c:"rem " "%~f0"
exit /b 0
