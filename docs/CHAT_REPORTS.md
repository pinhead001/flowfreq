# Reports through a chat

`flowfreq-mcp` is flowfreq's report tools as a [Model Context Protocol](https://modelcontextprotocol.io)
server. A chat client (Claude Code, Claude Desktop) starts it and calls its tools, and a
report takes shape over the conversation: find a gage, fit it, estimate at a site, weight
the two, write a few paragraphs, render the Word document.

The division of labor is the one `flowfreq.reporting` already has. Every table and figure
is computed by flowfreq from inputs stored in the report's specification. The model
contributes prose only, and that prose prints under "Draft narrative (written by a language
model; not reviewed)" until you say you have read it.

## Install

```bash
pip install "flowfreq[mcp]"      # the mcp package (2.3 or later, 2.x) and python-docx
flowfreq-mcp --version
```

The server speaks MCP on standard input and output, so it is started by the client and not
by hand. Run directly, it waits silently for a client.

## Connect a client

Give the client the full path of `flowfreq-mcp` inside the environment flowfreq is installed
in. A bare `flowfreq-mcp` works only if that environment is on the client's `PATH`.

**Claude Code**

```bash
claude mcp add flowfreq --env USGS_API_KEY=your-key -- D:\a\flowfreq\.venv\Scripts\flowfreq-mcp.exe --dir D:\reports
```

Add `--scope user` to have it in every project. `/mcp` inside Claude Code lists the server
and its sixteen tools.

**Claude Desktop** — in `claude_desktop_config.json` (Settings, Developer, Edit Config):

```json
{
  "mcpServers": {
    "flowfreq": {
      "command": "D:\\a\\flowfreq\\.venv\\Scripts\\flowfreq-mcp.exe",
      "args": ["--dir", "D:\\reports"],
      "env": {"USGS_API_KEY": "your-key"}
    }
  }
}
```

Restart the client after editing the file.

`--dir` is where reports are written and where saved reports are read from. Without it the
server uses `$FLOWFREQ_REPORT_DIR`, or `flowfreq-reports` in your home directory. Every
file name a tool takes is a plain name inside that directory; a name with a directory in it
is refused. `USGS_API_KEY` is optional (see the README): without one the USGS Water Data
API allows 1000 requests an hour.

## A conversation

> Find the Big Sandy River gage at Bruceton, Tennessee and run a flood-frequency analysis.
> Then estimate the regression at the gage from a drainage area of 205 square miles in
> region "HA4 DA only", weight the two, add a short summary, and give me the report as
> `big_sandy.docx`.

The model calls `find_gage`, `run_ffa`, `estimate_from_characteristics`,
`weight_gage_with_regression`, `get_section` (to read the numbers it will write about),
`edit_narrative` and `render_report`. The result is `big_sandy.docx` and `big_sandy.json`
in the output directory. The JSON is the specification as built, with the downloaded peaks
and the looked-up skew in it, so `flowfreq report big_sandy.json -o big_sandy.docx`
rebuilds the same report with no network and no chat.

## Tools

| Tool | What it does | Needs |
|---|---|---|
| `find_gage` | Search the shipped gage catalog by station number or words from the name | nothing |
| `run_ffa` | Add a Bulletin 17C section for a gage: record, method, parameters, quantiles with 90-percent confidence limits, frequency curve, stationarity screen, peaks | USGS (peaks; site location for the regional skew) |
| `estimate_ungaged` | Add a regional regression estimate at a point, basin delineated by StreamStats | StreamStats |
| `estimate_from_characteristics` | The same estimate from basin characteristics you supply | nothing |
| `weight_gage_with_regression` | Weight a gage section with a regression section by inverse variance (Bulletin 17C Appendix 9) | nothing |
| `add_text_section` | Add a section of prose only | nothing |
| `edit_narrative` | Set the prose after a section's computed content | nothing |
| `mark_narrative_reviewed` | Remove the draft label from one section's narrative | your say-so |
| `get_section` | A section as it will print, in Markdown | nothing |
| `report_outline` | Title block and the sections, with ids and narrative status | nothing |
| `set_report_details` | Title, subtitle, author, project, date | nothing |
| `remove_section`, `move_section` | Take a section out, or move it | nothing |
| `render_report` | Write `.docx` or `.md`, and the `.json` specification beside it | nothing |
| `save_report`, `load_report` | Carry a report between conversations as its `.json` specification | nothing |

Sections have ids (`s1`, `s2`, ...) that do not change when sections move and are never
reused. By default `weight_gage_with_regression` replaces the two sections it weights,
since the weighted section carries both analyses; `keep_parts` keeps them, and a part that
has a narrative is never replaced.

## What the server holds the model to

- **Numbers come from tools.** The server's instructions tell the model not to state a
  flow, a probability or a statistic a tool did not return, and not to compute one. That
  is an instruction and not an enforcement: read the narrative against the tables before
  marking it reviewed. The tables and figures themselves cannot be written by the model.
- **Narrative is paragraphs and bullet lists.** One with a table or a heading is refused,
  so a table in a report always came from flowfreq.
- **Narrative is a draft until you review it.** Editing a reviewed narrative makes it a
  draft again. `render_report` returns how many drafts remain.
- **Refusals are reported, not worked around.** A regulated record, a characteristic
  outside an equation's range, or no regional skew on file comes back as an error naming
  the reason. The `allow_*` arguments override a refusal, and the model is told to set
  one only when you ask.

One report is held for the life of the server process, which is one client session. Use
`save_report` (or `render_report`, which also writes the specification) before closing the
client, and `load_report` to pick the report up again. Tool calls are run one at a time,
in the order they arrive.

## From Python

`flowfreq.reporting.ReportSession` is the object behind the tools, usable without MCP:

```python
from flowfreq.reporting import ReportSession

report = ReportSession("Big Sandy River flood frequency", output_dir="reports")
gage = report.add_section("gaged_ffa", {"site_no": "03606500", "regional_skew_from_site": True})
print(report.section_markdown(gage["section_id"]))
report.set_narrative(gage["section_id"], "The fit weights the station skew with the regional skew.")
report.render("big_sandy.docx")
```

`flowfreq.mcp_server.build_server(session)` wraps a session in the tools, which is how the
tests drive it in process.

## Limits

- `estimate_ungaged` depends on StreamStats, which takes a minute or more and is sometimes
  unavailable. `estimate_from_characteristics` is the route that needs no service.
- There is no PDF writer; export one from Word.
- The section kinds are those of `flowfreq.reporting.SECTION_KINDS`. Transposition,
  low-flow and regime analyses are not report sections yet.
