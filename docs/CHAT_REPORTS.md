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
and its eighteen tools.

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
| `transpose_gage_to_site` | Carry a gage section's quantiles to an ungaged site on the same stream by the drainage-area ratio, with exponents you supply and cite | nothing |
| `check_gage_against_reference` | Compare a gage section's fit with USGS peakfq 8.1.0 on the same record, or with the numbers of another program's run that you supply | nothing |
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

`transpose_gage_to_site` leaves the gage's section where it is and puts the transposed one
after it. The drainage-area exponents are yours to give, from the regional regression report
for the area, with the report and table as `exponent_citation`: either one `exponent` for
every AEP, or `exponent_aeps` and `exponents` as the report's table has them. flowfreq ships
no table of exponents, and the server's instructions tell the model to ask for them rather
than recall them. An example request:

> Transpose the gage's quantiles to the bridge site, 150 square miles, with exponents 0.70,
> 0.66 and 0.62 at AEP 0.5, 0.1 and 0.01 from table 4 of the regional report (give the full
> citation here).

The section's table shows, for each AEP, the gage's flow, the exponent, the factor applied
and whether that exponent was published, interpolated between published AEPs or held beyond
them. An area ratio outside 0.5 to 1.5 is refused unless `allow_out_of_range`, and a
regulated site on one side only is refused unless `allow_regulation_mismatch`; either
override is written into the section's caveats.

`check_gage_against_reference` adds a section that says whether the gage's fit agrees with
a second computation of the same record, and tabulates the two side by side: fitted
parameters, quantiles and 90-percent confidence limits, each difference against its
tolerance (1 percent for quantiles and 2 percent for confidence limits unless you set
others). The gage's section stays. There are two sources for the second computation:

- **USGS peakfq 8.1.0**, when no `reference_*` argument is given. flowfreq sends the gage's
  peaks, thresholds and regional skew to peakfq's Fortran routine through its own
  extension, which has to be built on the machine the server runs on
  (`python build_fortran/build.py`, needing gfortran and meson; `docs/FORTRAN_UPLOAD.md`).
  Without it the tool refuses and says so. What peakfq returned is stored in the section,
  with the date, so the report rebuilds later on a machine without the extension. If the
  gage's inputs change afterwards the stored numbers are refused as belonging to another
  record.
- **Numbers you give**, from your own run of peakfq, HEC-SSP or another program:
  `reference_source` (the program, its version and where the numbers were taken from),
  `reference_aeps` with `reference_flows_cfs`, and optionally `reference_lower_cfs` and
  `reference_upper_cfs`. Only what you give is compared. The model is told to ask for
  these and not to supply them from memory.

Agreement here means the two programs compute the same thing from the same inputs. It
says nothing about whether the inputs suit the gage.

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
- The section kinds are those of `flowfreq.reporting.SECTION_KINDS`. Low-flow and regime
  analyses are not report sections yet, and neither are the other transfers in
  `flowfreq.transpose` (flow-duration, low-flow, between two gages, regression-weighted).
