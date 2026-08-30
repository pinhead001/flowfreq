# Vendored: peakfq (USGS)

Reference copy of the USGS `peakfq` R package, vendored so the Fortran EMA
implementation can be read, rebuilt, and compared against from this repository.

| | |
|---|---|
| Upstream package | `peakfq` |
| Version | 8.1.0 |
| License | CC0 (public domain dedication) — see `LICENSE.md` |
| Authors | Seth Siefken, Sophia Crouch, Shanna Blount (USGS) |
| Retrieved | 2026-08-30 |
| Reference | Bulletin 17C, doi.org/10.3133/tm4B5 |

Naming note: upstream the package is `peakfq`; this directory is `peakfqr` to match
the local reference tree and every existing reference in `TODO.md` and
`tests/peakfqsa/fixtures/`.

## What is here

- `src/` — the authoritative Fortran. `emafit.f` carries `emafitpr`, and with it
  `var_mom`, `EXPMOMCDERIV`, `DEXPECT`, `as_G_PRL_o` and the Inverse Modified
  Cholesky Gaussian Quadrature — none of which `hydrolib` implements yet, and which
  are what the open confidence-interval question turns on.
- `R/` — call conventions. `fortranWrappers.R` is the specification that
  `build_fortran/_emafort.pyf` was written against.
- `inst/testdata/` — byte-exact expected output, marked `-text` in `.gitattributes`.
- `tests/testthat/` — the upstream tests the fixtures in `tests/peakfqsa/fixtures/`
  were transcribed from.

Do not edit anything in this directory. It is a verbatim reference copy, and a change
here would silently invalidate every comparison made against it.
