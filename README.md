# Catchment-scale conservative solute transport: paper configurations

This branch contains input configurations only: mesh data, global and solver
settings, boundary concentration histories, channel paths and transport settings.
It contains no solver source code, executable, simulation results or NetCDF payloads.

## Source code and revisions

DRUtES source code is maintained at
https://github.com/michalkuraz/drutes-dev/tree/testing . The source snapshot
verified when publishing these configurations is
`2e0f1de7151dcebc553f0f373bcfd20847470f53`.
The `testing` branch can change; record a commit when reproducing a calculation.

This repository's `adencconfigs` branch contains two input collections:

- `original-inputs/`: nine archived input sets recovered from the actual paper
  simulations, retaining their original small-file contents and old input format.
- `current-inputs/`: six record-preserving conversions of Cases II and III to
  the consolidated `transport.conf` and `stabilization.conf` input format read by
  the source snapshot above. No physical or numerical parameter was changed.

Case I is deliberately not supplied as a current-version rerun: its published
results used contact-limited widths and channel-path directions. Later code uses
full transverse cell widths and mRM routing directions. Changing only input
filenames would not reproduce the historical Case I calculation.

For Case I, use the archived inputs with historical source commit
`810e5c9528ce72ad62c1531cec4ba962c63935f8`, available in the source repository's
`testing` history. Those runs predate mandatory restart routing. The restart
symlink supplied in these input sets is for data completeness and is not evidence
that the historical calculation used restart-derived directions.

Cases II and III originally used a private Linux build based on commit
`b83950795aa0c6fe7c4cd4f9df8bbf70ca6fdb75` plus the initial unported-component
filter. Public commit `25988326c1582e5a618404bed8a27ca6ff2ed2e3` contains that
filter and reads the original five-file format. Current commit `2e0f1de7`
includes input consolidation and coding-style changes. Previously documented
paired tests of the input-format conversion found identical saved fields in
short comparisons; publishing these files does not constitute a new full-duration
rerun or a guarantee of bitwise identity across compilers and computers.
All six current-format configurations also passed fresh initialization, routing
and conservative matrix-assembly checks during this publication step.

## Experiments

| Paper case | Input folder prefix | Start | Duration | Contaminated inlet | Release | Qmin |
|---|---|---|---|---|---|---|
| I | `case-I-october-rhine-pulse` | 2015-10-24 | 14 days | Rhine (101) | C=1 for 6 h, then 0 | 300 m3/s |
| II | `case-II-june-moselle-continuous` | 2015-06-01 | 100 days | Moselle (102) | Continuous C=1 | 40 m3/s |
| III | `case-III-june-moselle-pulse` | 2015-06-01 | 100 days | Moselle (102) | C=1 for 2 days, then 0 | 40 m3/s |

Every case has `galerkin`, `supg2`, and `supg2-shock1` variants. The other inlet
in Cases II and III supplies clean water. These are idealized releases with no
assigned environmental concentration scale, not equal-injected-mass experiments.
All cases use longitudinal dispersivity 200 m and transverse dispersivity 0.2 m.
Case I has four observation points; Cases II and III have seven. The 100-day
calculations retain dense observation times during the first two weeks.

## Download the shared NetCDF inputs

Download the three files from **https://drutes.org/~miguel/adencconfigs/**
into `shared-data/` at the root of this checkout:

```sh
git clone --branch adencconfigs --single-branch https://github.com/michalkuraz/drutes_configs.git
cd drutes_configs
curl -fL --retry 3 -C - -o shared-data/mRM_Fluxes_States.nc https://drutes.org/~miguel/adencconfigs/mRM_Fluxes_States.nc
curl -fL --retry 3 -C - -o shared-data/mRM_restart_001.nc https://drutes.org/~miguel/adencconfigs/mRM_restart_001.nc
curl -fL --retry 3 -C - -o shared-data/dem.nc https://drutes.org/~miguel/adencconfigs/dem.nc
cd shared-data
shasum -a 256 -c SHA256SUMS
cd ..
```

The public URL contains `~miguel` (a tilde). The path `/miguel/adencconfigs/`
without the tilde is not the same URL.

| File | Size (bytes) | Role |
|---|---:|---|
| `mRM_Fluxes_States.nc` | 281014845 | Time-dependent routed discharge used by transport |
| `mRM_restart_001.nc` | 488987890 | Original full restart input used in the June simulations; the transport reader uses its static routing fields only |
| `dem.nc` | 16944 | Shared terrain input |

Each configuration contains relative symlinks to this single shared copy.
Keep the checkout's folder structure and use a filesystem that supports symlinks.
Git tracks the links, not the downloaded NetCDF contents. Do not replace the
shared files with different datasets without rechecking the grid, routing and
hydraulic window. These are inputs to the transport model, not the complete
meteorological forcing and calibration package needed to reproduce mHM itself.

## Selecting and running a configuration

Build the appropriate source revision in the separate `drutes-dev` repository.
For a current-format example, change to a case directory that contains
`drutes.conf`, create a new empty `out` directory, and run the external executable:

```sh
cd current-inputs/case-III-june-moselle-pulse-supg2-shock1
test ! -e out && mkdir out && /absolute/path/to/drutes-dev/bin/drutes
```

Replace the executable path with your own build. For historical Case I, use
`original-inputs/case-I-october-rhine-pulse-<variant>` and source commit `810e5c9`.
The executable clears `out/*` at startup: preserve previous results before a rerun.
No simulations were launched while assembling this repository.

Current-format stabilization records are, in order: SUPG flag, SUPG multiplier,
residual-viscosity flag, residual-viscosity multiplier. The three variants are
`n,2,n,1`; `y,2,n,1`; and `y,2,y,1`, each value on its own line. The separate
historical files are retained only under `original-inputs/`; current source
rejects them, so do not copy them into a current-format configuration.

`INPUT_SHA256SUMS` covers all regular published configuration inputs. To verify:

```sh
shasum -a 256 -c INPUT_SHA256SUMS
```

The shared NetCDF inputs have their separate manifest in `shared-data/`.
The full published restart and the earlier compact routing extract have identical
values for all ten variables retained by the extract; the full archived input is
provided here to preserve the original data file, not to restore mHM water states.
Bank conditions, forcing, mesh and inlet concentration histories remain the
archived settings. Conservation diagnostics do not establish concentration
positivity, hydrological calibration, or accurate field-scale travel times.
