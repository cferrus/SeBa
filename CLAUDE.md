# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

SeBa requires only a C++ compiler. Build in two stages — the top-level step builds the shared libraries, then `dstar/` builds the executable:

```bash
make clean && make   # builds libsstar.a, libnode.a, libstd.a, librdc.a in their subdirs
cd dstar && make     # builds libdstar.a and the SeBa executable
```

The executable lands at `dstar/SeBa`. `Makefile.inc` is auto-generated from `Makefile.inc.conf` by substituting the absolute path; delete it to force regeneration.

There are no automated tests. To verify a build, run a quick single-system evolution:

```bash
cd dstar
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

Output goes to `SeBa.data` in the current directory.

## Running SeBa

Three modes:

| Mode | Example |
|---|---|
| Single system | `./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001` |
| Input file | `./SeBa -I SeBa_input.txt` (format: `a e M m z` per line) |
| Random population | `./SeBa -R -n 250000 -T 13500 -f 4 -m 0.95 -M 10` |

## Filtering output with rdc_SeBa

`rdc/rdc_SeBa` filters `SeBa.data` by stellar type and binary state. Must be built as part of the top-level `make`:

```bash
less SeBa.data | ./rdc/rdc_SeBa -f -p white_dwarf -s white_dwarf > SeBa_wdwd.data
```

Key flags: `-f` (first occurrence only), `-R` (full evolution), `-p`/`-s` (primary/secondary type), `-B` (binary type).

## Architecture

The codebase is a C++ class hierarchy rooted at `starbase` → `star`:

- **`star`** (`include/star/star.h`): abstract base for all evolving objects; holds the `node*` tree pointer and unit-conversion factors.
- **`single_star`** (`include/star/single_star.h`): concrete stellar evolution. Each stellar phase is a derived class (`main_sequence`, `hertzsprung_gap`, `sub_giant`, `horizontal_branch`, `super_giant`, `helium_star`, `helium_giant`, `white_dwarf`, `neutron_star`, `black_hole`, etc.) living in `sstar/starclass/`. Phase transitions work by replacing the `single_star` object on the `node`.
- **`double_star`** (`include/star/double_star.h`, `dstar/starclass/double_star.C`): wraps a binary system. Owns the orbital parameters (`semi`, `eccentricity`, `bin_type`) and drives timestep-coupled evolution of both components. Mass-transfer stability, common-envelope events, and tidal circularisation are handled here.
- **`node`** / **`dyn`** (`node/`, `include/`): the N-body tree structure inherited from Starlab. In SeBa only the tree bookkeeping is used; direct gravitational dynamics are not exercised.

### Directory map

```
include/          Global headers (star/, node-level)
sstar/starclass/  Single-star phase implementations (one .C per stellar type)
dstar/starclass/  double_star and binary support
dstar/stardyn/    Binary-state evolution logic (mass transfer, CE)
dstar/init/       Initialisation of double systems
dstar/util/       dstar-level utilities (dump, semi↔period conversion)
dstar/io/         Binary I/O
dstar/SeBa.C      Main entry point; CLI parsing, population loop
rdc/              Post-processing tools (rdc_SeBa and friends)
std/              Shared low-level utilities
node/dyn/         Kepler orbit routines used by double_star
starev.C          Standalone single-star evolution tool (top-level)
```

### Stellar type enum

Defined in `include/star/stellar_type.h`. The integer codes used in `SeBa.data` output map to: 3 = Main_Sequence, 5 = Hertzsprung_Gap, 6 = Sub_Giant, 7 = Horizontal_Branch, 8 = Super_Giant, 10 = Helium_Star, 11 = Helium_Giant, 12–14 = White_Dwarf variants, 18 = Neutron_Star, 19 = Black_Hole.

### Output files

- `SeBa.data`: one line per interesting event; columns described in `README.md`.
- `init.dat`: initial conditions for each system.
- `binev.data`: remnant formation events.

### Key compile-time flags

- `-DTOOLBOX`: activates the `main()` in each `.C` file (allows per-file standalone builds).
- `-DHAVE_CONFIG_H`: enables `config.h` inclusion.
- Debug verbosity is controlled by `#define REPORT_*` booleans near the top of `double_star.C` and related files — flip them to `true` for detailed trace output without recompiling the full tree.
