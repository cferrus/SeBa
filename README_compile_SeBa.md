# Installing & Compiling SeBa (`cferrus` fork)

This is a modified fork of [`amusecode/SeBa`](https://github.com/amusecode/SeBa).
See [`README.md`](README.md) for how it differs (9-column output, silent-by-default,
OpenMP) and [`CLAUDE.md`](CLAUDE.md) for build internals.

## Installing

Clone this fork:

```bash
git clone https://github.com/cferrus/SeBa.git
cd SeBa
```

The only prerequisite is a C++ compiler and `make`. There are no external
dependencies.

## Compiling

From the repo root:

```bash
make clean
make            # builds libsstar.a, libnode.a, libstd.a, librdc.a
cd dstar
make            # builds libdstar.a and the SeBa executable
```

This creates the executable at `dstar/SeBa`. You can rename it (e.g. to tag a
version) or move it anywhere. When you run it, an output file `SeBa.data` is
created in the current directory.

**Repeat these steps after changing the code.**

### OpenMP (automatic)

- **Linux / GCC**: `dstar/Makefile` detects `uname -s == Linux` and adds
  `-fopenmp`, giving multi-core parallelism automatically. Verify with:
  ```bash
  ldd dstar/SeBa | grep omp     # should show libgomp or similar
  ```
- **macOS / Apple Clang**: OpenMP is unavailable without `libomp`; the build
  falls back to correct single-threaded execution automatically.
- **macOS / Homebrew GCC**: force it with `CXX=g++-14 make`.

### Relinking after library changes

`dstar/SeBa` only relinks when `SeBa.C` itself is newer than the binary. After
editing any library file (e.g. `double_star.C`, `single_star.C`), the executable
will **not** rebuild on its own. Force a relink:

```bash
rm -f dstar/SeBa && cd dstar && make SeBa
```

### Makefile.inc regeneration

`Makefile.inc` is auto-generated from `Makefile.inc.conf` by substituting the
absolute repo path, and is not committed. Delete it to force regeneration:

```bash
rm -f Makefile.inc dstar/Makefile.inc && make clean && make
```

## Smoke test

```bash
cd dstar
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
wc -l SeBa.data   # must be 2 (initial + final state)
```
</content>
