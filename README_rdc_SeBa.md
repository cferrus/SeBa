# Reducing SeBa data: rdc_SeBa

> ## ⚠️ Not compatible with this fork's output
>
> `rdc_SeBa` **cannot read the 9-column `SeBa.data` produced by this fork** and is
> not used in this workflow. Its reader (`rdc/rdc_SeBa.C`, `read_SeBa_hist`) still
> expects upstream's **18-field** format (mass-transfer type, stellar radii,
> effective temperatures, core masses). Fed this fork's
> [9-column output](README.md#output-format), every field shifts by one — `time`
> is read into the mass-transfer slot, `semi` into `time`, and so on — so it
> emits **misaligned garbage without printing any error**.
>
> **Do not pipe this fork's `SeBa.data` through `rdc_SeBa`.** Filter it yourself
> instead — it is just 9 whitespace-separated columns with exactly two lines per
> binary (initial state at T=0, final state at T=end), which is trivial to slice
> in Python/pandas or awk. See the [filtering example](#filtering-without-rdc_seba)
> below.
>
> Fixing `rdc_SeBa` would require rewriting `read_SeBa_hist` (and the matching
> `operator<<` / `put_single_reverse` writers) for the 9-column layout. This has
> not been done because the tool is not part of this workflow. The original usage
> is preserved below for reference only.

---

## Filtering without rdc_SeBa

Because the output is 9 columns and exactly two lines per binary, common filters
are a few lines of Python. Column order (see [`README.md`](README.md#output-format)):

```
identity  bin_type  time  semi  ecc  prim_type  prim_mass  sec_type  sec_mass
```

Example — pull the final state of all double white dwarf (types 12–14) systems:

```python
WD = {12, 13, 14}
finals = {}
for line in open("SeBa.data"):
    c = line.split()
    bid, t = c[0], float(c[2])
    if bid not in finals or t >= float(finals[bid][2]):
        finals[bid] = c
for c in finals.values():
    if int(c[5]) in WD and int(c[7]) in WD:
        print(" ".join(c))
```

Swap the `WD` set / column indices for other selections (binary type in col 2,
masses in cols 7 and 9, etc.). Binary-type integers and stellar-type integers are
tabulated in [`README.md`](README.md#output-format).

Note this selects on stellar type alone, so it also catches *disrupted* pairs
(binary type 7) that are no longer bound. To keep only bound systems, additionally
require the binary type — e.g. add `int(c[1]) == 2 and` (2 = Detached) to the
`if`.

---

## Original rdc_SeBa usage (reference only — does not work with this fork)

`rdc_SeBa` filtered upstream's full event-by-event `SeBa.data`, picking out the
binary systems of interest. It is still built by the top-level `make`, but see the
warning above before trying to use it.

```bash
less SeBa.data | ./rdc_SeBa -f -p white_dwarf -s white_dwarf > SeBa_wdwd.data
```

Filter parameters:

```
-f  first occasion — first matching line for a given binary
-R  full evolution — every matching line
-p  stellar type of primary star
-s  stellar type of secondary star
-B  binary type
```

Stellar type strings: `any, proto_star, planet, brown_dwarf, main_sequence,
hertzsprung_gap, sub_giant, horizontal_branch, super_giant, helium_star,
helium_giant, white_dwarf, neutron_star, black_hole` (short forms `ps, pl, bd, ms,
hg, gs, hb, sg, he, gh, hd, cd, od, ns, bh`).

Binary type strings (from `extract_binary_type_string`): `synchronized,
strong_encounter, detached, semi_detached, contact, common_envelope,
double_spiral_in, merged, disrupted, spiral_in`.

Other range constraints: `-a/-A` separation, `-m/-M` primary mass, `-n/-N`
secondary mass, `-q/-Q` mass ratio, `-e/-E` eccentricity.
</content>
