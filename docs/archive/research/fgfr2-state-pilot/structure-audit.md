> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# FGFR2 paired-conformation structure audit

## Recommendation

Use the published **unphosphorylated FGFR2 starting-coordinate pair** from the
2023 PNAS molecular-dynamics study. Both preparations contain the same engineered
C491A expression substitution. The active preparation derives from 2PVF; the
autoinhibited preparation is a Rosetta homology model built from FGFR1 3KY2.
Describe the comparison as **structure-informed prediction from prepared
conformation proxies**, not an experimental FGFR2 active/inactive pair, free-energy
calculation, or mutation-specific simulation.

This selection was made from structural provenance before fitting the activation
models. No assay outcomes, mutant MD archives, trajectories, or original final
labels were accessed in this audit.

## Why the obvious crystal pair is misleading

The original 2007 work and PDB titles present 2PSQ as unphosphorylated/inactive and
2PVF as phosphorylated/active. The 2017 mechanistic study subsequently demonstrated
that 2PSQ has an unusually open activation loop. Crystal contacts and a sulfate
ion favor features of the active conformation. NMR measurements support a more
closed activation loop in solution. The 2023 study explicitly used a homologous
FGFR1 template because an authentically autoinhibited FGFR2 crystal was unavailable.
[2007 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2094128/),
[2017 study](https://elifesciences.org/articles/21137),
[2023 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9974468/).

A 2024 computational structural survey labels 1GJO inactive and favors it as an
unphosphorylated reference. That is useful for cataloging structures, but does not
establish a genuinely autoinhibited activation loop. In our coordinate audit,
1GJO's Y656-to-I707 C-alpha distance is 30.9 Angstrom, compared with 13.9 Angstrom
in the published autoinhibited model. We therefore did not use a categorical
"inactive" label from a structural survey as sufficient mechanistic evidence.
[2024 survey](https://pmc.ncbi.nlm.nih.gov/articles/PMC11050683/).

## Experimental candidates inspected

All numbers below come from the downloaded RCSB entry JSON, polymer entity JSON,
and mmCIF coordinates. Counts are per selected author chain, counting residues
with C-alpha atoms; they are not the deposited construct lengths.

| Entry | Protein and selected chain | Resolution | Observed C-alpha residues | Main limitation |
|---|---|---:|---:|---|
| [1GJO](https://www.rcsb.org/structure/1GJO) | FGFR2 A | 2.4 A | 280 | WT and unphosphorylated, but not demonstrated to capture authentic autoinhibition; sulfate present. |
| [2PSQ](https://www.rcsb.org/structure/2PSQ) | FGFR2 A | 2.4 A | 287 | Partially active crystal conformation; sulfate and lattice confounding. |
| [2PVF](https://www.rcsb.org/structure/2PVF) | FGFR2 A | 1.8 A | 288 | C491A; phosphotyrosines 586, 656, 657; ATP analog, magnesium and substrate peptide. |
| [3RI1](https://www.rcsb.org/structure/3RI1) | FGFR2 A | 2.1 A | 285 | WT but ARQ069 inhibitor bound; not an unperturbed inactive control. |
| [3KY2](https://www.rcsb.org/structure/3KY2) | FGFR1 A | 2.7 A | 307 | Different isoform, C488A/C584S, and missing 488–489; supports inactive template modeling, not direct FGFR2 substitution. |

The exact FGFR2 reference used by this repository matches every observed 1GJO,
2PSQ and 3RI1 chain-A C-alpha identity at the author residue number. 2PVF matches
287/288 identities after PTR-to-Y canonicalization; its sole sequence mismatch is
C491A. There is no author-number offset to the assay's FGFR2 reference.

Observed internal gaps: 1GJO A 504–506 and 581–595; 2PVF A 509–511 and 589–596;
2PSQ A 583–593; 3RI1 A 583–595. Activation-loop residues are observed in these
chains, so their different conformations are not simply missing-loop imputation.

## Selected published preparations

The [original laboratory data page](https://yzhang.hpc.nyu.edu/MD_data_FGFR/) links
separate archives for WT and mutant simulations. Only the **first PDB member** of
each WT archive was streamed and retained. Reading stopped immediately after that
member, before any trajectory. The complete 909 MB and 520 MB archives were not
downloaded. Their first members were:

- `Autoinhibited_WT/FGFR2_wt_atp-3ky2.pdb`, saved as
  `structures/published-inactive-wt-first.pdb`.
- `Active_WT/FGFR2_wt_atp-2pvf.pdb`, saved as
  `structures/published-active-wt-first.pdb`.

The inactive preparation has 310 residues, author positions 459–768; the active
preparation has 299, positions 467–765. Both use blank author chain ID, preserve
C491A, canonicalize HIE/HID histidines to H, contain ATP/Mg, and contain no
phosphotyrosine. Their shared sequence is identical apart from their different
terminal coverage; both differ from the assay reference only at position 491.
The active preparation has rebuilt loops absent from its parent 2PVF coordinates.

Two coordinate-only sanity checks support distinct state geometries:

| Check | Inactive prepared model | Active prepared model |
|---|---:|---:|
| F645 chi1, degrees | +58.23 | -67.89 |
| Y656 CA to I707 CA, Angstrom | 13.88 | 29.29 |
| Y657 CA to I707 CA, Angstrom | 12.78 | 27.97 |

These checks validate the intended conformational contrast, not biological
accuracy or state occupancy. The 2023 paper emphasizes dynamics and state
populations; one prepared structure per state cannot represent that ensemble.

## Feature extraction requirements

1. Use `state-spec.json` for paths, member hashes, identity mapping and exclusions.
2. Restrict residue geometry to the common coordinate range and common atom
   identities. Strip ligands, water and hydrogens from protein-context features.
3. Exclude position 491 as a query and structural context position because neither
   preparation has the assay's WT cysteine there.
4. The predeclared common mask additionally excludes 509–511 and 581–596, which
   contain rebuilt or flexible loops. This union with 491 applies to both query
   and context residues in both states and was fixed before feature extraction
   or scoring.
5. Represent missing/unsupported queries explicitly. Never silently copy geometry
   from a neighboring residue or treat a missing contact as evidence of absence.
6. Retain the same regression model, split, training labels and feature controls in
   the matched baseline. Any new improvement is conditional on these structural
   assumptions and remains developmental evidence.

`structure-audit.json` records coordinate counts, sequence checks, structural
sanity checks and source links. The two archive-audit JSON files record public
Drive IDs, member names, byte counts and SHA-256 hashes. Full-archive hashes are
unavailable because the archives were intentionally not read in full.
