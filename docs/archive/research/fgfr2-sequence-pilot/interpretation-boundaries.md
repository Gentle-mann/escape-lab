> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# What the frozen sequence pilot can establish

This pilot tests whether ESM2-35M contextual sequence features, paired with a
small supervised head, improve retrieval of higher-activation FGFR2 variants
across held-out sequence blocks in the existing development assay. It does not
test whether we can physically simulate kinase activation.

A positive result would support another developmental step. It would not
establish a molecular mechanism, active-state occupancy, catalytic rate, drug
resistance, or external generalization. A negative result would reject this
particular representation, compression, head and training-data combination. It
would not rule out all sequence models or mechanistic approaches.

The model must distinguish disruption magnitude from functional direction. A
mutation can disrupt a protein while decreasing activity; destabilizing an
inhibitory interaction can instead increase activity. PreMode explicitly
separates distance from WT from biochemical direction and learns phenotype
direction with supervision. This motivates signed mutant-minus-WT features,
but does not guarantee that the fixed low-dimensional head will recover the
relevant direction from 26 high-activation examples.

Source: [PreMode, Nature Communications 2025](https://www.nature.com/articles/s41467-025-62318-4).

If this fixed pilot is negative, an independently specified next hypothesis is
to contrast consistently mapped active and inactive structural environments.
Candidate features could compare contact networks, packing, solvent exposure,
and WT/mutant compatibility in each state. Frozen inverse-folding scores could
provide a learned compatibility feature. These would remain proxies; they are
not binding energies, free-energy differences, or simulated growth rates.

Paired-state kinase free-energy research provides the biological rationale:
activation can reflect a shifted balance between active and inactive states.
However, the study also explains that adopting an active conformation is not
sufficient to establish increased catalytic activity. Construct sequence,
phosphorylation, ligands and missing residues would all need auditing before
interpreting a structural contrast.

Sources: [Gizzio et al., Nature Communications 2024](https://www.nature.com/articles/s41467-024-50812-0),
[Reveguk and Simonson, Protein Science 2024](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4918).

DFG orientation and active/inactive classification are distinct. One crystal or
predicted structure is not a distribution over conformations. No state-based
features, new models, or external outcomes are added to the present pilot.
