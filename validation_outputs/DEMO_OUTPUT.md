# GhostScope audit

```text
GhostScope 1.0.0
================

Positive topology control
GhostScope 1.0.0
================
Evidence    Tier A — Conditional multi-taxon topological support
Loci        2535523
Unrooted-6  signs=70/70  excluded=70/70  p_union=5e-22
Compatible  Ghost pulse within the certified family; Unexamined histories (multiple pulses, mixtures, ancestral structure, gene-tree bias)
Next        Check backbone, candidate coverage and gene-tree error with independent inference and held-out null simulations.
Scope       Evidence tier is conditional on its stated candidate family and observational assumptions; no tier proves the unique biological history.
Fixture     Illustrative high-power demonstration, not a replay or claim of detection in the 500-locus paper controls.

Ambiguous depth control
GhostScope 1.0.0
================
Evidence    Tier E — Depth signal unsupported after path calibration
Loci        100
Triplet     a=0.5  b=0.2  wedge=True
Depth       raw=0.009  B-clock=-0.002
Compatible  Ghost pulse; Direct ingroup pulse; Rate-heterogeneous pure tree; Unexamined histories
Flag        SCOPE_ROUTING: OPEN_TOPOLOGY_WEDGE
Flag        SCOPE_ROUTING: TWO_RATE_PURE_TREE_MOMENT_COMPATIBLE
Flag        SIGN_FLIP
Next        Add a recipient sister and a pure outgroup for rooted four-leaf contrasts; sample outer ladder clades for six-leaf contrasts.
Scope       Evidence tier is conditional on its stated candidate family and observational assumptions; no tier proves the unique biological history.

Paired-network trap
GhostScope 1.0.0
================
Source      SAVED_AUDIT_SUMMARY_NOT_RAW_DRAWS
Evidence    Tier D — Multiple histories / stage-specific diagnostics
Tail splits 0
Paired      ['synthetic', 100, 'default', None, 'single'] SINGLE_GRAPH_BF_ILLUSION_DETECTED
lnL delta   ghost - inflow: -2.0
Flag        SINGLE_GRAPH_BF_ILLUSION_DETECTED
Compatible  Multiple candidate histories remain unless independently compared.
Next        Inspect both windows, rare-tail visits and all TI nodes before model selection.
Scope       Recomputed counts/decisions from saved summaries. Original convergence diagnostics and means are reported, not independently recomputed from unavailable 500k trajectories. epsilon / estimated posterior mass(phi<epsilon), assuming Uniform(0,1). A finite-interval mass ratio, not an exact point-null Savage-Dickey BF or a between-network BF. Zero hits do not establish zero probability. Chain means of lnL are descriptive, not marginal likelihoods.
```

Full statistics and assumptions: `report.json`.
