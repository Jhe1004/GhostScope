# GhostScope

Audit topology, distance and competing-network evidence for ghost introgression.
Python ≥3.9; runtime dependencies: NumPy and SciPy. CLI and Python API. MIT.

**Scope and finite-genome power notice:** `GhostScope` is an identifiability-routing
and competing-network audit toolkit, **not** an uncalibrated black-box detector of
unique biological history. Asymptotic algebraic identifiability ($L \to \infty$)
does not guarantee full-competitor rejection at finite locus counts: in the
project's 500-locus controls, rooted four-leaf ghost replicates consistently
show $\hat{\Delta} > 0$ ($2/5$ latent-genealogy and $1/3$ 500-bp NJ replicates
passing the heuristic 0.015 GTEE tolerance envelope), whereas unrooted six-leaf
deep-ghost replicates reject $0\text{--}8 / 70$ direct ingroup cells under the
strict intersection-union criterion and are retained as `Tier D` (compatible
histories). See `validation_outputs/VALIDATION_CERTIFICATE.json` and
`DEVELOPMENT_REVIEW.md` for full validation boundaries.

## Quickstart

```sh
python -m pip install .
ghostscope demo
ghostscope trees -g genes.nwk --taxa A,S,B,C --outgroup O
```

Without installation:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m ghostscope.cli demo
```

Commands print a compact terminal report. Add `-o audit/` to save complete
`report.json`, `report.tsv` and `report.md`. Nothing is written without `-o`,
except the explicitly requested `make-bpp` suite and `reproduce-paper` outputs.
Python's own bytecode caching is controlled by `PYTHONDONTWRITEBYTECODE=1`.

The positive demo uses **2,535,523 synthetic topology observations**, represented
as compact counts and sampled with a fixed seed from the specified exact six-leaf
model. It demonstrates a high-power positive route; it does not claim that the
paper's 500-locus controls pass. The other examples show an ambiguous depth signal
and two competing graphs with large within-graph interval ratios.

## Inputs and the four audit steps

1. **Topology and candidate scope.** Start from rooted gene trees with a pure
   outgroup, or unrooted trees for a six-leaf ladder. Check the backbone before
   interpreting four-leaf or six-leaf contrasts.
2. **Distance and observational checks.** Inspect topology-conditioned full-path
   distances, six-moment compatibility, sequence saturation and optional
   multi-individual contrasts. A compatible rate-variable pure tree is a routing
   result; it does not reject a ghost model.
3. **Competing fits.** Generate and run comparable BPP graphs. A high conditional
   phi estimate or a large within-graph interval ratio does not select that graph.
4. **Numerical audit.** Check both chain windows, rare-tail visits and every TI
   node. Retain unresolved alternatives when the required gate fails.

### Gene trees and scan

```sh
ghostscope trees -g genes.nwk.gz -s '((((A,S),B),C),O);' \
  --recipient A --outgroup O --imap imap.txt -o audit/
ghostscope scan -g genes.nwk -s species.nwk --outgroup O --imap imap.txt
```

`--taxa A,S,B,C` gives ordered ladder roles directly. Six-leaf order is
`recipient,sister,outer1,outer2,outer3,outer4`. `--rooted` declares already rooted
input when no outgroup is supplied. Arbitrary Newick root placement is not treated
as evidence of biological rooting. Multiple monophyletic outgroups are accepted
as `--outgroup O1,O2`; nonmonophyletic and missing-outgroup loci are counted.

Newick supports branch lengths in scientific notation, support labels, quoted
names, comments, gzipped streams, tree-file directories and the project's
`genealogies.jsonl`. Named one-line TREE statements in author PhyloNet/NEXUS files
are also supported; numeric NEXUS TRANSLATE blocks are rejected explicitly.

`--imap` has two columns, `sample species`; BPP `species^individual` suffix
identifiers are supported. Alternatively, `--taxon-delim _` maps `A_1` to `A` and
`--taxon-delim '^'` maps `A^individual` to `A`. Exact labels already present in the
species tree are preserved. Loci with multiple copies of a mapped species are
reported and excluded from single-copy topology analysis. No representative is
silently chosen. Full paths survive induced-tree pruning; zero internal branches
are unresolved, regardless of bracket order.

Scan examines every ingroup tip. Recipients lacking the certified terminal
ladder get an explicit unresolved route, not a fabricated topology result.
Rankings are exploratory; scanning is not adjusted for selection across species.
For n>6, every supplied outer representative must pass its own full 70-candidate
reduction. Other backbones and multi-pulse networks are outside that certificate.

### Alignments and optional NJ

```sh
ghostscope seqs -i loci/ -s '((((A,S),B),C),O);' \
  --recipient A --outgroup O --taxon-delim _ --nj-trees
ghostscope seqs -i data.bpp --format bpp \
  --haplotypes-a 'A^one,A^two' --taxon-b 'B^one' --taxon-c 'C^one'
```

FASTA, sequential/interleaved PHYLIP and multi-block BPP are accepted. JC69 uses
only pairwise called A/C/G/T sites. Saturated or uncalled pairs are never assigned
an arbitrary finite distance. The output separates unconditional JC69 existence
from the finite-locus conditional NJ stream. Nonpositive NJ internal branches
are collapsed. Bootstrap blocks are loci; physical linkage needs coordinates and
is not automatically corrected.

Sequence multi-copy statistics include the raw squared mismatch contrast and its
conditional-iid-site debiasing. They are **measurement proxies**, not certified
latent genealogical distances. Tier B therefore cannot be earned from sequence
proxies alone. For known latent tree distances, the explicit route is:

```sh
ghostscope trees -g latent.nwk --haplotypes-a A1,A2 --taxon-b B --taxon-c C \
  --latent-genealogies --assume-exchangeable
```

The assumptions must be justified, not inferred from sample names. Survival
curvature outputs are exploratory: the six-point theorem requires a common
breakpoint-free demographic window and mode constraints that files alone do not
establish. Three-role site-pattern ratios are provided; the software does not
pretend to reproduce the official HyDe variance test.

### MCMC and competing graphs

```sh
ghostscope mcmc --chains run1/chain.mcmc.txt run2/chain.mcmc.txt
ghostscope mcmc --ghost g1.txt g2.txt --inflow i1.txt i2.txt --outflow o1.txt o2.txt
ghostscope mcmc --summary-json initial.json continuation.json
ghostscope mcmc --ti-dirs ti_run1/ ti_run2/
```

Chains must have identical headers and lengths; columns follow BPP names such as
`phi:5<-4:H<-R`. At least two independent chains are needed. Full and recent-half
windows use the project `R8_TWO_TIER_20261108_v1` formulas and thresholds. Fewer
than 40 draws per window are refused. Strict and operational gates stay separate.
The policy is less stringent than a claim of complete posterior convergence.

The primary reported ratio is epsilon / estimated mass(phi<epsilon), assuming a
Uniform(0,1) prior. This is a **finite-interval mass ratio**, not an exact
point-null Savage–Dickey Bayes factor or a between-network Bayes factor. Zero
hits are Monte Carlo nonresolution. Clopper–Pearson bounds are labelled
**iid-only**; correlated-chain coverage is not guaranteed. Time-block hit counts
and indicator ESS are reported. Posterior mean lnL differences are descriptive.

Saved summaries are explicitly marked as such. Later records replace earlier
records with the same ID. Initial cumulative eligibility is not relabelled as a
dual-window result. The R8 audit TSV schema is also accepted; integer hit counts
are reconstructed from its interval ratios only when they can be recovered
unambiguously, and the provenance is marked. Unsupported tables are rejected.
TI audits check all supplied nodes and both windows;
weights and endpoint coverage must also be established before integration.

See the [BPP manual](https://bpp.github.io/bpp-manual/bpp-4-manual/) and
[Stan's diagnostics reference](https://mc-stan.org/rstan/reference/Rhat.html).

### Generate BPP controls

```sh
ghostscope make-bpp -s '((A,B),C);' --recipient A --sister B --out-clade C \
  --seqfile data.bpp --imap imap.txt --nloci 100 -o suite/
sh suite/run_and_audit.sh
```

Generation does not start BPP. The suite has Ghost, Inflow and Outflow × default,
locusrate and clock2 × two independent seeds: 18 controls. Roles match the project
plant baseline: Ghost targets A; Inflow is C→B; Outflow is B→C. Full larger
backbones are retained when the anchors induce ((A,B),C). Both Ghost `tau-parent`
annotations are `yes`; the direct donor annotation is `no`. Two `yes` annotations
do **not** encode two independent bidirectional migration probabilities.

BPP ≥4.8 requires `jobname = chain`, yielding `chain.mcmc.txt`. Obsolete
`outfile/mcmcfile` options are not emitted. Existing nonempty output suites are
never overwritten. Sample counts and `nloci` are checked against the real input.
Input paths must not contain whitespace because BPP's file-name parser is
version-dependent. The shipped priors match the historical plant control and
must be reviewed for other studies. Ordinary untagged input labels are normalized
into a copied BPP alignment and suffix map inside the suite; the original input
is preserved and `NORMALIZED_INPUT_MAP.tsv` records every rename. Native BPP
labels/maps are reused without an unnecessary sequence copy. The backbone must
be fully bifurcating, and the rate-model suite needs at least two loci. The generated script runs sequentially; adapt
scheduling and threads for the available machine.

## Evidence tiers

| Tier | What it reports | What it does not establish |
| --- | --- | --- |
| A | Conditional rooted-four or unrooted-six topological support in the stated candidate family | A unique history among all networks; a calibrated GTEE false-positive rate |
| B | Negative latent multi-haplotype contrast under the canonical assumptions | Validity of squared sequence-distance proxies or exchangeability |
| C | Quasi-clock consistency screen and a positive calibrated depth contrast | Proof of a molecular clock or absence of ancestral structure |
| D | Multiple histories remain compatible, or the required topology lacks power | Rejection of a ghost history |
| E | A particular signal is unsupported, or numerical eligibility fails | Biological absence of ghost ancestry |

Open topology wedges, positive six-moment pure-tree fits and the multi-haplotype
confounded cone are **SCOPE_ROUTING**, never automatic rejection of ghosts.
Type-III warnings indicate a possible stem/crown ambiguity requiring additional
coupling and chronology conditions; clade size alone is not a KL=0 proof.

Rooted-four output includes Delta, I1, I5, bootstrap SE/CI and a zero-null Wald
test. The default 0.015 GTEE floor is a declared heuristic. Both point estimates
must lie within the floor, both CIs must overlap its interval, Delta must be
positive with CI>0 and one-sided normal p<0.01, and the dominance ratio must be≥3.
This resolves the directive's inconsistent point/CI rule conservatively; it has
not been calibrated for arbitrary gene-tree bias.

Six-leaf output retains **all 70 linear invariants**, including the eight-term
invariant for (0,a,e). A sign count alone cannot promote evidence: every
candidate-specific null must be excluded in the certified direction. The reported
maximum directional p is an intersection-union test of the candidate union under
the individual asymptotic normal approximations. Sparse cells, unresolved-locus
conditioning and systematic inferred-tree error remain limitations.

## Python API and validation

```python
from ghostscope import audit_trees, audit_chains
from ghostscope.io_parsers import read_trees
result = audit_trees(list(read_trees('genes.nwk')),
                     taxa=['A', 'S', 'B', 'C'], outgroups=['O'])
```

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -v
# From the root of PhyloIdentifiability_Public_Archive_SciAdv_v1:
PYTHONPATH=09_GhostScope_Toolkit_v1/src python3 -m ghostscope.cli reproduce-paper --project-root .
```

When installed standalone via `pip install ghostscope-1.0.0-py3-none-any.whl` outside the archive, `reproduce-paper` requires `--project-root` pointing to the root of `PhyloIdentifiability_Public_Archive_SciAdv_v1` (where `06_Empirical_Data_Audits_R5/` and `08_Higher_Order_4to6_and_ArbitraryN_R8/` reside). All replay outputs are written to `validation_outputs/` (`292 / 292` checks: `289` deterministic replay assertions across simulation controls, five empirical systems, and MCMC/TI traces + `3` live SHA-256 integrity checks verifying all frozen archive files against `PUBLIC_ARCHIVE_SHA256.tsv` via `tools/protect.py`). It returns a nonzero exit status if any replay assertion or frozen-file hash check fails. See `VALIDATION_REPORT_20261111.md` and `validation_outputs/BPP_COMPATIBILITY.json` for the control outcomes, empirical comparisons, and the historical `36 / 36` `BPP v4.8.7` control-file startup verification record.
