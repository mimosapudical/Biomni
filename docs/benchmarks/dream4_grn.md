# DREAM4 perturbational direct-network reconstruction

## Research question and scope

Can an agent reconstruct **direct, directed gene regulatory edges** from noisy,
controlled single-gene perturbations, rather than merely identify every downstream
expression response? A chain A -> B -> C can make A's knockout change C without
a direct A -> C edge. This benchmark measures unsigned direct topology, not
activation/repression signs, clinical causality, or factual biomedical recall.

This is a controlled **in-silico** benchmark built from one existing dataset,
not a newly collected biological dataset, a new inference algorithm, an official
Biomni-E2 benchmark, or a replication of the complete original DREAM scoring.

## Candidate selection (2026-10-17)

| Candidate | Strength | Reason for selection/exclusion |
| --- | --- | --- |
| DREAM4, 5 size-10 + 5 size-100 networks | Complete simulator direct-edge truth, labeled interventions, small pinned archive | Selected for a controlled direct-versus-indirect test; old and synthetic, so external validity is limited |
| CausalBench / Replogle Perturb-seq | Real K562/RPE1 interventional single-cell data; important contemporary comparison | Not selected: incomplete direct-edge truth, cell/guide-level splitting and licensing require additional decisions; held-out expression changes alone do not establish direct edges |
| BEELINE single-cell networks | Established systematic GRN methodology with synthetic and experimental data | Not selected: largely observational/pseudotime inputs, mixed reference definitions; a different information regime |
| Modern perturbation-prediction suites / Systema | Expose mean-baseline and metric confounds in real single-cell prediction | Not selected: expression prediction is not direct topology recovery; useful for interpreting limitations |

Repository-wide inspection found only HLE and Lab Bench under `biomni/task`;
Eval1 covers GWAS gene prioritization and related reasoning rather than this
matrix-to-directed-network task. Inspection covered all 210 pull requests,
362 issue/PR records, 15 latest commits, genomics/genetics/systems biology,
descriptions, environment manifests and the evaluation modules. An outstanding
TruthSeq contribution already validates predicted relationships against Replogle
knockdowns. This task deliberately does not reimplement that data/tool integration.
No existing DREAM4 task/integration was found in that search scope.

## Dataset, pin and licensing

- Original data: DREAM4 In Silico Network Challenge, simulated by GeneNetWeaver.
- Distribution: Paul Shannon's Bioconductor `DREAM4` package, version 1.33.0.
- Fixed mirror commit: `d901bd60f4c621ad2df0244931c60a9f3780fed8`.
- [Package DESCRIPTION and declared GPL license](https://github.com/bioconductor-source/DREAM4/blob/d901bd60f4c621ad2df0244931c60a9f3780fed8/DESCRIPTION).
- [Fixed archive download](https://raw.githubusercontent.com/bioconductor-source/DREAM4/d901bd60f4c621ad2df0244931c60a9f3780fed8/inst/extdata/lightlyProcessedDownloadedData.tar.gz).
- Archive SHA256: `1dede2b7578c552da2521c2ca31ebbcdaa04245eb7e54084e217d46731a871b2`.

The package is deprecated but its fixed Git commit remains addressable. The task
checks the bytes' digest before parsing, reads exact regular tar members without
extracting them, and rejects missing members, duplicate labels, incomplete gold
labels, bad gene order and nonfinite data. Missing or corrupted data is an
environment/input error and is never scored as an agent failure. Alternate digests
are explicitly labeled `custom:<sha256>` in every evaluation result.

The source distribution declares GPL without a version qualifier. Treat the data
under that upstream license; this integration does not relicense it as Apache-2.0.
Download and use it separately; do not bundle its archive or derived expression
exports in Biomni distributions under an Apache-only notice. The new Python code
uses Biomni's license and does not incorporate GPL source code. The tiny fixture
is independently authored, contains a 10-gene chain with isolated genes, and is
not a subset of DREAM4. The full data retain all five size-100 networks and all
9,900 non-self candidate edges per network: no easy-edge subsampling.

## Frozen protocol: `dream4-direct-ko-kd-v1`

| Item | Definition |
| --- | --- |
| Development/reference | All five size-10 networks; labels may be used for development only |
| Evaluation | All five size-100 networks, kept fixed; no tuning with their labels |
| Input | Ordered gene names; one wildtype vector; n-by-n knockout and knockdown expression matrices |
| Intervention meaning | Row i perturbs gene i; columns are measured genes. KO removes transcription, KD halves transcription rate; noisy steady state |
| Allowed data/tools | Exported evidence only, local statistical software and code execution |
| Excluded evidence | Gold edges, evaluator state/source archive, external networks/databases, internet lookup, time series, multifactorial and double-KO data |
| Output | `{"scores": [[...], ...]}`; finite numeric n-by-n matrix; rows sources, columns targets, zero diagonal; higher means stronger evidence |
| Candidate universe | All ordered non-self pairs; same universe for every method |
| Primary metric | Average precision (AP), then equal-weight macro average across networks |
| Secondary metrics | AUROC, positive prevalence, precision/recall at fixed k=n, fraction of top-k false-positive mass that is indirect-only reachability |
| Invalid answers | AP=0 and included in macro denominator; format-valid fraction is reported separately |
| Batch mismatch | Missing/extra batch entries raise an orchestration error; use explicit `None` for an unanswered task |

There is no separately learned training set: this is evidence-conditioned network
inference. The size-10 reference set is for method/prompt development. Evaluation
interventions are inputs for a different underlying network, not random rows of a
shared expression matrix. Development and evaluation use the same simulator and
may share structural motifs, so this is not evidence of distribution-independent
generalization. `network_ids` allows explicit smoke subsets and always appears in
results; only all five size-100 networks count as the full benchmark.

### Metric and uncertainty rationale

For decreasing score thresholds, AP = sum over thresholds of
(recall increment × precision). It is the noninterpolated, threshold-grouped
`sklearn.metrics.average_precision_score`, **not trapezoidal PR-AUC** and not the
original DREAM challenge's integrated score or p-value ranking. All-equal scores
give AP=π, the positive-edge prevalence. Finite random rankings have a slightly
different expected AP, so R2 supplies an empirical repeated-seed null.

AUROC = P(score_positive > score_negative) + 0.5 P(tie). It does not directly
measure the fraction of reported edges that are real. In particular:

precision = TPR × π / [TPR × π + FPR × (1 − π)].

At small π, a modest FPR can overwhelm true edges despite an attractive AUROC.
Therefore AP is primary and every valid per-network result includes prevalence.
See Saito and Rehmsmeier and BEELINE below. Scores are rankings, not claimed
probabilities: Brier/calibration scores would have no justified meaning here.

The top-k budget is k=n, chosen without gold edge counts or model outcomes. If a
score group straddles k, each edge receives its expected inclusion weight under
uniform tie breaking. This avoids choosing gene-index order as a hidden tiebreaker.
Indirect-only false positives have a directed path in gold but no direct edge.
This is a descriptive error decomposition, not a second positive class. When no
false positive exists, its fraction is `null` rather than an invented zero rate.

Edges sharing sources, targets and paths are dependent. Do not bootstrap individual
edges as independent observations. R3 reports a **descriptive** percentile bootstrap
of the five network APs (seed 0, 2,000 resamples). Five potentially related simulator
networks are too few for a reliable general-purpose 95% coverage claim. R2's seed
SD describes random-ranking variance on fixed data, not biological uncertainty.

Under a stable local linear model x = Bx + u, a total response can involve
(I − B)⁻¹ = I + B + B² + ... when the spectral radius of B is below one. Direct
edges correspond to B, whereas perturbations propagate through higher powers.
This explains the chain diagnostic; it is **not** a claim that the noisy nonlinear
DREAM simulator is globally linear, or that its direct graph is uniquely
identifiable from these finite steady-state observations. Inactive regulators,
weak effects, saturation, feedback and cancellation can hide direct edges.

### Leakage boundary

`get_example` exports only public expression evidence and returns no answer field.
Exporting deterministic inputs refuses to overwrite differing existing files.
The archive must be outside the public input directory. **This is not a sandbox**:
the A1 executor can otherwise read local files and use the network. For a scored
agent run, launch A1 in a separate process/container exposing only that public
directory; keep the archive, task instance, gold labels and result analysis outside
its filesystem and disable network access. Pass only `example['prompt']` to A1,
then move the response out for evaluation. Do not run an adversarial evaluation by
giving an agent the evaluator in the same REPL.

DREAM4 is public and old; pretrained contamination cannot be excluded. Renaming
genes does not prove decontamination. Future agent claims need evidence ablations
and contamination controls before being treated as causal reasoning evidence.

## Reproduction

Existing Biomni E1 provides NumPy and scikit-learn; no dependency is added. For
standalone task tests, install NumPy, scikit-learn, Pydantic 2 and pytest in the
development environment. The task never downloads anything during construction.

```python
from biomni.task.dream4_grn import dream4_grn

task = dream4_grn('/evaluator/dream4.tar.gz', '/public/inputs')
example = task.get_example(0)
# Agent runs separately with only /public mounted.
# responses is a list of five parsed JSON objects/JSON strings/output models or None.
metrics = task.evaluate(responses)
```

From the repository root, after downloading the fixed archive separately:

```bash
python -m pytest tests/test_dream4_grn.py -q
python -m tutorials.benchmarks.run_dream4_grn --round 1 --inputs /tmp/dream4_r1/inputs --output /tmp/benchmark_eval_r1_NEW.json
python -m tutorials.benchmarks.run_dream4_grn --round 2 --archive /evaluator/dream4.tar.gz --inputs /tmp/dream4_r2/inputs --output /tmp/benchmark_eval_r2_NEW.json
python -m tutorials.benchmarks.run_dream4_grn --round 3 --archive /evaluator/dream4.tar.gz --inputs /tmp/dream4_r3/inputs --output /tmp/benchmark_eval_r3_NEW.json
python -m tutorials.benchmarks.run_dream4_grn --round 4 --archive /evaluator/dream4.tar.gz --inputs /tmp/dream4_r4/inputs --output /tmp/benchmark_eval_r4_NEW.json
```

Every output filename must be fresh. R4 is an environment/readiness record, not a
fake or mocked agent run. It never makes a live call. Actual A1 runs must follow
the isolation boundary above and record exact model/version, prompt, repeated
trial seeds, wall time and raw responses. If runtime/access are unavailable,
report `not_run`, never model failure. MCZ in R3 uses the core median-corrected
standardization on the allowed one WT vector plus n KO profiles (sample SD).
It does not reproduce the original paper's additional WT observations from time
series, resampling, or ensemble pipelines, and its scores must not be compared
as an exact reproduction of original challenge results.

## Related research and methodological sources

1. Schaffter, Marbach and Floreano (2011), **GeneNetWeaver**, Bioinformatics.
   [DOI](https://doi.org/10.1093/bioinformatics/btr373). Simulator truth, systematic
   error profiling and the direct-versus-indirect distinction.
2. Pratapa et al. (2020), **BEELINE**, Nature Methods.
   [DOI](https://doi.org/10.1038/s41592-019-0690-6). Controlled GRN benchmarking,
   realistic limitations, precision-recall evaluation and comparable information.
3. Chevalley et al. (2025), **CausalBench Challenge**, CLeaR/PMLR 275.
   [Paper](https://proceedings.mlr.press/v275/chevalley25a.html). Modern real
   perturbational network inference; the selected synthetic task is complementary
   and cannot support the same real-cell claims.
4. **Systema** (2026 volume), Nature Biotechnology.
   [DOI](https://doi.org/10.1038/s41587-025-02777-8). Separating perturbation-specific
   effects from systematic variation; reinforces the need for trivial baselines.
5. Greenfield et al. (2010), **DREAM4: Combining Genetic and Dynamic Information**,
   PLOS ONE. [DOI](https://doi.org/10.1371/journal.pone.0013397). Source for MCZ core
   and its limitations; a baseline reference, not a new method claim.
6. Saito and Rehmsmeier (2015), PLOS ONE.
   [DOI](https://doi.org/10.1371/journal.pone.0118432). Authoritative methodology
   for precision-recall versus ROC under class imbalance.

## Authorship significance

This adds a reusable, reproducible evaluation asset and an error analysis focused
on a scientifically real distinction. It is more than an API wrapper, but a
single integration of an old synthetic dataset, without live agent validation or
real-cell external validity, does **not** by itself establish equivalence to
Biomni's “10+ significant & integrated tool contributions or equivalent” standard.
That decision belongs to maintainers and depends on actual integration and impact.

The smallest coherent follow-up is to establish an agent-level study on this same
frozen task: repeated isolated A1 trials, matched-information classical baselines,
and an evidence-removal/control analysis that separates reconstruction from
memorization. If useful to maintainers, the same causal-evaluation contribution
family can then be externally validated on one properly licensed real perturbation
dataset, with an explicitly weaker total-effect estimand rather than falsely
claiming complete direct-edge ground truth. Those are follow-ups, not unrelated
features padded into this PR.
