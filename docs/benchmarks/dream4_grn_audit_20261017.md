# DREAM4 contribution audit — 2026-10-17

## Status and exact provenance

Implementation and offline evaluation are complete. **No PR has been opened.**
The complete root-level pytest invocation remains environment-blocked, and a live
Biomni attempt has not been executed. This branch is a reviewable research pilot,
not an unqualified all-green or authorship-ready submission.

- Upstream `main`: `400c1f366b96a35ca253e13c9b06c5076af41d65`, fetched and checked again against the remote.
- Fork: `mimosapudical/Biomni`.
- Branch: `research/perturbation-grn-20261017`.
- Tested implementation commit: `ba11d1a7bb2e86afad227a63686b5f175472f588`.
- GitHub implementation commit: `ce1a5b3ba197849414f1707a39619d8738c9637f`. Its complete Git tree is identical to the local tested implementation (`4210e4b42d17ad26fe5a7703fdb1a87ad4ffd6e4`); only commit metadata differs because GitHub API credentials were used to publish it.
- Result records additionally contain SHA256 fingerprints of the exact task,
  audit script and synthetic fixture. Their implementation commit is the one above;
  later commits only add these audit/result records. Obtain the final result commit
  from the branch head, rather than inventing a self-referential SHA inside this file.
- Dataset commit: `d901bd60f4c621ad2df0244931c60a9f3780fed8`.
- Dataset archive SHA256: `1dede2b7578c552da2521c2ca31ebbcdaa04245eb7e54084e217d46731a871b2`.

## Initial checklist recheck

- [x] Research contribution: controlled direct-GRN reconstruction, with an indirect-path error diagnostic.
- [x] Read README, CONTRIBUTION, DETAILS, every task module, relevant tool/description code, Eval1, dependency/environment manifests and pre-commit config.
- [x] Upstream has no `.github` workflow directory or PR template at the recorded SHA; do not invent CI commands.
- [x] Overlap search: all 210 PRs, all 362 issue/PR records, recent 15 commits and repository code. Related outstanding Replogle validation work excluded from scope.
- [x] Dataset candidates, pin, license boundary, frozen evidence/label split, full candidate universe and leakage limitations recorded in [protocol](dream4_grn.md).
- [x] Reused `base_task.py` lines 1–18 and `lab_bench.py` lines 65–115 conventions; did not copy their multiple-choice evaluator or global RNG mutation.
- [x] All five required task methods, plus the existing iterator/prompt interface.
- [x] Unit fixtures separate from full scientific data; no network-dependent unit tests.
- [x] Focused tests, offline suite, compilation, installed-wheel smoke, real-data load, all 14 configured hook checks by the equivalent execution described below, and diff check pass.
- [x] Original root-pytest and pre-commit launcher failures reproduced on untouched upstream and disclosed, not silently skipped as passes.
- [x] Full implementation diff reread; existing task/tool methods compared again. No duplicate GRN integration, unrelated refactor or added dependency.
- [x] Fresh filenames for each experiment/audit; no numerical experimental result overwritten.
- [ ] Root-level `python -m pytest -q` is not green in this environment: missing optional MCP runtime.
- [ ] Live isolated A1 evaluation is not available: missing agent runtime, credentials and custom endpoint. Conditional R4 is recorded as `not_run`, not model failure.
- [ ] Do not open an upstream PR or claim full readiness while those environment/validation gaps remain unresolved.

## Files and functions

| File | Change |
| --- | --- |
| `biomni/task/dream4_grn.py` | `_table` 27–39; `_gold` 42–62; `_top_weights` 65–70; `dream4_grn.__init__` 85–154; `__len__` 156–157; `get_example` 159–162; `get_iterator` 164–166; `output_class` 168–175; `evaluate` 177–235 |
| `tests/test_dream4_grn.py` | 34 offline cases: load determinism, pins, missing/corrupt data, valid/malformed answers, hand metrics, directionality, ties, no overwrite, MCZ core correctness |
| `tests/fixtures/dream4_grn/chain.json` | Independently authored 10-gene chain and seven isolated genes; not original DREAM4 data |
| `tutorials/benchmarks/run_dream4_grn.py` | `median_corrected_zscore`, synthetic sanity packing, seeded nulls, descriptive network bootstrap, environment record, exclusive result writes |
| `docs/benchmarks/dream4_grn.md` | Scientific protocol, licensing/pin, literature, equations, uncertainty, leakage and authorship assessment |
| `README.md` | Task discoverability link |
| `docs/benchmarks/dream4_grn_results/` | Exact R1–R4 results and validation receipts |
| This audit | Review status, outcomes, commands and draft PR material |

No existing tool, description, evaluator, agent, dependency specification or
original task implementation is modified. Full data and gold edges are absent
from both the commit and built wheel.

## Experimental results and interpretation

| Round | Hypothesis and exact result | Interpretation / smallest next step |
| --- | --- | --- |
| R1 | Two positives at ranks 1 and 3: AP=0.8333333333333333; AUROC=0.9943181818181819 | Matches AP=5/6 and AUROC=175/176; G1->G3 is an indirect false edge in the fixture. Keep metric frozen, test nulls. |
| R2 | Constant macro AP=0.020686868686868688; every AUROC=0.5. Random seeds 0..19: mean macro AP=0.021822142323550106, seed SD=0.0006752452445599101 | Chance behavior is sensible; seed SD is conditional algorithmic variance, not a biological CI. Proceed to a label-free standard baseline. |
| R3 | MCZ core macro AP=0.3582888775018058; descriptive network-bootstrap interval=[0.25376360365601036, 0.4390549367141487] | Above empirical random reference, while indirect-path false positives remain. Do not tune metric or baseline after these outcomes. |
| R4 | `not_run`; runtime=false, provider credentials=false, custom endpoint=false | No observation of agent performance; no claim of agent failure/success. Provision isolated runtime before agent-level claims. |

All final experiments were rerun from the recorded implementation commit with
new `benchmark_eval_rN_20261017b.json` names; initial runs are retained separately.
Each final JSON includes the hypothesis, exact result, interpretation, biological
error mode, uncertainty, 3–4 literature/methodology sources, statistical principle
and smallest next step. Raw JSON snapshots precede repository formatting;
formatting is checked for identical parsed numerical content.

| 100-gene network | MCZ AP | MCZ AUROC | Precision@100 | Indirect fraction among MCZ false positives | Indirect fraction among all nonedges |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.4876513620629756 | 0.9080520081522755 | 0.63 | 1.0 | 0.04761415055532722 |
| 2 | 0.38818273840568995 | 0.7969725758281285 | 0.72 | 0.8571428571428571 | 0.04331157393016285 |
| 3 | 0.37953194676831414 | 0.8294260161957225 | 0.60 | 0.95 | 0.22462648119526085 |
| 4 | 0.3661602986909084 | 0.8422129654041643 | 0.61 | 0.9743589743589743 | 0.2172566828362067 |
| 5 | 0.16991804158114115 | 0.7652604738527987 | 0.32 | 0.9264705882352942 | 0.07829401462861886 |

The indirect error concentration is larger than its availability among nonedges.
This is a descriptive mechanistic check; dependent edges and only five simulator
networks preclude an independent-edge significance argument. The local stable
linear analogy is total response B+B²+... versus direct adjacency B. The actual
simulator is noisy and nonlinear: direct identifiability is not guaranteed.

Relevant sources are [GeneNetWeaver](https://doi.org/10.1093/bioinformatics/btr373),
[BEELINE, Nature Methods](https://doi.org/10.1038/s41592-019-0690-6),
[CausalBench, CLeaR 2025](https://proceedings.mlr.press/v275/chevalley25a.html),
[Systema, Nature Biotechnology](https://doi.org/10.1038/s41587-025-02777-8),
[MCZ](https://doi.org/10.1371/journal.pone.0013397), and
[precision-recall methodology](https://doi.org/10.1371/journal.pone.0118432).
These justify the task's interpretation and limitations, not a new theorem or
claim to outperform modern real-cell perturbation methods.

## Validation commands and exact outcomes

The [machine-readable validation receipt](dream4_grn_results/audit_r12_20261017_validation.json)
contains complete stdout and exit codes for the 16 checks below. The
[environment reproduction log](dream4_grn_results/audit_r12_20261017_environment.txt)
contains both failures on unchanged upstream and the feature-absence check.

| Command | Exact outcome |
| --- | --- |
| `python -m pytest tests/test_dream4_grn.py -q` | Initial focused suite: 32 passed in 1.98s; final suite adds two MCZ baseline cases |
| `python -m pytest tests -q` | 34 passed in 1.04s, exit 0 |
| `python -m compileall -q biomni tutorials/benchmarks tests` | exit 0 |
| `git diff --cached --check` | exit 0 |
| `python -m pre_commit run HOOK --all-files` | All 13 non-Biome hooks passed, exit 0 each |
| `npm exec --yes --package @biomejs/biome@2.3.10 -- biome format --files-ignore-unknown=true --no-errors-on-unmatched tests/fixtures/dream4_grn/chain.json docs/benchmarks/dream4_grn_results` | Exact configured Biome version; initial check: 5 files, no fixes, exit 0; final expanded receipts checked again |
| `python -m build --wheel --no-isolation --outdir /workspace/scratch/7860330af65a/dream4_build` | Built `biomni-0.0.8-py3-none-any.whl`, exit 0 |
| `python -m pip install --no-deps --target /workspace/scratch/7860330af65a/dream4_wheel_install WHEEL` | exit 0; task import from installed wheel and real development-data evaluation passed: AP=0.16666666666666666, valid fraction=1.0, AUROC=0.5 |
| Real pinned-data load, both splits | 5 development networks of size 10 and 5 evaluation networks of size 100 loaded and exported successfully |
| Feature tests copied outside checkout and run with `PYTHONPATH` set to pristine upstream | Expected collection failure: upstream lacks `biomni.task.dream4_grn`; the new feature is absent there, not a pre-existing metric regression |
| `python -m pytest -q`, both implementation and pristine upstream | Both exit 2 at collection: optional MCP tutorial imports missing `mcp`; neither is presented as a successful full-repository run |
| `python -m pre_commit run --all-files`, pristine upstream | Launcher exits 3 installing Biome's Git package: npm EALLOWGIT. Same failure initially on implementation branch. Registered Biome 2.3.10 binary checks passed instead; upstream config was not modified. |

The 13 hook IDs above are `ruff-check`, `ruff-format`, `detect-private-key`,
`check-ast`, `end-of-file-fixer`, `mixed-line-ending`, `trailing-whitespace`,
`check-case-conflict`, `check-added-large-files`, `check-toml`, `check-yaml`,
`check-merge-conflict`, and `no-commit-to-branch`.

## Remaining limits and authorship assessment

The asset is reproducible and scientifically interpretable, but this is an old
synthetic dataset with only five evaluation networks, public gold and potential
pretraining contamination, a restricted KO/KD evidence regime, and no live agent
study. The source data are separately GPL-declared; no original data are bundled.
Actual agent execution requires a sandbox that hides gold and blocks external
retrieval; a prompt instruction alone is insufficient.

This is a useful benchmark contribution, **not yet enough evidence for a major
authorship-equivalent claim**. The minimal coherent continuation is a repeated,
isolated Biomni evaluation plus matched-information/evidence-removal controls,
followed by real perturbation external validation if maintainers want that
contribution family. It should not be padded with unrelated tool wrappers.

## Draft PR title and body — not submitted

Title: **Add a pinned perturbational DREAM4 benchmark for direct GRN reconstruction**

Body:

> Biomni's task collection does not currently measure directed regulatory network
> reconstruction from controlled perturbation evidence. This adds a DREAM4 task
> that distinguishes direct edges from indirect downstream effects, with five
> size-10 development networks and five size-100 evaluation networks.
>
> The task verifies a fixed dataset commit/SHA256, exports evidence without gold,
> validates dense directed score matrices, and reports macro average precision,
> AUROC, fixed-budget precision/recall and indirect-path false positives. It adds
> independently authored offline fixtures, 34 tests, baseline audits and a frozen
> scientific protocol. No new dependency or existing tool/agent behavior changes;
> original data are downloaded separately under their upstream GPL declaration.
>
> On the fixed evaluation set, random rankings over 20 seeds achieve mean macro
> AP 0.02182, versus 0.35829 for an MCZ core baseline. This is a synthetic topology
> benchmark, not a new algorithm, real-cell causal validation, or an official E2
> designation. No live Biomni result is claimed.
>
> Validation: 34 offline tests, compilation, installed-wheel smoke and equivalent
> configured hook checks pass. Full root pytest remains blocked by a missing MCP
> dependency in an unchanged tutorial; the failure reproduces on pristine upstream.
> The normal pre-commit launcher also has a reproduced npm Git-install environment
> blocker; the exact Biome version was checked through its registry binary.
>
> This description remains a draft. Do not submit until the outstanding environment
> checks are resolved and the readiness statement can be updated with real results.
