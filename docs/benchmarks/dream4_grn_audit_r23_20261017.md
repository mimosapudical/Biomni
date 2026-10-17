# DREAM4 runtime follow-up audit — 2026-10-17, round 23

This supersedes the runtime status in the original audit; it does not overwrite
any previous experiment or alter the frozen task, baseline, or metric.

## Provenance and scope

- Upstream main remains `400c1f366b96a35ca253e13c9b06c5076af41d65` after fetching.
- Tested fork branch: `research/perturbation-grn-20261017`.
- Tested implementation/result commit: `528de05be101d33a7bb6460ce5de7be1f60c13b7`.
- No task, agent, dependency manifest, or existing test code changed this round.
- A refreshed overlap search found no DREAM4 benchmark contribution. The existing
  perturbation-validation tool proposal remains a different scope.

## Completed environment repairs

A separate validation virtual environment now contains MCP, LangChain,
LangGraph, OpenAI/Anthropic adapters, pytest-asyncio and SOCKS transport support.
`from biomni.agent.a1 import A1` succeeds; `python -m pip check` reports no broken
requirements. This does not claim the complete biomedical software environment
has been installed.

Using an npm 10 installation outside the repository resolves the earlier Biome
Git-package installation error. The unchanged repository pre-commit configuration
now runs normally: all 14 hooks pass. No replacement hook or skipped Biome check
is used for this result.

## Exact validation outcomes

| Command | Result |
| --- | --- |
| `python -m pytest tests/test_dream4_grn.py -q` | 34 passed |
| `python -m compileall -q biomni tutorials/benchmarks tests` | Exit 0 |
| `python -m pip check` | No broken requirements found |
| `python -m pre_commit run --all-files` with npm 10 on PATH | All 14 configured hooks passed |
| `git diff --check` | Exit 0 |
| `python -m pytest -q` in validation environment | 34 passed, 1 pre-existing async tutorial failure |
| Same default command on untouched upstream | Same async tutorial failure |
| `python -m pytest -q --asyncio-mode=auto` | 35 passed, with the qualification below |
| Same auto-mode command on untouched upstream | 1 passed, with the same qualification |

The MCP tutorial catches connection exceptions and returns `False` instead of
asserting or raising. Pytest's reported pass in auto mode therefore does **not**
establish an operational MCP server. An explicit verbose execution printed a
server import/connection failure followed by `PASSED`. We neither fix this
unrelated upstream tutorial in the benchmark PR nor present that result as a
successful live integration. The default async-mode failure is also preserved,
not hidden by the alternative invocation.

## Live-model gate

A real provider preflight through Biomni's `get_llm()` now reaches authentication
and raises:

```text
TypeError: Anthropic authentication failed: no API key or authorization
credentials were provided.
```

This is an environment/authentication failure, not a model or benchmark failure.
The refreshed R4 record reports `runtime_available=true`, credentials and custom
endpoint absent, and `status=not_run`. It is still a readiness record, not an
agent trajectory. No AP, biological error, or model variance can be inferred.

The machine-readable receipts are
[runtime validation](dream4_grn_results/audit_r23_20261017_runtime_validation.json)
and [R4 readiness](dream4_grn_results/benchmark_eval_r4_20261017c.json).

## Interpretation and smallest next action

Hypothesis: provisioning the Python/tooling runtime would remove the local
validation blockers. The import, package consistency and exact pre-commit results
support that limited hypothesis. They do not establish usable model access or
successful evidence-based GRN inference.

The scientific interpretation remains the frozen one in the protocol:
[BEELINE](https://doi.org/10.1038/s41592-019-0690-6),
[CausalBench](https://proceedings.mlr.press/v275/chevalley25a.html), and
[Systema](https://doi.org/10.1038/s41587-025-02777-8) motivate comparable evidence,
explicit baselines and careful perturbational evaluation. Infrastructure
missingness supplies no scientific performance observation. No metric or
statistical claim is changed after seeing these environment outcomes.

The next required input is a usable model provider configured securely in the
execution environment: for example `ANTHROPIC_API_KEY`, or an appropriately
configured alternative supported by Biomni. Do not paste secrets into a public
issue, commit or PR. Once available, run isolated evidence-only A1 trials, retain
raw responses and score them with the unchanged evaluator. Validate isolation
before reporting a scored trajectory; runtime import alone does not verify it.

## Checklist and submission decision

- [x] Latest upstream SHA and fork head checked.
- [x] Existing implementation and overlap rechecked; no duplicate implementation.
- [x] Agent import and tooling dependencies provisioned outside repository.
- [x] Focused tests, compilation, exact hooks and diff check pass.
- [x] Full pytest behavior checked against unchanged upstream and disclosed.
- [x] Fresh records preserve earlier outcomes and original numerical baselines.
- [ ] Usable model authentication and isolated A1 evaluation.
- [ ] All requested submission gates met; upstream PR opened.

No upstream PR is opened while the requested live-model gate is unresolved.
The PR title/body in the original audit remain drafts. No new authorship claim
is warranted by installing the runtime.
