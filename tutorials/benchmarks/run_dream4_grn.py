"""Reproduce DREAM4 baseline audits without downloading data or calling an LLM.

Run from the repository root with ``python -m tutorials.benchmarks.run_dream4_grn``.
Every result is opened exclusively: reruns need a fresh --output filename.
"""

import argparse
import hashlib
import importlib.util
import io
import json
import os
import platform
import subprocess
import tarfile
from importlib.metadata import version
from pathlib import Path

import numpy as np
from biomni.task.dream4_grn import dream4_grn

SOURCES = {
    "gnw": "https://doi.org/10.1093/bioinformatics/btr373",
    "beeline": "https://doi.org/10.1038/s41592-019-0690-6",
    "mcz": "https://doi.org/10.1371/journal.pone.0013397",
    "precision_recall": "https://doi.org/10.1371/journal.pone.0118432",
    "causalbench": "https://proceedings.mlr.press/v275/chevalley25a.html",
    "systema": "https://doi.org/10.1038/s41587-025-02777-8",
}


def median_corrected_zscore(evidence):
    """MCZ core score (Greenfield et al. 2010), without the ensemble/resampling pipeline.

    For each target, estimate its reference by the median across WT + single KOs,
    and its scale by the sample standard deviation over those same observations.
    """
    ko = np.asarray(evidence["knockouts"])
    observations = np.vstack([evidence["wildtype"], ko])
    center = np.median(observations, axis=0)
    scale = observations.std(axis=0, ddof=1)
    score = np.divide(np.abs(ko - center), scale, out=np.zeros_like(ko), where=scale > 0)
    np.fill_diagonal(score, 0)
    return score


def _synthetic_archive(path):
    """Package only the authored fixture; no original DREAM4 bytes are bundled."""
    fixture = json.loads((Path(__file__).resolve().parents[2] / "tests/fixtures/dream4_grn/chain.json").read_text())
    genes = fixture["genes"]
    matrices = {
        "wildtype": [fixture["wildtype"]],
        "knockouts": fixture["knockouts"],
        "knockdowns": fixture["knockdowns"],
    }
    files = {
        name + ".tsv": "\t".join(genes) + "\n" + "\n".join("\t".join(map(str, row)) for row in rows) + "\n"
        for name, rows in matrices.items()
    }
    positives = {tuple(edge) for edge in fixture["positive_edges"]}
    files["goldStandard.tsv"] = "".join(
        f"{a}\t{b}\t{int((a, b) in positives)}\n" for a in genes for b in genes if a != b
    )
    # Exclusive creation protects previous experimental inputs too.
    with path.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name, text in files.items():
            payload = text.encode()
            info = tarfile.TarInfo("lightlyProcessedDownloadedData/insilico_size10_1/" + name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize_runs(runs):
    values = [r["macro_average_precision"] for r in runs]
    return {
        "mean_macro_ap": float(np.mean(values)),
        "min_macro_ap": min(values),
        "max_macro_ap": max(values),
        "seed_sd_macro_ap": float(np.std(values, ddof=1)) if len(values) > 1 else None,
    }


def run(args):
    rng = np.random.default_rng(0)
    source_keys = {
        1: ["gnw", "beeline", "precision_recall"],
        2: ["precision_recall", "beeline", "systema"],
        3: ["mcz", "gnw", "causalbench", "systema"],
        4: ["causalbench", "beeline", "systema"],
    }
    root = Path(__file__).resolve().parents[2]
    source_paths = [
        "biomni/task/dream4_grn.py",
        "tutorials/benchmarks/run_dream4_grn.py",
        "tests/fixtures/dream4_grn/chain.json",
    ]
    result = {
        "round": args.round,
        "environment": {
            "python": platform.python_version(),
            "numpy": version("numpy"),
            "sklearn": version("scikit-learn"),
        },
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "working_tree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()),
        "source_sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in source_paths},
        "sources": {key: SOURCES[key] for key in source_keys[args.round]},
    }
    if args.round == 1:
        digest = _synthetic_archive(args.inputs.parent / (args.output.stem + "_toy.tar.gz"))
        task = dream4_grn(
            args.inputs.parent / (args.output.stem + "_toy.tar.gz"), args.inputs, "development", (1,), digest
        )
        s = np.zeros((10, 10))
        s[0, 1], s[0, 2], s[1, 2] = 3, 2, 1
        result["hypothesis"] = "Hand-ranked direct edges at ranks 1 and 3 give AP=5/6; an indirect chain edge is false."
        result["evaluation"] = task.evaluate([{"scores": s.tolist()}])
        actual = result["evaluation"]["macro_average_precision"]
        if not np.isclose(actual, 5 / 6):
            raise AssertionError(f"Metric sanity failed: {actual}")
        result.update(
            interpretation="Supports metric/orientation sanity, not biological validity; one authored noiseless chain.",
            biological_error_mode="G1 affects G3 through G2 without directly regulating G3.",
            uncertainty="Exact finite example, no population confidence interval.",
            principle="AP sums precision at positive recall increments; AUROC counts correctly ordered positive-negative pairs.",
            next_change="Keep metric frozen; run constant/random baselines on all evaluation networks.",
        )
        return result
    if args.archive is None:
        raise ValueError("--archive is required for rounds 2, 3, 4")
    task = dream4_grn(args.archive, args.inputs)
    if args.round == 4:
        runtime = all(importlib.util.find_spec(name) is not None for name in ["langchain_core", "langgraph"])
        credentials = any(
            bool(os.getenv(name)) for name in ["OPENAI_API_KEY", "ANTHROPIC_API_KEY", "CUSTOM_MODEL_API_KEY"]
        )
        endpoint = bool(os.getenv("CUSTOM_MODEL_BASE_URL"))
        result.update(
            hypothesis="Biomni can perform evidence-based direct GRN reconstruction through its normal agent loop.",
            status="not_run",
            reason="This baseline audit does not make live calls; run A1 in an isolated evidence-only environment.",
            runtime_available=runtime,
            provider_credentials_present=credentials,
            custom_endpoint_present=endpoint,
            dataset_revision=task.revision,
            dataset_sha256=task.dataset_sha256,
            interpretation="No agent-performance conclusion is possible. Missing runtime/access is an environment blocker.",
            biological_error_mode="Not observed: no live biological inference was executed.",
            uncertainty="Agent performance and variance are unmeasured.",
            principle="Missing measurements cannot be imputed as failures or successes.",
            next_change="Provision an isolated Biomni environment and run fixed-model repeated trials without changing scoring.",
        )
        return result
    evidence = [json.loads(Path(example["input_path"]).read_text()) for example in task.get_iterator()]
    if args.round == 2:
        constant = [{"scores": np.zeros((100, 100)).tolist()} for _ in task.get_iterator()]
        constant_result = task.evaluate(constant)
        random_runs = []
        for seed in range(20):
            seed_rng = np.random.default_rng(seed)
            answers = []
            for _ in task.get_iterator():
                s = seed_rng.random((100, 100))
                np.fill_diagonal(s, 0)
                answers.append({"scores": s.tolist()})
            random_runs.append({"seed": seed, **task.evaluate(answers)})
        result.update(
            hypothesis="Constant rank scores yield prevalence AP and 0.5 AUROC; random ranks recover no useful network.",
            constant=constant_result,
            random_summary=summarize_runs(random_runs),
            random_runs=random_runs,
            interpretation="Finite random-ranking AP is not exactly prevalence; empirical seeds provide the null reference.",
            biological_error_mode="Random edges ignore direction, transcriptional activity and downstream structure.",
            uncertainty="20 seed SD measures algorithmic randomness conditional on these five networks, not biological sampling.",
            principle="Precision = TPR*pi / (TPR*pi + FPR*(1-pi)); sparse positives can coexist with deceptively high AUROC.",
            next_change="Run the frozen, label-free MCZ core baseline; do not tune it on evaluation labels.",
        )
    else:
        answers = [{"scores": median_corrected_zscore(e).tolist()} for e in evidence]
        evaluation = task.evaluate(answers)
        # Networks, not dependent edges, are the resampling units. Five clusters => descriptive only.
        ap = np.array([n["average_precision"] for n in evaluation["networks"]])
        samples = rng.choice(ap, size=(2000, len(ap)), replace=True).mean(axis=1)
        result.update(
            hypothesis="The published MCZ core recovers direct edges above null but retains indirect-effect errors.",
            evaluation=evaluation,
            descriptive_network_bootstrap_95_percent=np.quantile(samples, [0.025, 0.975]).tolist(),
            interpretation="Compare with R2 and per-network diagnostics; MCZ is not a uniquely identifying causal estimator.",
            biological_error_mode="Knockout responses propagate through intermediates; inactive/weak regulators may be missed.",
            uncertainty="Descriptive percentile bootstrap over five networks, seed=0, B=2000; no guaranteed coverage or human generalization.",
            principle="Under a local linear approximation, total response includes B+B^2+...; direct adjacency is B alone.",
            next_change="Freeze code/metric and run the agent when runtime access is available.",
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round", type=int, choices=[1, 2, 3, 4], required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Result exists: use a fresh result/log filename")
    args.inputs.parent.mkdir(parents=True, exist_ok=True)
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Saved R{args.round}: {args.output}")


if __name__ == "__main__":
    main()
