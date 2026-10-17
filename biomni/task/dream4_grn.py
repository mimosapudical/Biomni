"""Pinned DREAM4 perturbational direct-network reconstruction (offline evaluation).

Full data are downloaded separately, never by the task. See
docs/benchmarks/dream4_grn.md for the frozen protocol and data license.
"""

import csv
import hashlib
import io
import json
import tarfile
from pathlib import Path

import numpy as np

from biomni.task.base_task import base_task

DATA_REVISION = "d901bd60f4c621ad2df0244931c60a9f3780fed8"
DATA_SHA256 = "1dede2b7578c552da2521c2ca31ebbcdaa04245eb7e54084e217d46731a871b2"
DATA_URL = (
    "https://raw.githubusercontent.com/bioconductor-source/DREAM4/"
    f"{DATA_REVISION}/inst/extdata/lightlyProcessedDownloadedData.tar.gz"
)
PROTOCOL = "dream4-direct-ko-kd-v1"


def _table(archive, member, n, rows):
    """Read only the named regular member; never extract an archive to disk."""
    info = archive.getmember(member)
    if not info.isfile():
        raise ValueError(f"Expected a regular data file: {member}")
    with archive.extractfile(info) as stream:
        records = list(csv.reader(io.StringIO(stream.read().decode("utf-8")), delimiter="\t"))
    if not records or records[0] != [f"G{i + 1}" for i in range(n)]:
        raise ValueError(f"Unexpected gene order: {member}")
    values = np.asarray(records[1:], dtype=float)
    if values.shape != (rows, n) or not np.isfinite(values).all():
        raise ValueError(f"Invalid expression matrix: {member}")
    return values


def _gold(archive, member, n):
    info = archive.getmember(member)
    if not info.isfile():
        raise ValueError(f"Expected a regular gold file: {member}")
    with archive.extractfile(info) as stream:
        records = list(csv.reader(io.StringIO(stream.read().decode("utf-8")), delimiter="\t"))
    genes = {f"G{i + 1}": i for i in range(n)}
    edges = {}
    for source, target, label in records:
        if source not in genes or target not in genes or source == target or label not in {"0", "1"}:
            raise ValueError(f"Invalid gold edge: {member}")
        pair = (genes[source], genes[target])
        if pair in edges:
            raise ValueError(f"Duplicate gold edge: {member}")
        edges[pair] = int(label)
    if len(edges) != n * (n - 1) or not 0 < sum(edges.values()) < len(edges):
        raise ValueError(f"Gold must label every non-self edge with both classes: {member}")
    labels = np.zeros((n, n), dtype=int)
    for pair, label in edges.items():
        labels[pair] = label
    return labels


def _top_weights(scores, k):
    """Expected top-k inclusion under uniform tie breaking at the boundary."""
    boundary = np.partition(scores, len(scores) - k)[len(scores) - k]
    above = scores > boundary
    tied = scores == boundary
    return above.astype(float) + tied * ((k - above.sum()) / tied.sum())


class dream4_grn(base_task):
    """One directed network per example; all candidate edges must receive scores.

    Args:
        path: Local pinned data archive, not an extracted directory.
        input_dir: Public evidence export directory, mount this alone into the agent.
        split: ``development`` = five size-10 networks; ``evaluation`` = five size-100 networks.
        network_ids: Optional explicit subset of 1..5 for smoke tests; reported in results.
        dataset_sha256: Explicit alternate digest for synthetic fixtures/custom replications.
            Alternate digests are always marked ``custom:...``, never as official DREAM4.
    """

    def __init__(self, path, input_dir, split="evaluation", network_ids=None, dataset_sha256=DATA_SHA256):
        if split not in {"development", "evaluation"}:
            raise ValueError("split must be development or evaluation")
        ids = tuple(range(1, 6)) if network_ids is None else tuple(network_ids)
        if not ids or len(set(ids)) != len(ids) or any(type(i) is not int or i not in range(1, 6) for i in ids):
            raise ValueError("network_ids must be unique integers in 1..5")
        ids = tuple(sorted(ids))
        path, input_dir = Path(path).resolve(), Path(input_dir).resolve()
        if path.is_relative_to(input_dir):
            raise ValueError("Evaluator archive must be outside the public input directory")
        # Hash and parse the same bytes: no unverified download and no file-change race.
        payload = path.read_bytes()
        actual_hash = hashlib.sha256(payload).hexdigest()
        if actual_hash != dataset_sha256:
            raise ValueError(f"Dataset SHA256 mismatch: expected {dataset_sha256}, got {actual_hash}")
        self.revision = DATA_REVISION if actual_hash == DATA_SHA256 else f"custom:{actual_hash}"
        self.dataset_sha256 = actual_hash
        self.split, self.network_ids = split, ids
        self.n = 10 if split == "development" else 100
        self._labels, self._examples = [], []
        evidence_items = []
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
            names = archive.getnames()
            if len(names) != len(set(names)):
                raise ValueError("Duplicate archive members")
            for i in ids:
                name = f"insilico_size{self.n}_{i}"
                prefix = f"lightlyProcessedDownloadedData/{name}/"
                try:
                    wt = _table(archive, prefix + "wildtype.tsv", self.n, 1)
                    ko = _table(archive, prefix + "knockouts.tsv", self.n, self.n)
                    kd = _table(archive, prefix + "knockdowns.tsv", self.n, self.n)
                    labels = _gold(archive, prefix + "goldStandard.tsv", self.n)
                except KeyError as exc:
                    raise ValueError(f"Missing required data member for {name}: {exc}") from exc
                evidence_items.append(
                    {
                        "network": name,
                        "genes": [f"G{j + 1}" for j in range(self.n)],
                        "wildtype": wt[0].tolist(),
                        "knockouts": ko.tolist(),
                        "knockdowns": kd.tolist(),
                    }
                )
                self._labels.append(labels)
        input_dir.mkdir(parents=True, exist_ok=True)
        for evidence in evidence_items:
            target = input_dir / f"{PROTOCOL}_{self.dataset_sha256[:12]}_{evidence['network']}.json"
            content = json.dumps(evidence, sort_keys=True, allow_nan=False) + "\n"
            # Repeated deterministic loads are safe; differing existing evidence is never overwritten.
            try:
                with target.open("x", encoding="utf-8") as stream:
                    stream.write(content)
            except FileExistsError:
                if target.read_text(encoding="utf-8") != content:
                    raise ValueError(f"Existing evidence differs; choose a fresh input_dir: {target}") from None
            prompt = (
                f"Infer the direct gene regulatory network from {target}. This is simulated DREAM4 evidence. "
                "Use only this evidence and local analysis software; no internet, external networks, gold labels, "
                "or evaluator files. Gene order is given by 'genes'. Rows of knockouts/knockdowns are single-gene "
                "perturbations in this same order; columns are measured genes. Wildtype is the unperturbed state. "
                "A knockout removes a gene's transcription; a knockdown halves its transcription rate. "
                "Expression is noisy. Downstream responses may be indirect: predict direct edges, not all effects. "
                f"Return JSON {{'scores': <{self.n} by {self.n} numeric matrix>}} with double-quoted keys. "
                "scores[i][j] ranks evidence for gene i directly regulating gene j. Higher means stronger; "
                "scores need not be probabilities. Supply every entry, including a zero diagonal. Self edges "
                "are excluded. The primary score is average precision over all directed non-self edges. "
                f"Secondary precision/recall use a fixed top-{self.n} budget with fractional boundary ties."
            )
            self._examples.append({"prompt": prompt, "input_path": str(target), "network": evidence["network"]})

    def __len__(self):
        return len(self._examples)

    def get_example(self, index=0):
        if type(index) is not int or not 0 <= index < len(self):
            raise IndexError("Example index out of range")
        return dict(self._examples[index])

    def get_iterator(self):
        for i in range(len(self)):
            yield self.get_example(i)

    def output_class(self):
        from pydantic import BaseModel, ConfigDict, Field

        class NetworkOutput(BaseModel):
            model_config = ConfigDict(extra="forbid", strict=True, allow_inf_nan=False)
            scores: list[list[float]] = Field(description="Dense source-by-target rank scores; zero diagonal")

        return NetworkOutput

    def evaluate(self, response):
        """Score a list of dicts, output models, JSON strings, or None (one per network).

        Invalid answers receive AP=0 and remain in the macro denominator.
        Missing/extra network answers raise an orchestration error, rather than dropping tasks.
        Dataset/environment errors are raised during loading, never counted as agent failures.
        """
        from sklearn.metrics import average_precision_score, roc_auc_score

        if not isinstance(response, (list, tuple)) or len(response) != len(self):
            raise ValueError("Provide exactly one response per example, in iterator order")
        mask = ~np.eye(self.n, dtype=bool)
        results = []
        for example, labels, answer in zip(self._examples, self._labels, response, strict=True):
            item = {"network": example["network"], "valid": False, "average_precision": 0.0}
            try:
                if isinstance(answer, str):
                    answer = json.loads(answer)
                parsed = self.output_class().model_validate(answer, from_attributes=True)
                scores = np.asarray(parsed.scores, dtype=float)
                if scores.shape != (self.n, self.n) or not np.isfinite(scores).all():
                    raise ValueError("scores must be a finite n-by-n matrix")
                if np.any(np.diag(scores) != 0):
                    raise ValueError("Self scores must be zero")
            except (ValueError, TypeError, OverflowError) as exc:
                item["error"] = str(exc)
                results.append(item)
                continue
            y, s = labels[mask], scores[mask]
            weights = _top_weights(s, self.n)
            true_positives = float(weights @ y)
            # Transitive reachability is a diagnostic only; never use it as direct-edge truth.
            reach = labels.astype(bool).copy()
            for j in range(self.n):
                reach |= reach[:, j, None] & reach[None, j, :]
            indirect = (reach & ~labels.astype(bool))[mask]
            false_positives = float(weights @ (1 - y))
            item.update(
                valid=True,
                average_precision=float(average_precision_score(y, s)),
                auroc=float(roc_auc_score(y, s)),
                positive_edges=int(y.sum()),
                candidate_edges=int(len(y)),
                prevalence=float(y.mean()),
                precision_at_n=true_positives / self.n,
                recall_at_n=true_positives / int(y.sum()),
                indirect_fp_fraction=float(weights @ indirect) / false_positives if false_positives else None,
            )
            results.append(item)
        return {
            "protocol": PROTOCOL,
            "dataset_revision": self.revision,
            "dataset_sha256": self.dataset_sha256,
            "split": self.split,
            "network_ids": list(self.network_ids),
            "macro_average_precision": float(np.mean([r["average_precision"] for r in results])),
            "valid_fraction": sum(r["valid"] for r in results) / len(self),
            "networks": results,
        }
