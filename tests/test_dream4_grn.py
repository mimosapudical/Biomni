"""Offline regression tests; fixture is authored synthetic data, not a DREAM4 subset."""

import hashlib
import io
import json
import socket
import tarfile
from pathlib import Path

import numpy as np
import pytest
from biomni.task.dream4_grn import DATA_REVISION, DATA_SHA256, DATA_URL, _top_weights, dream4_grn
from sklearn.metrics import average_precision_score
from tutorials.benchmarks.run_dream4_grn import median_corrected_zscore

FIXTURE = Path(__file__).parent / "fixtures/dream4_grn/chain.json"


def build_archive(path, missing=None, transform=None):
    fixture = json.loads(FIXTURE.read_text())
    genes = fixture["genes"]
    matrices = {
        "wildtype": [fixture["wildtype"]],
        "knockouts": fixture["knockouts"],
        "knockdowns": fixture["knockdowns"],
    }
    files = {}
    for name, matrix in matrices.items():
        files[name + ".tsv"] = "\t".join(genes) + "\n" + "\n".join("\t".join(map(str, row)) for row in matrix) + "\n"
    positives = {tuple(edge) for edge in fixture["positive_edges"]}
    files["goldStandard.tsv"] = "".join(
        f"{a}\t{b}\t{int((a, b) in positives)}\n" for a in genes for b in genes if a != b
    )
    if transform:
        transform(files)
    with tarfile.open(path, "w:gz") as archive:
        for name, content in files.items():
            if name == missing:
                continue
            payload = content.encode()
            info = tarfile.TarInfo("lightlyProcessedDownloadedData/insilico_size10_1/" + name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def task(tmp_path, monkeypatch):
    def offline(*args, **kwargs):
        raise AssertionError("Unit tests must not use the network")

    monkeypatch.setattr(socket, "socket", offline)
    path = tmp_path / "fixture.tar.gz"
    digest = build_archive(path)
    return dream4_grn(path, tmp_path / "inputs", "development", (1,), dataset_sha256=digest)


def answer(task, value=0.0):
    scores = np.full((task.n, task.n), value)
    np.fill_diagonal(scores, 0)
    return {"scores": scores.tolist()}


def test_deterministic_load_and_no_gold_export(task, tmp_path):
    original = task.get_example()
    evidence_path = Path(original["input_path"])
    before = evidence_path.read_bytes()
    repeated = dream4_grn(tmp_path / "fixture.tar.gz", tmp_path / "inputs", "development", (1,), task.dataset_sha256)
    assert repeated.get_example() == original
    assert evidence_path.read_bytes() == before
    evidence = json.loads(before)
    assert set(evidence) == {"network", "genes", "wildtype", "knockouts", "knockdowns"}
    assert "answer" not in original
    assert len(task) == len(list(task.get_iterator())) == 1
    assert task.get_prompt_from_input(0) == original["prompt"]
    assert task.evaluate([answer(task)])["dataset_revision"].startswith("custom:")


def test_pin_and_default_checksum(tmp_path):
    assert DATA_REVISION == "d901bd60f4c621ad2df0244931c60a9f3780fed8"
    assert DATA_SHA256 == "1dede2b7578c552da2521c2ca31ebbcdaa04245eb7e54084e217d46731a871b2"
    assert f"/{DATA_REVISION}/" in DATA_URL
    path = tmp_path / "fixture.tar.gz"
    build_archive(path)
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        dream4_grn(path, tmp_path / "inputs")
    assert not (tmp_path / "inputs").exists()


def test_missing_archive(tmp_path):
    with pytest.raises(FileNotFoundError):
        dream4_grn(tmp_path / "absent.tar.gz", tmp_path / "inputs")


@pytest.mark.parametrize("missing", ["knockouts.tsv", "goldStandard.tsv", "wildtype.tsv", "knockdowns.tsv"])
def test_missing_member(tmp_path, missing):
    path = tmp_path / "fixture.tar.gz"
    digest = build_archive(path, missing=missing)
    with pytest.raises(ValueError, match="Missing required data"):
        dream4_grn(path, tmp_path / "inputs", "development", (1,), digest)


@pytest.mark.parametrize("bad", [None, "not json", "[]", {}, {"scores": [[1.0]]}, {"scores": [[float("nan")]]}])
def test_malformed_answers_remain_in_denominator(task, bad):
    result = task.evaluate([bad])
    assert result["macro_average_precision"] == 0
    assert result["valid_fraction"] == 0
    assert result["networks"][0]["error"]


def test_valid_answer_and_hand_calculated_metrics(task):
    # Two true edges at ranks 1 and 3: AP=(1+2/3)/2=5/6.
    a = answer(task)
    a["scores"][0][1] = 3.0
    a["scores"][0][2] = 2.0  # Indirect G1 -> G3 is NOT a direct edge.
    a["scores"][1][2] = 1.0
    result = task.evaluate([a])["networks"][0]
    assert result["average_precision"] == pytest.approx(5 / 6)
    assert result["auroc"] == pytest.approx(175 / 176)
    assert result["precision_at_n"] == pytest.approx(0.2)
    assert result["recall_at_n"] == 1
    assert result["indirect_fp_fraction"] == pytest.approx(1 / 8)
    assert result["candidate_edges"] == 90
    assert result["prevalence"] == pytest.approx(2 / 90)
    model = task.output_class().model_validate(a)
    assert task.evaluate([model])["macro_average_precision"] == pytest.approx(5 / 6)
    assert task.evaluate([json.dumps(a)])["macro_average_precision"] == pytest.approx(5 / 6)


def test_constant_scores_and_fractional_ties(task):
    result = task.evaluate([answer(task)])["networks"][0]
    assert result["average_precision"] == pytest.approx(2 / 90)
    assert result["auroc"] == 0.5
    assert result["precision_at_n"] == pytest.approx(2 / 90)
    assert result["recall_at_n"] == pytest.approx(10 / 90)


def test_perfect_and_reverse_orientation(task):
    a = answer(task)
    a["scores"][0][1] = a["scores"][1][2] = 1.0
    result = task.evaluate([a])["networks"][0]
    assert result["average_precision"] == result["auroc"] == 1
    reverse = {"scores": np.asarray(a["scores"]).T.tolist()}
    assert task.evaluate([reverse])["macro_average_precision"] < 0.1


@pytest.mark.parametrize("value", [float("inf"), float("-inf"), float("nan"), "1.0", True])
def test_invalid_score_types(task, value):
    a = answer(task)
    a["scores"][0][1] = value
    assert task.evaluate([a])["valid_fraction"] == 0


def test_extra_field_diagonal_and_wrong_counts(task):
    a = answer(task)
    a["extra"] = "forbidden"
    assert task.evaluate([a])["valid_fraction"] == 0
    a = answer(task)
    a["scores"][0][0] = 1.0
    assert task.evaluate([a])["valid_fraction"] == 0
    for responses in [[], [a, a], "not a batch"]:
        with pytest.raises(ValueError, match="exactly one response"):
            task.evaluate(responses)


def test_no_evidence_overwrite(task, tmp_path):
    p = Path(task.get_example()["input_path"])
    p.write_text("tampered")
    with pytest.raises(ValueError, match="Existing evidence differs"):
        dream4_grn(tmp_path / "fixture.tar.gz", tmp_path / "inputs", "development", (1,), task.dataset_sha256)
    assert p.read_text() == "tampered"


def test_archive_cannot_be_in_agent_mount(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        dream4_grn(tmp_path / "data.tar.gz", tmp_path)


def test_input_validation(task, tmp_path):
    for index in [-1, 1, None, "0"]:
        with pytest.raises(IndexError):
            task.get_example(index)
    for ids in [[], [1, 1], [6], [True], [1.0]]:
        with pytest.raises(ValueError, match="network_ids"):
            dream4_grn(tmp_path / "missing", tmp_path / "inputs", network_ids=ids)
    with pytest.raises(ValueError, match="split"):
        dream4_grn(tmp_path / "missing", tmp_path / "inputs", split="other")


@pytest.mark.parametrize("corruption", ["header", "nan", "gold_duplicate", "gold_missing", "gold_self"])
def test_corrupt_dataset_is_environment_error(tmp_path, corruption):
    def corrupt(files):
        if corruption == "header":
            files["wildtype.tsv"] = files["wildtype.tsv"].replace("G1", "Gwrong", 1)
        elif corruption == "nan":
            files["wildtype.tsv"] = files["wildtype.tsv"].replace("1.0", "nan", 1)
        elif corruption == "gold_duplicate":
            files["goldStandard.tsv"] += "G1\tG2\t1\n"
        elif corruption == "gold_missing":
            files["goldStandard.tsv"] = files["goldStandard.tsv"].split("\n", 1)[1]
        else:
            files["goldStandard.tsv"] += "G1\tG1\t1\n"

    path = tmp_path / "bad.tar.gz"
    digest = build_archive(path, transform=corrupt)
    with pytest.raises(ValueError):
        dream4_grn(path, tmp_path / "inputs", "development", (1,), digest)


def test_tie_weights_permutation_invariant():
    rng = np.random.default_rng(0)
    for _ in range(25):
        s = rng.integers(0, 4, 90).astype(float)
        p = rng.permutation(90)
        w = _top_weights(s, 10)
        assert w.sum() == pytest.approx(10)
        np.testing.assert_allclose(_top_weights(s[p], 10), w[p])


def test_fixture_is_not_accidentally_trivial(task):
    # A response-to-knockout ranking confuses a chain's indirect effect with a direct edge.
    evidence = json.loads(Path(task.get_example()["input_path"]).read_text())
    score = np.abs(np.asarray(evidence["knockouts"]) - evidence["wildtype"])
    np.fill_diagonal(score, 0)
    result = task.evaluate([{"scores": score.tolist()}])
    assert result["macro_average_precision"] == pytest.approx(2 / 3)
    assert result["macro_average_precision"] == average_precision_score(
        task._labels[0][~np.eye(10, dtype=bool)], score[~np.eye(10, dtype=bool)]
    )


def test_mcz_core_hand_computable():
    # Target 2 has observations [WT=1, KO1=0, KO2=0]: median=0, sample SD=1/sqrt(3).
    # Target 1 has [WT=1, KO1=0, KO2=1]: median=1, so the off-diagonal score is zero.
    e = {"wildtype": [1.0, 1.0], "knockouts": [[0.0, 0.0], [1.0, 0.0]]}
    np.testing.assert_allclose(median_corrected_zscore(e), np.zeros((2, 2)))
    # With a third KO that leaves target 2 unchanged: [1,0,0,1], median=0.5, SD=sqrt(1/3).
    e = {"wildtype": [1.0, 1.0, 1.0], "knockouts": [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]}
    assert median_corrected_zscore(e)[0, 1] == pytest.approx(0.5 / np.sqrt(1 / 3))


def test_mcz_constant_target_is_finite():
    e = {"wildtype": [1.0, 0.0], "knockouts": [[0.0, 0.0], [1.0, 0.0]]}
    assert np.isfinite(median_corrected_zscore(e)).all()
