"""Tests for analysis.py — regression, correlation, normalization."""

import numpy as np
import pandas as pd
import pytest

from src.analysis import (
    linear_regression,
    compute_pearson_correlation,
    normalize_minmax,
    normalize_zscore,
    pathway_selectivity_score,
)


# ── linear_regression ─────────────────────────────────────────────────────────

def test_linear_regression_perfect_fit():
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]  # y = 2x
    result = linear_regression(x, y)
    assert abs(result["slope"] - 2.0) < 1e-10
    assert abs(result["intercept"]) < 1e-10
    assert abs(result["r2"] - 1.0) < 1e-10


def test_linear_regression_returns_required_keys():
    result = linear_regression([1, 2, 3], [1, 2, 3])
    for key in ("slope", "intercept", "r2", "p_value", "std_err", "equation"):
        assert key in result


def test_linear_regression_equation_string():
    result = linear_regression([0, 1], [0, 1])
    assert "R²" in result["equation"]
    assert "y =" in result["equation"]


def test_linear_regression_r2_range():
    x = np.arange(10)
    y = np.arange(10) + np.random.default_rng(0).normal(0, 0.5, 10)
    result = linear_regression(x, y)
    assert 0.0 <= result["r2"] <= 1.0


# ── compute_pearson_correlation ───────────────────────────────────────────────

def test_pearson_correlation_diagonal_is_one():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [3, 2, 1], "c": [1, 3, 2]})
    corr = compute_pearson_correlation(df)
    for col in corr.columns:
        assert abs(corr.loc[col, col] - 1.0) < 1e-10


def test_pearson_ignores_non_numeric():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "label": ["x", "y", "z"]})
    corr = compute_pearson_correlation(df)
    assert "label" not in corr.columns
    assert "a" in corr.columns


# ── normalize_minmax ──────────────────────────────────────────────────────────

def test_minmax_range():
    df = pd.DataFrame({"x": [0.0, 5.0, 10.0], "y": [-1.0, 0.0, 1.0]})
    normed = normalize_minmax(df)
    assert abs(normed["x"].min()) < 1e-10
    assert abs(normed["x"].max() - 1.0) < 1e-10


def test_minmax_ignores_non_numeric():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "label": ["a", "b", "c"]})
    normed = normalize_minmax(df)
    assert "label" not in normed.columns


# ── normalize_zscore ──────────────────────────────────────────────────────────

def test_zscore_mean_and_std():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0, 5.0]})
    normed = normalize_zscore(df)
    assert abs(normed["x"].mean()) < 1e-10
    assert abs(normed["x"].std() - 1.0) < 1e-1   # ddof=1, close enough


# ── pathway_selectivity_score ─────────────────────────────────────────────────

def test_selectivity_score_ratio():
    barriers = pd.DataFrame(
        {"sys_A": [0.8, 0.4], "sys_B": [1.0, 0.5]},
        index=["nrr_barrier", "her_barrier"],
    )
    scores = pathway_selectivity_score(barriers)
    assert abs(scores["sys_A"] - 2.0) < 1e-10   # 0.8 / 0.4
    assert abs(scores["sys_B"] - 2.0) < 1e-10   # 1.0 / 0.5


def test_selectivity_score_missing_row_raises():
    barriers = pd.DataFrame({"sys_A": [0.8]}, index=["nrr_barrier"])
    with pytest.raises(ValueError, match="her_barrier"):
        pathway_selectivity_score(barriers)
