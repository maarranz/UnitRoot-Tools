import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from unitroot_tools import select_lags

SERIES = ["rw", "rw_negma", "ar07", "trend_ar08", "hetero_rw", "mic_boundary"]


@pytest.mark.parametrize("name", SERIES)
@pytest.mark.parametrize("trend", ["n", "c", "ct"])
@pytest.mark.parametrize("K", [4, 12])
@pytest.mark.parametrize("method", ["maic", "mbic", "aic", "bic"])
def test_reviewed_common_sample_diagnostics(stage2, name, trend, K, method):
    data, saved, _, _ = stage2
    expected = [r for r in saved if (r["series"], r["trend"], int(r["K"]), r["variant"]) ==
                (name, trend, K, "PQ_full")]
    result = select_lags(data[name], trend=trend, max_lags=K, method=method)
    assert len(expected) == len(result.candidates) == K + 1
    assert result.selection_nobs == len(data[name]) - K - 1
    assert (result.sample_start, result.sample_stop) == (K + 1, len(data[name]))
    scores = []
    for candidate, row in zip(result.candidates, expected, strict=True):
        p, n = int(row["p"]), int(row["n"])
        aic = math.log(float(row["sigma2"])) + 2*p/n
        bic = math.log(float(row["sigma2"])) + math.log(n)*p/n
        assert candidate.lags == p and candidate.nobs == n
        assert candidate.df_resid == n - p - 1
        # Absolute 1e-8 preserves the reviewed numerical convention.
        assert_allclose(
            [candidate.level_coefficient, candidate.variance, candidate.ssr,
             candidate.level_sum_squares, candidate.tau, candidate.maic, candidate.mbic,
             candidate.aic, candidate.bic],
            [float(row["beta"]), float(row["sigma2"]), float(row["sigma2"])*n,
             float(row["level_ss"]), float(row["tau"]), float(row["MAIC"]),
             float(row["MBIC"]), aic, bic], rtol=0, atol=1e-8,
        )
        score = {"aic": aic, "bic": bic, "maic": float(row["MAIC"]), "mbic": float(row["MBIC"])}[method]
        scores.append((score, p))
        assert_allclose(candidate.criterion_value, score, rtol=0, atol=1e-8)
    ordered = sorted(scores)
    assert result.selected_lags == ordered[0][1]
    assert_allclose(result.criterion_value, ordered[0][0], rtol=0, atol=1e-8)
    assert_allclose(result.runner_up_gap, ordered[1][0]-ordered[0][0], rtol=0, atol=1e-8)


@pytest.mark.parametrize("name", SERIES)
@pytest.mark.parametrize("trend", ["c", "ct"])
@pytest.mark.parametrize("K", [4, 12])
def test_saved_development_dfgls_and_r_discrepancy(stage2, name, trend, K):
    data, _, development, r_results = stage2
    result = select_lags(data[name], trend=trend, max_lags=K)
    native = [r for r in development if (r["series"], r["trend"], r["test"], r["variant"]) ==
              (name, trend, "MAIC_DFGLS", f"K{K}_full")][0]
    assert result.selected_lags == int(native["lags"])
    assert_allclose(result.criterion_value, float(native["criterion_value"]), rtol=0, atol=1e-8)
    r = [r for r in r_results if (r["series"], r["trend"], r["test"], r["variant"]) ==
         (name, trend, "MAIC", f"two_step_K{K}_unscaled")][0]
    if (name, trend, K) == ("mic_boundary", "c", 4):
        assert result.selected_lags == 3 and int(r["lags"]) == 2
    else:
        assert result.selected_lags == int(float(r["lags"]))


def test_boundary_penalty_and_distinct_conventional_choice(stage2):
    y = stage2[0]["mic_boundary"]
    modified = select_lags(y, max_lags=4)
    conventional = select_lags(y, max_lags=4, method="aic")
    assert modified.selected_lags == 3
    assert conventional.selected_lags != modified.selected_lags
    v = y-y.mean()
    assert_allclose(modified.candidates[0].level_sum_squares, np.dot(v[4:-1], v[4:-1]), rtol=0, atol=1e-12)
    for row in modified.candidates:
        missing_penalty = 2/59 * row.level_coefficient**2 * v[4]**2 / row.variance
        saved = [r for r in stage2[1] if (r["series"],r["trend"],r["K"],r["p"],r["variant"]) ==
                 ("mic_boundary","c","4",str(row.lags),"PQ_bootUR")][0]
        assert_allclose(row.maic-float(saved["MAIC"]), missing_penalty, rtol=0, atol=1e-12)


@pytest.mark.parametrize("name", SERIES)
@pytest.mark.parametrize("trend", ["n", "c", "ct"])
@pytest.mark.parametrize("K", [4, 12])
@pytest.mark.parametrize("method", ["maic", "mbic"])
def test_reviewed_bootur_modified_choices(stage2, name, trend, K, method):
    data, _, _, saved = stage2
    r = select_lags(data[name], trend=trend, max_lags=K, method=method)
    expected = [row for row in saved if
                (row["series"], row["trend"], row["test"], row["variant"]) ==
                (name, trend, method.upper(), f"two_step_K{K}_unscaled")][0]
    if (name, trend, K, method) == ("mic_boundary", "c", 4, "maic"):
        assert (r.selected_lags, int(expected["lags"])) == (3, 2)
    else:
        assert r.selected_lags == int(expected["lags"])
