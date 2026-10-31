"""Strict independent score checks complement the historical raw-fit budgets."""
import importlib.util
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('verify_deterministic_scaled',Path(__file__).resolve().parents[1]/'validation/verify_deterministic_scaled.py')
validator=importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)

@pytest.mark.parametrize('key',['y','positive','negative'])
@pytest.mark.parametrize('p',[0,1,4])
@pytest.mark.parametrize('criterion',['aic','bic'])
def test_saved_scaled_scores_against_decimal(key,p,criterion):
    validator.verify_case(key,p,criterion)
