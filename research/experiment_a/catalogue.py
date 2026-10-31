"""Exact 100-cell Experiment A catalogue and directly addressable PCG streams."""
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np

MASTER_SEED = 0x20261028B6A84FE0E1A0000000000001
RNG_VERSION = 'experiment-a-rng-v1'
PHASES = ('validation', 'benchmark', 'pilot', 'main')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def exact(coefficient, radical=1):
    coefficient = F(coefficient)
    factor = 1
    for k in range(2, math.isqrt(radical)+1):
        while radical % (k*k) == 0:
            radical //= k*k
            factor *= k
    coefficient *= factor
    if coefficient == 0:
        radical = 1
    return {'coefficient': [coefficient.numerator, coefficient.denominator], 'radical': radical}


def resolved(value):
    return float(F(*value['coefficient'])) * math.sqrt(value['radical'])


def catalogue():
    cells = {}
    for T in (100, 250, 500):
        drifts = [('zero', exact(0)), ('local_c0.5', exact(F(1, 2*T), T)),
                  ('local_c1', exact(F(1, T), T)), ('local_c2', exact(F(2, T), T)),
                  ('fixed_0.1', exact(F(1, 10))), ('fixed_0.5', exact(F(1, 2)))]
        roots = [('fixed_0.90', F(9, 10)), ('fixed_0.95', F(19, 20)),
                 ('fixed_0.99', F(99, 100)), ('local_c5', 1-F(5, T)),
                 ('local_c10', 1-F(10, T))]

        def add(rho, a, alias):
            law = dict(experiment='A_core', T=T, p=0, rho=[rho.numerator, rho.denominator],
                       intercept=a, sigma=[1, 1], y0=[0, 1], dtype='float64')
            key = digest(law)
            if key not in cells:
                cells[key] = dict(id=key, law=law, rho=float(rho), intercept=resolved(a),
                                  is_null=rho == 1, aliases=[])
            cells[key]['aliases'].append(alias)
        for label, a in drifts:
            add(F(1), a, dict(family='null', drift=label, root='unit'))
        for root_label, rho in roots:
            for drift_label, a in drifts:
                add(rho, a, dict(family='I', root=root_label, drift=drift_label))
            for m in (0, 1, 5):
                add(rho, exact((1-rho)*m), dict(family='M', root=root_label, mean=m))
    return sorted(cells.values(), key=lambda c: (c['law']['T'], c['rho'], c['intercept'], c['id']))


def innovations(cell, replicate, phase='pilot', master=MASTER_SEED):
    if phase not in PHASES or not isinstance(replicate, int) or isinstance(replicate, bool) or replicate < 0:
        raise ValueError('Only nonnegative replicate IDs in validation/benchmark/pilot/main are supported')
    if not 0 <= master < 2**128:
        raise ValueError('master seed must be an unsigned 128-bit integer')
    key = [RNG_VERSION, cell['id'], phase, replicate, 'innovations']
    words = np.frombuffer(hashlib.sha256(canonical(key).encode()).digest(), dtype='<u4').tolist()
    entropy = [(master >> (32*j)) & 0xffffffff for j in range(4)] + words
    rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(entropy)))
    return rng.standard_normal(cell['law']['T'], dtype=np.float64)


def generate(cell, z):
    T = cell['law']['T']
    z = np.asarray(z, dtype=np.float64)
    if z.shape != (T,) or not np.isfinite(z).all():
        raise ValueError('Exactly T finite standard innovations are required')
    y = np.empty(T+1, dtype=np.float64)
    y[0] = 0.
    for t in range(1, T+1):
        y[t] = cell['rho'] * y[t-1] + cell['intercept'] + z[t-1]
    return y


def array_hash(a):
    return hashlib.sha256(np.asarray(a, dtype='<f8').tobytes()).hexdigest()
