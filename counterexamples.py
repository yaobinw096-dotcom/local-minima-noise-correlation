"""Check exact counterexamples and correlation shape, without Monte Carlo."""
from fractions import Fraction as F
from pathlib import Path
import json

from verify_threshold import cycle, graph, poly_equal
from verify_noise import exact_noise, noise_formula


def normalize(coefficients):
    variance = sum(coefficients)
    assert variance > 0
    return [coefficient / variance for coefficient in coefficients]


def value(coefficients, p):
    return sum(coefficient * p**k for k, coefficient in enumerate(coefficients))


def circular_ladder(length):
    edges = [(i, (i + 1) % length) for i in range(length)]
    edges += [(length + i, length + (i + 1) % length) for i in range(length)]
    edges += [(i, length + i) for i in range(length)]
    return graph(2 * length, edges)


def generalized_petersen(length):
    edges = [(i, (i + 1) % length) for i in range(length)]
    edges += [(length + i, length + (i + 2) % length) for i in range(length)]
    edges += [(i, length + i) for i in range(length)]
    return graph(2 * length, edges)


if __name__ == '__main__':
    records = []
    expected = {
        'C4': [F(0), F(0), F(3, 5), F(2, 5)],
        'C5': [F(0), F(0), F(1, 2), F(1, 2)],
    }
    for name, n in [('C4', 4), ('C5', 5)]:
        actual = exact_noise(cycle(n))
        correlation = normalize(actual)
        assert poly_equal(correlation, expected[name]), name
        records.append({
            'graph': name,
            'method': 'exact linear-extension count',
            'variance': str(sum(actual)),
            'correlation_coefficients': [str(x) for x in correlation],
            'correlation_at_half': str(value(correlation, F(1, 2))),
            'passed': True,
        })

    for length in (9, 10, 11, 20):
        for name, adj, prediction in [
            ('circular_ladder', circular_ladder(length),
             [F(0), F(0), F(11, 20), F(9, 25), F(9, 100)]),
            ('generalized_petersen', generalized_petersen(length),
             [F(0), F(0), F(1, 2), F(2, 5), F(1, 10)]),
        ]:
            polynomial, moments, d = noise_formula(adj)
            correlation = normalize(polynomial)
            assert poly_equal(correlation, prediction), (name, length)
            assert correlation[0] == correlation[1] == 0
            assert all(x >= 0 for x in correlation)
            assert sum(correlation) == 1
            assert correlation[2] > 0 and correlation[-1] > 0
            for p in (F(1, 10), F(1, 2), F(9, 10)):
                assert 0 < value(correlation, p) < p*p
                ratio_derivative = sum(
                    (k-2)*x*p**(k-3)
                    for k, x in enumerate(correlation) if k >= 3
                )
                assert ratio_derivative > 0
            records.append({
                'graph': name,
                'length': length,
                'n': len(adj),
                'd': d,
                'method': 'exact graph invariants and covariance identity',
                'S2': moments[2],
                'correlation_coefficients': [str(x) for x in correlation],
                'passed': True,
            })

    destination = Path(__file__).resolve().parent / 'counterexamples.json'
    destination.write_text(json.dumps(records, indent=2), encoding='utf-8')
    print('All 10 exact counterexample checks passed.')
    print('C4 Corr(1/2)=1/5; C5 Corr(1/2)=3/16.')
    print('Their difference is 1/80 and does not shrink with repeated components.')
