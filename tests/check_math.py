"""Independent complex-rotation and Taylor-series neutrino checks. Standard library only."""

import cmath
import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PHASE = 2.5338653580781976
POTENTIAL = math.sqrt(2) * 8.961877245622253e-38 * 6.02214076e23
RUNTIME_POINTS = [
    [0.59, 0.15, 0.78, 1.2, 7.5e-5, 2.5e-3, 2.5, 295.0, 1.4, 0.01],
    [0.61, 0.14, 0.81, -0.7, 7.4e-5, -2.4e-3, 1.1, 800.0, 1.9, 0.02],
    [0.0, 0.0, 0.0, 0.0, 7.5e-5, 2.5e-3, 2.5, 295.0, 0.0, 0.01],
]


def multiply(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)
    ]


def adjoint(a):
    return [[complex(a[j][i]).conjugate() for j in range(3)] for i in range(3)]


def diagonal(xs):
    return [[xs[i] if i == j else 0 for j in range(3)] for i in range(3)]


def pmns(t12, t13, t23, delta):
    c12, c13, c23 = map(math.cos, [t12, t13, t23])
    s12, s13, s23 = map(math.sin, [t12, t13, t23])
    r12 = [[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]]
    r13 = [
        [c13, 0, s13 * cmath.exp(-1j * delta)],
        [0, 1, 0],
        [-s13 * cmath.exp(1j * delta), 0, c13],
    ]
    r23 = [[1, 0, 0], [0, c23, s23], [0, -s23, c23]]
    return multiply(r23, multiply(r13, r12))


def amplitude(x, density=None, anti=False):
    t12, t13, t23, delta, dm21, dm31, energy, length = x[:8]
    u = pmns(t12, t13, t23, -delta if anti else delta)
    if density is None:
        phases = [cmath.exp(-1j * PHASE * length * d / energy) for d in (0, dm21, dm31)]
        return multiply(multiply(u, diagonal(phases)), adjoint(u))
    h = multiply(multiply(u, diagonal([0, dm21, dm31])), adjoint(u))
    h[0][0] += (-1 if anti else 1) * 2e9 * energy * POTENTIAL * density
    # Independent 3x3 complex Taylor expansion, rather than a real-block Padé solve.
    a = [[-1j * PHASE * length * h[i][j] / energy for j in range(3)] for i in range(3)]
    result = diagonal([1, 1, 1])
    term = diagonal([1, 1, 1])
    for order in range(1, 100):
        term = [[v / order for v in row] for row in multiply(term, a)]
        result = [[result[i][j] + term[i][j] for j in range(3)] for i in range(3)]
        if max(abs(v) for row in term for v in row) < 1e-17:
            break
    return result


def probabilities(s):
    return [abs(v) ** 2 for row in s for v in row]


def runtime_values(x):
    t12, t13, _, delta, dm21, dm31, _, length, density, lightest = x
    vacuum = probabilities(amplitude(x))
    matter = probabilities(amplitude(x, density))
    first, second = list(x), list(x)
    first[7], second[7] = 0.25 * length, 0.75 * length
    layers = probabilities(
        multiply(amplitude(second, 0.5 * density), amplitude(first, density))
    )
    masses = [math.sqrt(lightest**2 + max(0, -dm31) + d) for d in (0, dm21, dm31)]
    row = pmns(t12, t13, 0, delta)[0]
    beta = math.sqrt(sum(abs(u) ** 2 * m * m for u, m in zip(row, masses)))
    majorana = abs(
        sum(
            u * u * m * cmath.exp(1j * p)
            for u, m, p in zip(row, masses, [0, 0.8, -0.6])
        )
    )
    baseline = math.sqrt((6371 + 20) ** 2 - 6371**2 * (1 - 0.6**2)) + 6371 * 0.6
    half_inner, half_outer = math.sqrt(3), math.sqrt(8) - math.sqrt(3)
    return (
        vacuum
        + matter
        + layers
        + masses
        + [beta, majorana, baseline, half_outer, half_inner, half_inner, half_outer]
    )


def posterior(theta13):
    if not 0.1 <= theta13 <= 0.2:
        return -math.inf
    score = math.log(10)
    for energy, observed in zip([0.45, 0.6, 0.9, 1.2], [38, 48, 32, 21]):
        point = [0.59, theta13, 0.78, 1.2, 7.5e-5, 0.0025, energy, 295]
        probability = abs(amplitude(point, 1.4)[0][1]) ** 2
        score += (
            math.lgamma(1001) - math.lgamma(observed + 1) - math.lgamma(1001 - observed)
        )
        score += observed * math.log(probability) + (1000 - observed) * math.log1p(
            -probability
        )
    return score


def decode(value):
    if isinstance(value, dict):
        if "data" in value:
            data = [decode(x) for x in value["data"].values()]
            return data[0] if value.get("shape") == [] else data
        return {key: decode(x) for key, x in value.items()}
    if isinstance(value, list):
        return [decode(x) for x in value]
    return float(value) if isinstance(value, str) else value


def flatten(value):
    return (
        [item for v in value for item in flatten(v)]
        if isinstance(value, list)
        else [value]
    )


def evaluate(request):
    result = subprocess.run(
        ["node", str(ROOT / "tests/evaluate.cjs")],
        input=json.dumps(request),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stderr)
    return decode(json.loads(result.stdout))


def near(actual, expected, tolerance=3e-10):
    assert math.isfinite(actual) and math.isclose(
        actual, expected, rel_tol=tolerance, abs_tol=tolerance
    ), (actual, expected)


def main():
    source = (ROOT / "tests/runtime.flatppl").read_text()
    for point in RUNTIME_POINTS:
        actual = evaluate(
            {
                "source": source.replace(
                    "x = elementof(cartpow(reals, 10))", f"x = {point}"
                ),
                "bindings": [
                    "v",
                    "m",
                    "l",
                    "mass",
                    "beta",
                    "majorana",
                    "baseline",
                    "path",
                ],
            }
        )
        actual[-1] = actual[-1]["columns"]["length"]
        for value, reference in zip(
            flatten(actual), runtime_values(point), strict=True
        ):
            near(value, reference)
    # Strong asymmetric CP test distinguishes P(mu<-e) from P(e<-mu).
    p = "a.parameters(pi/4, pi/4, pi/4, pi/2, 1, 2)"
    length = 2 * math.pi / (3 * PHASE)
    values = evaluate(
        [
            f"a.vacuum_probabilities(1,{length},{p},false)[2,1]",
            f"a.vacuum_probabilities(1,{length},{p},false)[1,2]",
            f"a.vacuum_probabilities(1,{length},{p},true)[2,1]",
            "a.atmospheric_baseline(1,6371,1e-10)",
            "a.atmospheric_baseline(-1,6371,20)",
            "a.atmospheric_baseline(0,6371,0)",
        ]
    )
    expected = [
        0.9280293267718458,
        0.009470673228154102,
        0.009470673228154102,
        1e-10,
        12762,
        0,
    ]
    for value, reference in zip(values, expected, strict=True):
        near(value, reference, 2e-13)
    # The matter potential is resonant for neutrinos, and has the opposite sign for anti-neutrinos.
    d, energy, length = 7.5e-5, 1.0, 2000.0
    density = (d / 2) / (2e9 * energy * POTENTIAL)
    p = f"a.parameters(pi/6,0,0,0,{d},0.0025)"
    values = evaluate(
        [
            f"a.matter_probabilities(1,{length},{p},{density},{anti})[2,1]"
            for anti in ("false", "true")
        ]
    )
    expected = [
        math.sin(PHASE * length * math.sqrt(3) * d / 4) ** 2,
        (3 / 7) * math.sin(PHASE * length * math.sqrt(7) * d / 4) ** 2,
    ]
    for value, reference in zip(values, expected, strict=True):
        near(value, reference)
    point = RUNTIME_POINTS[1]
    p = "a.parameters(" + ",".join(map(str, point[:6])) + ")"
    first, second = list(point), list(point)
    first[7], second[7] = 0.25 * point[7], 0.75 * point[7]
    reference = multiply(amplitude(second, 0.5 * point[8]), amplitude(first, point[8]))
    expression = f"a.layered_amplitude({point[6]},{p},[{first[7]},{second[7]}],[{point[8]},{0.5 * point[8]}],false)"
    result = evaluate([expression])[0]
    for part in ("re", "im"):
        expected = [
            getattr(v, "real" if part == "re" else "imag")
            for row in reference
            for v in row
        ]
        for value, target in zip(result[part], expected, strict=True):
            near(value, target)
    source = (ROOT / "examples/playground.flatppl").read_text()
    for theta13 in (0.12, 0.15, 0.18):
        actual = evaluate(
            {
                "source": source
                + f"\nscore = logdensityof(Pi_post, (theta13 = {theta13}))\n",
                "bindings": ["score"],
            }
        )[0]
        near(actual, posterior(theta13), 2e-10)
    print(
        "Independent complex propagation, CP, matter sign, mass and geometry checks passed."
    )


if __name__ == "__main__":
    main()
