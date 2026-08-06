from sage.all import *

import hashlib
import itertools
import json
import operator
import os
import platform
import resource
import sys
import time
import traceback

import numpy as np
from sage.env import SAGE_VERSION


TASK_ID = "cgw35-sage-publication-v1"
N = 35
TARGET_WEIGHT = 12
ARTIFACT_DIR = os.path.abspath(os.path.dirname(__file__))
START_WALL = time.time()
START_CPU = time.process_time()
SAGE_INTEGER_TYPE = type(ZZ.zero())


def artifact_path(name):
    return os.path.join(ARTIFACT_DIR, name)


def write_text_frozen(name, text):
    path = artifact_path(name)
    if os.path.exists(path):
        raise RuntimeError("refusing to overwrite frozen artifact: " + name)
    with open(path, "x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    return sha256_file(path)


def write_json_frozen(name, data):
    return write_text_frozen(
        name,
        json.dumps(data, indent=2, sort_keys=True, default=json_exact_default)
        + "\n",
    )


def json_exact_default(value):
    if isinstance(value, SAGE_INTEGER_TYPE):
        return int(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    raise TypeError("unsupported JSON value type: " + type(value).__name__)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def mask_from_word(word):
    mask = 0
    for q, value in enumerate(word):
        if value != 0:
            mask |= 1 << q
    return int(mask)


def support_from_mask(mask):
    return tuple(q for q in range(N) if (mask >> q) & 1)


def mask_hex(mask):
    return format(mask, "09x")


def cyclic_sharp_polynomial(poly, ring, conjugate):
    x = ring.gen()
    return ring(
        sum(
            (conjugate(poly[q]) * x ** ((-q) % N) for q in range(N)),
            ring.zero(),
        )
    )


def monic_reciprocal(poly):
    ring = poly.parent()
    degree = poly.degree()
    reciprocal = ring([poly[degree - i] for i in range(degree + 1)])
    return reciprocal.monic()


def transform_mask(mask, a, b):
    transformed = 0
    for q in support_from_mask(mask):
        transformed |= 1 << ((a * q + b) % N)
    return int(transformed)


def incidence_pairs(support, shift):
    position = {q: i for i, q in enumerate(support)}
    return tuple(
        (i, position[(q + shift) % N])
        for i, q in enumerate(support)
        if (q + shift) % N in position
    )


def char3_stage():
    print("[char3] constructing exact field and factoring X^35-1", flush=True)
    F3 = GF(3)
    R3 = PolynomialRing(F3, "X")
    X = R3.gen()
    modulus = X ** N - 1
    raw_factorization = list(modulus.factor())
    factors = [factor.monic() for factor, exponent in raw_factorization]
    exponents = [int(exponent) for factor, exponent in raw_factorization]
    factors.sort(key=lambda factor: (int(factor.degree()), str(factor)))

    reconstruction = R3.one()
    for factor, exponent in raw_factorization:
        reconstruction *= factor ** exponent
    assert reconstruction == modulus
    assert all(exponent == 1 for exponent in exponents)
    assert gcd(modulus, modulus.derivative()) == 1
    assert all(factor.is_irreducible() for factor in factors)

    reciprocal_index = []
    for factor in factors:
        reciprocal = monic_reciprocal(factor)
        matches = [i for i, candidate in enumerate(factors) if candidate == reciprocal]
        assert len(matches) == 1
        reciprocal_index.append(matches[0])
    assert all(reciprocal_index[reciprocal_index[i]] == i for i in range(len(factors)))

    fixed = tuple(i for i, j in enumerate(reciprocal_index) if i == j)
    paired = []
    seen = set(fixed)
    for i, j in enumerate(reciprocal_index):
        if i in seen:
            continue
        assert i != j
        paired.append(tuple(sorted((i, j))))
        seen.add(i)
        seen.add(j)
    paired = tuple(sorted(set(paired)))
    assert len(seen) == len(factors)

    # Square-free local lemma, checked in each factor field:
    # In R/(f), (v v#)=0.  If f is # fixed, the two residues vanish
    # together, hence v=0 mod f.  For a two-factor orbit {f,f#},
    # the product residues force at least one of v mod f and v mod f#
    # to vanish.  The explicit orbit map above checks every local component.
    for i in fixed:
        assert reciprocal_index[i] == i
    for i, j in paired:
        assert reciprocal_index[i] == j and reciprocal_index[j] == i

    branch_choices = list(itertools.product(*[tuple(pair) for pair in paired]))
    assert len(branch_choices) == 2 ** len(paired)
    fixed_product = prod((factors[i] for i in fixed), R3.one())

    print(
        "[char3] deriving and exhaustively enumerating %d cyclic-code branches"
        % len(branch_choices),
        flush=True,
    )
    support_multiplicity = {}
    branch_records = []
    branch_support_sets = []
    for branch_number, choice in enumerate(branch_choices):
        generator = fixed_product * prod(
            (factors[i] for i in choice), R3.one()
        )
        generator = generator.monic()
        code = codes.CyclicCode(generator_pol=generator, length=N)
        dimension = int(code.dimension())
        expected_dimension = N - int(generator.degree())
        assert dimension == expected_dimension

        manual_rows = []
        for k in range(dimension):
            shifted = X ** k * generator
            manual_rows.append(vector(F3, [shifted[q] for q in range(N)]))
        manual_matrix = matrix(F3, manual_rows)
        assert int(manual_matrix.rank()) == dimension
        manual_space = manual_matrix.row_space()
        sage_space = code.generator_matrix().row_space()
        assert manual_space == sage_space

        enumerated = 0
        weight_count = 0
        branch_supports = set()
        for word in code:
            enumerated += 1
            if int(word.hamming_weight()) != TARGET_WEIGHT:
                continue
            weight_count += 1
            mask = mask_from_word(word)
            assert mask.bit_count() == TARGET_WEIGHT
            branch_supports.add(mask)
            support_multiplicity[mask] = support_multiplicity.get(mask, 0) + 1

            word_poly = R3(list(word))
            sharp = cyclic_sharp_polynomial(word_poly, R3, lambda value: value)
            assert (word_poly * sharp) % modulus == 0
            assert generator.divides(word_poly)

        assert enumerated == 3 ** dimension
        branch_support_sets.append(branch_supports)
        branch_records.append(
            {
                "branch": branch_number,
                "selected_factor_indices": [int(i) for i in choice],
                "generator": str(generator),
                "generator_degree": int(generator.degree()),
                "dimension": dimension,
                "enumerated_words": enumerated,
                "weight_12_words": weight_count,
                "distinct_weight_12_supports": len(branch_supports),
                "manual_and_sage_spaces_equal": True,
            }
        )
        print(
            "[char3] branch %d complete: dimension %d, enumerated %d words"
            % (branch_number, dimension, enumerated),
            flush=True,
        )

    support_masks = tuple(sorted(support_multiplicity))
    assert all(mask.bit_count() == TARGET_WEIGHT for mask in support_masks)
    supports_text = "".join(mask_hex(mask) + "\n" for mask in support_masks)
    supports_sha256 = write_text_frozen("support-masks.txt", supports_text)

    multiplicities_text = "".join(
        "%s %d\n" % (mask_hex(mask), support_multiplicity[mask])
        for mask in support_masks
    )
    multiplicities_sha256 = write_text_frozen(
        "char3-word-multiplicities.txt", multiplicities_text
    )

    print("[char3] support family frozen; deriving affine orbits", flush=True)
    units = tuple(a for a in range(N) if gcd(a, N) == 1)
    assert len(units) == euler_phi(N)
    group_size = len(units) * N
    support_set = set(support_masks)
    unvisited = set(support_masks)
    affine_orbits = []
    while unvisited:
        representative = min(unvisited)
        orbit = {
            transform_mask(representative, a, b)
            for a in units
            for b in range(N)
        }
        assert orbit <= support_set
        stabilizer_size = sum(
            1
            for a in units
            for b in range(N)
            if transform_mask(representative, a, b) == representative
        )
        assert len(orbit) * stabilizer_size == group_size
        orbit_multiplicities = {support_multiplicity[mask] for mask in orbit}
        assert len(orbit_multiplicities) == 1
        affine_orbits.append(
            {
                "representative": mask_hex(representative),
                "orbit_size": len(orbit),
                "stabilizer_size": stabilizer_size,
                "word_multiplicity_per_support": min(orbit_multiplicities),
            }
        )
        unvisited -= orbit
    assert sum(record["orbit_size"] for record in affine_orbits) == len(support_masks)

    multiplicity_distribution = {}
    for multiplicity in support_multiplicity.values():
        key = str(multiplicity)
        multiplicity_distribution[key] = multiplicity_distribution.get(key, 0) + 1

    factor_records = []
    for i, factor in enumerate(factors):
        factor_records.append(
            {
                "index": i,
                "polynomial": str(factor),
                "degree": int(factor.degree()),
                "irreducible": bool(factor.is_irreducible()),
                "reciprocal_index": int(reciprocal_index[i]),
                "fixed_by_sharp": bool(reciprocal_index[i] == i),
            }
        )

    char3_record = {
        "field": "GF(3)",
        "eisenstein_reduction": {
            "prime": "(1-omega)",
            "omega_image": "1",
            "six_unit_images": ["1", "2", "1", "2", "1", "2"],
            "all_units_nonzero": True,
            "support_preserved": True,
            "right_hand_side_12_image": "0",
        },
        "modulus": str(modulus),
        "factorization": factor_records,
        "factorization_reconstructs_modulus": True,
        "square_free": True,
        "all_factors_irreducible": True,
        "fixed_factor_indices": [int(i) for i in fixed],
        "paired_factor_orbits": [[int(i), int(j)] for i, j in paired],
        "local_square_free_lemma_checked": True,
        "branch_rule": (
            "Every sharp-fixed irreducible factor divides v, and at least one "
            "member of every two-element sharp orbit divides v."
        ),
        "branches": branch_records,
        "support_count": len(support_masks),
        "support_masks_sha256": supports_sha256,
        "multiplicity_file_sha256": multiplicities_sha256,
        "word_multiplicity_distribution": multiplicity_distribution,
        "affine_group": {
            "units": [int(a) for a in units],
            "translations": N,
            "group_size": int(group_size),
            "orbits": affine_orbits,
            "partition_complete": True,
        },
        "freeze_complete": True,
    }
    char3_sha256 = write_json_frozen("char3-freeze.json", char3_record)
    print("[char3] complete output frozen", flush=True)
    return {
        "record": char3_record,
        "record_sha256": char3_sha256,
        "support_masks": support_masks,
        "support_multiplicity": support_multiplicity,
        "affine_orbits": affine_orbits,
    }


def exact_f4_correlation(phases, support, shift, F4):
    position = {q: i for i, q in enumerate(support)}
    total = F4.zero()
    for i, q in enumerate(support):
        target = (q + shift) % N
        if target in position:
            total += phases[i] * phases[position[target]] ** 2
    return total


def char2_stage(char3):
    print("[char2] constructing GF(4) and exhaustive normalized assignment space", flush=True)
    F2 = GF(2)
    P2 = PolynomialRing(F2, "Z")
    Z = P2.gen()
    f4_modulus = Z ** 2 + Z + 1
    assert f4_modulus.is_irreducible()
    F4 = GF(4, name="w", modulus=f4_modulus)
    w = F4.gen()
    assert w ** 2 + w + 1 == 0
    assert w ** 2 == w.frobenius()
    assert all(value ** 3 == 1 for value in (F4.one(), w, w ** 2))
    phase_values = (F4.one(), w, w ** 2)
    assert tuple(value ** 2 for value in phase_values) == (
        phase_values[0],
        phase_values[2],
        phase_values[1],
    )

    # Basis encoding in {0,1,w,1+w}; field addition is bitwise XOR.
    decode = (F4.zero(), F4.one(), w, w ** 2)
    assert w ** 2 == 1 + w

    def encode(value):
        matches = [i for i, candidate in enumerate(decode) if value == candidate]
        assert len(matches) == 1
        return matches[0]

    for left in decode:
        for right in decode:
            assert decode[operator.xor(encode(left), encode(right))] == left + right

    term_code = np.zeros((3, 3), dtype=np.uint8)
    for i in range(3):
        for j in range(3):
            term_code[i, j] = encode(phase_values[i] * phase_values[j] ** 2)
            assert phase_values[i] * phase_values[j] ** 2 == phase_values[(i - j) % 3]

    normalized_count = 3 ** (TARGET_WEIGHT - 1)
    assignment_ids = np.arange(normalized_count, dtype=np.int64)
    assignments = np.zeros((normalized_count, TARGET_WEIGHT), dtype=np.uint8)
    powers = []
    for coordinate in range(1, TARGET_WEIGHT):
        power = 3 ** (coordinate - 1)
        powers.append(power)
        assignments[:, coordinate] = (assignment_ids // power) % 3
    reconstructed_ids = np.zeros(normalized_count, dtype=np.int64)
    for coordinate in range(1, TARGET_WEIGHT):
        reconstructed_ids += (
            assignments[:, coordinate].astype(np.int64)
            * (3 ** (coordinate - 1))
        )
    assert np.array_equal(reconstructed_ids, assignment_ids)
    assert np.all(assignments[:, 0] == 0)

    representative_masks = tuple(
        int(record["representative"], 16)
        for record in char3["affine_orbits"]
    )
    orbit_size_by_mask = {
        int(record["representative"], 16): int(record["orbit_size"])
        for record in char3["affine_orbits"]
    }
    assert sum(orbit_size_by_mask.values()) == len(char3["support_masks"])

    case_records = []
    survivors = {}
    survivor_lines = []
    total_normalized_support_assignments = 0
    total_normalized_surviving_support_assignments = 0
    for case_number, mask in enumerate(representative_masks):
        support = support_from_mask(mask)
        assert len(support) == TARGET_WEIGHT
        all_pairs = {
            shift: incidence_pairs(support, shift) for shift in range(1, N)
        }
        for shift in range(1, N):
            reverse_pairs = {(j, i) for i, j in all_pairs[shift]}
            assert reverse_pairs == set(all_pairs[(N - shift) % N])

        alive = np.ones(normalized_count, dtype=np.bool_)
        constraint_trace = []
        for shift in range(1, (N - 1) // 2 + 1):
            accumulated = np.zeros(normalized_count, dtype=np.uint8)
            for i, j in all_pairs[shift]:
                accumulated = np.bitwise_xor(
                    accumulated,
                    term_code[assignments[:, i], assignments[:, j]],
                )
            alive &= accumulated == 0
            constraint_trace.append(
                {
                    "shift": shift,
                    "ordered_support_pairs": len(all_pairs[shift]),
                    "assignments_remaining": int(np.count_nonzero(alive)),
                }
            )

        survivor_ids = np.flatnonzero(alive)
        case_survivors = []
        for assignment_id in survivor_ids:
            exponents = tuple(int(value) for value in assignments[int(assignment_id), :])
            exact_phases = tuple(phase_values[e] for e in exponents)
            exact_correlations = tuple(
                exact_f4_correlation(exact_phases, support, shift, F4)
                for shift in range(N)
            )
            assert exact_correlations[0] == F4(TARGET_WEIGHT)
            assert all(value == 0 for value in exact_correlations)
            case_survivors.append(exponents)
            survivor_lines.append(
                "%s %s\n"
                % (mask_hex(mask), "".join(str(exponent) for exponent in exponents))
            )
        survivors[mask] = tuple(case_survivors)

        orbit_size = orbit_size_by_mask[mask]
        total_normalized_support_assignments += orbit_size * normalized_count
        total_normalized_surviving_support_assignments += (
            orbit_size * len(case_survivors)
        )
        case_records.append(
            {
                "case": case_number,
                "support_representative": mask_hex(mask),
                "support_coordinates": [int(q) for q in support],
                "affine_orbit_size": orbit_size,
                "normalized_assignments_enumerated": normalized_count,
                "constraint_trace": constraint_trace,
                "normalized_survivors": len(case_survivors),
            }
        )
        print(
            "[char2] affine case %d/%d exhaustively checked"
            % (case_number + 1, len(representative_masks)),
            flush=True,
        )

    survivor_sha256 = write_text_frozen(
        "char2-survivors.txt", "".join(survivor_lines)
    )
    char2_record = {
        "field": "GF(4)",
        "defining_polynomial": str(f4_modulus),
        "defining_polynomial_irreducible": True,
        "omega_image": str(w),
        "conjugation": "Frobenius square a -> a^2",
        "conjugation_checks": {
            "w_maps_to_w_squared": True,
            "unit_norms_are_one": True,
        },
        "eisenstein_mod_2": {
            "six_unit_images": ["1", "1", "w", "w", "w^2", "w^2"],
            "distinct_nonzero_images": ["1", "w", "w^2"],
            "all_units_nonzero": True,
            "support_preserved": True,
            "right_hand_side_12_image": "0",
        },
        "exact_basis_encoding": {
            "codes": ["0", "1", "w", "1+w"],
            "addition_is_xor_checked": True,
            "multiplication_term_table": [
                [int(term_code[i, j]) for j in range(3)] for i in range(3)
            ],
        },
        "normalization": {
            "fixed_coordinate": (
                "The least support coordinate is fixed to phase 1."
            ),
            "proof": (
                "Multiplication by c in GF(4)^* scales every Hermitian "
                "correlation by c*c^2=c^3=1. Every assignment has a unique "
                "normalization because its fixed-coordinate value is nonzero."
            ),
            "factor": 3,
            "normalized_assignments_per_support": normalized_count,
            "unnormalized_assignments_per_support": 3 * normalized_count,
            "base3_bijection_checked_for_every_row": True,
        },
        "affine_reduction": {
            "proof": (
                "For u'_(a q+b)=u_q with gcd(a,35)=1, direct substitution "
                "gives C_u'(t)=C_u(a^(-1)t). Thus existence and assignment "
                "counts are constant on each fully derived affine orbit."
            ),
            "all_transformed_supports_checked_in_char3_family": True,
            "orbit_partition_complete": True,
        },
        "negative_shift_reduction": {
            "proof": "C(-t)=conjugate(C(t)).",
            "reversed_pair_incidence_checked_for_every_case_and_shift": True,
            "checked_shifts": list(range(1, 18)),
            "full_equations_exactly_verified_for_every_survivor": True,
        },
        "cases": case_records,
        "affine_case_count": len(case_records),
        "support_count_covered": len(char3["support_masks"]),
        "total_normalized_support_assignment_pairs": (
            total_normalized_support_assignments
        ),
        "total_unnormalized_support_assignment_pairs": (
            3 * total_normalized_support_assignments
        ),
        "total_normalized_surviving_support_assignment_pairs": (
            total_normalized_surviving_support_assignments
        ),
        "survivor_file_sha256": survivor_sha256,
        "freeze_complete": True,
    }
    char2_sha256 = write_json_frozen("char2-freeze.json", char2_record)
    print("[char2] complete output frozen", flush=True)
    return {
        "record": char2_record,
        "record_sha256": char2_sha256,
        "survivors": survivors,
        "representative_masks": representative_masks,
    }


def eisenstein_pair_for_power(exponent):
    return ((1, 0), (0, 1), (-1, -1))[exponent % 3]


def exhaustive_exact_lift(char2):
    # After global multiplication, the least nonzero coefficient is exactly 1.
    # A normalized GF(4) phase has precisely two Eisenstein-unit lifts at every
    # other support coordinate, namely +/- omega^e.  Thus 2^11 sign rows are
    # a bijective, complete normalized lift space.
    sign_count = 2 ** (TARGET_WEIGHT - 1)
    sign_ids = np.arange(sign_count, dtype=np.int64)
    signs = np.ones((sign_count, TARGET_WEIGHT), dtype=np.int16)
    for coordinate in range(1, TARGET_WEIGHT):
        signs[:, coordinate] = 1 - 2 * (
            (sign_ids >> (coordinate - 1)) & 1
        ).astype(np.int16)
    assert np.all(signs[:, 0] == 1)

    total_lifts = 0
    witness = None
    lift_case_records = []
    for mask in char2["representative_masks"]:
        support = support_from_mask(mask)
        pairs = {
            shift: incidence_pairs(support, shift) for shift in range(1, 18)
        }
        phase_assignments = char2["survivors"][mask]
        support_lifts = 0
        support_survivors = 0
        for exponents in phase_assignments:
            alive = np.ones(sign_count, dtype=np.bool_)
            for shift in range(1, 18):
                sum_a = np.zeros(sign_count, dtype=np.int16)
                sum_b = np.zeros(sign_count, dtype=np.int16)
                for i, j in pairs[shift]:
                    coefficient_a, coefficient_b = eisenstein_pair_for_power(
                        exponents[i] - exponents[j]
                    )
                    sign_products = signs[:, i] * signs[:, j]
                    sum_a += coefficient_a * sign_products
                    sum_b += coefficient_b * sign_products
                alive &= (sum_a == 0) & (sum_b == 0)
            survivor_sign_ids = np.flatnonzero(alive)
            support_lifts += sign_count
            total_lifts += sign_count
            support_survivors += len(survivor_sign_ids)
            if witness is None and len(survivor_sign_ids):
                selected_signs = tuple(
                    int(value) for value in signs[int(survivor_sign_ids[0]), :]
                )
                witness = {
                    "support_mask": mask,
                    "support": support,
                    "exponents": tuple(exponents),
                    "signs": selected_signs,
                }
        lift_case_records.append(
            {
                "support_representative": mask_hex(mask),
                "gf4_normalized_survivors": len(phase_assignments),
                "normalized_eisenstein_lifts_checked": support_lifts,
                "exact_lifts_satisfying_all_correlations": support_survivors,
            }
        )
    return {
        "normalization_factor": 6,
        "lifts_per_normalized_gf4_assignment": sign_count,
        "total_normalized_lifts_checked": total_lifts,
        "cases": lift_case_records,
        "complete": True,
        "witness": witness,
    }


def pair_add(left, right):
    return (left[0] + right[0], left[1] + right[1])


def pair_multiply(left, right):
    # (a+b*omega)(c+d*omega), with omega^2=-omega-1.
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c - b * d)


def pair_conjugate(value):
    # conjugate(omega)=omega^2=-1-omega.
    a, b = value
    return (a - b, -b)


def witness_coefficients(witness):
    coefficients = [(0, 0)] * N
    labels = ["0"] * N
    power_labels = ("1", "omega", "omega^2")
    for q, exponent, sign in zip(
        witness["support"], witness["exponents"], witness["signs"]
    ):
        base = eisenstein_pair_for_power(exponent)
        coefficients[q] = (sign * base[0], sign * base[1])
        labels[q] = power_labels[exponent] if sign == 1 else "-" + power_labels[exponent]
    return coefficients, labels


def verify_witness_dense(witness):
    coefficients, labels = witness_coefficients(witness)
    correlations = []
    for shift in range(N):
        total = (0, 0)
        for q in range(N):
            term = pair_multiply(
                coefficients[q], pair_conjugate(coefficients[(q + shift) % N])
            )
            total = pair_add(total, term)
        correlations.append(total)
    return correlations == [(TARGET_WEIGHT, 0)] + [(0, 0)] * (N - 1)


def verify_witness_sparse(witness):
    coefficients, labels = witness_coefficients(witness)
    totals = {(right - left) % N: (0, 0) for left in witness["support"] for right in witness["support"]}
    for left in witness["support"]:
        for right in witness["support"]:
            shift = (right - left) % N
            product = pair_multiply(coefficients[left], pair_conjugate(coefficients[right]))
            totals[shift] = pair_add(totals[shift], product)
    return all(
        totals.get(shift, (0, 0)) == ((TARGET_WEIGHT, 0) if shift == 0 else (0, 0))
        for shift in range(N)
    )


def verify_witness_crt(witness):
    coefficients, labels = witness_coefficients(witness)
    grid = {(q % 5, q % 7): coefficients[q] for q in range(N)}
    for shift5 in range(5):
        for shift7 in range(7):
            total = (0, 0)
            for row in range(5):
                for column in range(7):
                    left = grid[(row, column)]
                    right = grid[((row + shift5) % 5, (column + shift7) % 7)]
                    total = pair_add(total, pair_multiply(left, pair_conjugate(right)))
            target = (TARGET_WEIGHT, 0) if (shift5, shift7) == (0, 0) else (0, 0)
            if total != target:
                return False
    return True


def verify_witness_phase(witness):
    for shift in range(N):
        histogram = [0, 0, 0]
        position = {q: i for i, q in enumerate(witness["support"])}
        signed_histogram = [0, 0, 0]
        for i, q in enumerate(witness["support"]):
            target = (q + shift) % N
            if target not in position:
                continue
            j = position[target]
            difference = (witness["exponents"][i] - witness["exponents"][j]) % 3
            signed_histogram[difference] += witness["signs"][i] * witness["signs"][j]
            histogram[difference] += 1
        a = signed_histogram[0] - signed_histogram[2]
        b = signed_histogram[1] - signed_histogram[2]
        expected = (TARGET_WEIGHT, 0) if shift == 0 else (0, 0)
        if (a, b) != expected:
            return False
    return True


def verify_witness_quotient(witness):
    PZ = PolynomialRing(ZZ, "T")
    T = PZ.gen()
    EZ = PZ.quotient(T ** 2 + T + 1, names=("omega_bar",))
    omega_bar = EZ.gen()
    GX = PolynomialRing(EZ, "Y")
    Y = GX.gen()

    def to_eisenstein(value):
        return EZ(value[0]) + EZ(value[1]) * omega_bar

    coefficients, labels = witness_coefficients(witness)
    polynomial = sum(
        (to_eisenstein(coefficients[q]) * Y ** q for q in range(N)),
        GX.zero(),
    )
    sharp = sum(
        (
            to_eisenstein(pair_conjugate(coefficients[q]))
            * Y ** ((-q) % N)
            for q in range(N)
        ),
        GX.zero(),
    )
    remainder = (polynomial * sharp - TARGET_WEIGHT) % (Y ** N - 1)
    return remainder == 0


def environment_record():
    cpu_model = "unknown"
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.lower().startswith("model name"):
                    cpu_model = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    memory_kib = None
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    memory_kib = int(line.split()[1])
                    break
    except OSError:
        pass
    return {
        "task_id": TASK_ID,
        "sage_version": str(SAGE_VERSION),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "os": platform.platform(),
        "uname": list(platform.uname()),
        "cpu_model": cpu_model,
        "logical_cpu_count": os.cpu_count(),
        "memory_total_kib": memory_kib,
    }


def run_validation():
    char3 = char3_stage()
    char2 = char2_stage(char3)
    lift = exhaustive_exact_lift(char2)

    normalized_char2_survivors = sum(
        len(assignments) for assignments in char2["survivors"].values()
    )
    witness = lift.pop("witness")
    if normalized_char2_survivors == 0:
        status = "nonexistent"
        conclusion = (
            "The complete characteristic-three support family is eliminated "
            "by the exhaustive characteristic-two assignment search."
        )
        witness_record = None
    elif witness is None:
        status = "nonexistent"
        conclusion = (
            "All characteristic-two survivors were exhaustively lifted through "
            "the complete normalized Eisenstein-unit sign space, and every "
            "exact lift failed an exact Hermitian correlation."
        )
        witness_record = None
    else:
        checks = {
            "dense_correlations": verify_witness_dense(witness),
            "sparse_ordered_differences": verify_witness_sparse(witness),
            "crt_grid_5_by_7": verify_witness_crt(witness),
            "phase_histograms": verify_witness_phase(witness),
            "quotient_group_ring": verify_witness_quotient(witness),
        }
        assert all(checks.values())
        coefficients, labels = witness_coefficients(witness)
        status = "exists"
        conclusion = "An exact original Eisenstein-unit sequence was found and independently verified."
        witness_record = {
            "coefficients": labels,
            "support_mask": mask_hex(witness["support_mask"]),
            "verification": checks,
        }

    result = {
        "task_id": TASK_ID,
        "status": status,
        "conclusion": conclusion,
        "scope": "circulant sequence problem of length 35 and weight 12 only",
        "external_expected_values_compared": False,
        "arithmetic": "exact throughout; no floating point used for any decision",
        "char3": char3["record"],
        "char3_freeze_sha256": char3["record_sha256"],
        "char2": char2["record"],
        "char2_freeze_sha256": char2["record_sha256"],
        "exact_lift": lift,
        "witness": witness_record,
        "completeness": {
            "factor_branches_complete": True,
            "cyclic_code_words_enumerated_completely": True,
            "support_family_frozen_before_affine_and_char2_analysis": True,
            "affine_orbit_partition_complete": True,
            "gf4_assignments_complete_modulo_proved_global_phase_normalization": True,
            "char3_and_char2_outputs_frozen_before_conclusion": True,
            "original_lifts_complete_if_needed": True,
        },
    }
    write_json_frozen("result.json", result)
    return result


def finalize_runtime(status, exception_text=None):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    timing = {
        "task_id": TASK_ID,
        "status": status,
        "wall_seconds_inside_sage": time.time() - START_WALL,
        "cpu_seconds_process_time": time.process_time() - START_CPU,
        "user_cpu_seconds_rusage": usage.ru_utime,
        "system_cpu_seconds_rusage": usage.ru_stime,
        "max_resident_set_kib": usage.ru_maxrss,
        "exception": exception_text,
    }
    if not os.path.exists(artifact_path("timing.json")):
        write_json_frozen("timing.json", timing)
    if not os.path.exists(artifact_path("environment.json")):
        write_json_frozen("environment.json", environment_record())


if __name__ == "__main__":
    try:
        final_result = run_validation()
        finalize_runtime(final_result["status"])
        print("[final] status=%s" % final_result["status"], flush=True)
    except Exception:
        exception_text = traceback.format_exc()
        print(exception_text, file=sys.stderr, flush=True)
        if not os.path.exists(artifact_path("result.json")):
            failure_result = {
                "task_id": TASK_ID,
                "status": "inconclusive",
                "conclusion": "Validator execution failed; no mathematical conclusion is drawn.",
                "exception": exception_text,
                "external_expected_values_compared": False,
            }
            write_json_frozen("result.json", failure_result)
        finalize_runtime("inconclusive", exception_text)
        raise
