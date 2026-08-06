// Independent exact checker for the generalized length-35, weight-12 mu_6 theorem.
//
// It imposes no prescribed zero coordinates.  The program:
//   1. verifies the relevant F3 factorization and irreducibility;
//   2. enumerates both 3^12 cyclic codes and all weight-12 words;
//   3. proves that the resulting 420 supports form one affine orbit;
//   4. directly checks the four-shift F4 contradiction on one representative.
//
// No floating point, SAT solver, timeout, or heuristic is used.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Poly = std::vector<int>;  // coefficients low to high
using Mask = std::uint64_t;

static int f3_add(int a, int b) { return (a + b) % 3; }
static int f3_mul(int a, int b) { return (a * b) % 3; }
static int f3_inv(int a) {
    a %= 3;
    if (a == 1) return 1;
    if (a == 2) return 2;
    throw std::runtime_error("F3 inverse of zero");
}
static int f3_neg(int a) { return (3 - (a % 3)) % 3; }

static int f4_add(int a, int b) { return a ^ b; }
static int f4_mul(int a, int b) {
    const int a0 = a & 1;
    const int a1 = (a >> 1) & 1;
    const int b0 = b & 1;
    const int b1 = (b >> 1) & 1;
    const int c0 = (a0 & b0) ^ (a1 & b1);
    const int c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1);
    return c0 | (c1 << 1);
}
static int f4_inv(int a) {
    if (a == 1) return 1;
    if (a == 2) return 3;
    if (a == 3) return 2;
    throw std::runtime_error("F4 inverse of zero");
}
static int f4_ratio(int a, int b) { return f4_mul(a, f4_inv(b)); }

static Poly trim(Poly p) {
    while (p.size() > 1 && p.back() == 0) p.pop_back();
    if (p.empty()) p.push_back(0);
    return p;
}

static bool poly_equal(const Poly& left, const Poly& right) {
    return trim(left) == trim(right);
}

static Poly poly_add(const Poly& left, const Poly& right) {
    Poly result(std::max(left.size(), right.size()), 0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        const int a = i < left.size() ? left[i] : 0;
        const int b = i < right.size() ? right[i] : 0;
        result[i] = f3_add(a, b);
    }
    return trim(result);
}

static Poly poly_sub(const Poly& left, const Poly& right) {
    Poly negative = right;
    for (int& value : negative) value = f3_neg(value);
    return poly_add(left, negative);
}

static Poly poly_mul(const Poly& left, const Poly& right) {
    Poly result(left.size() + right.size() - 1, 0);
    for (std::size_t i = 0; i < left.size(); ++i) {
        for (std::size_t j = 0; j < right.size(); ++j) {
            result[i + j] = f3_add(result[i + j], f3_mul(left[i], right[j]));
        }
    }
    return trim(result);
}

static std::pair<Poly, Poly> poly_divmod(Poly dividend, const Poly& divisor_raw) {
    Poly divisor = trim(divisor_raw);
    dividend = trim(dividend);
    if (divisor.size() == 1 && divisor[0] == 0) {
        throw std::runtime_error("polynomial division by zero");
    }
    Poly quotient(dividend.size() >= divisor.size() ? dividend.size() - divisor.size() + 1 : 1, 0);
    while (!(dividend.size() == 1 && dividend[0] == 0) && dividend.size() >= divisor.size()) {
        const std::size_t degree = dividend.size() - divisor.size();
        const int coefficient = f3_mul(dividend.back(), f3_inv(divisor.back()));
        quotient[degree] = f3_add(quotient[degree], coefficient);
        for (std::size_t i = 0; i < divisor.size(); ++i) {
            dividend[degree + i] = f3_add(
                dividend[degree + i], f3_neg(f3_mul(coefficient, divisor[i]))
            );
        }
        dividend = trim(dividend);
    }
    return {trim(quotient), trim(dividend)};
}

static Poly poly_mod(const Poly& polynomial, const Poly& modulus) {
    return poly_divmod(polynomial, modulus).second;
}

static Poly poly_gcd(Poly left, Poly right) {
    while (!(right.size() == 1 && right[0] == 0)) {
        Poly remainder = poly_mod(left, right);
        left = right;
        right = remainder;
    }
    if (left.size() == 1 && left[0] == 0) return left;
    const int scale = f3_inv(left.back());
    for (int& value : left) value = f3_mul(scale, value);
    return trim(left);
}

static Poly poly_pow_mod(Poly base, std::uint64_t exponent, const Poly& modulus) {
    Poly result{1};
    base = poly_mod(base, modulus);
    while (exponent) {
        if (exponent & 1U) result = poly_mod(poly_mul(result, base), modulus);
        base = poly_mod(poly_mul(base, base), modulus);
        exponent >>= 1U;
    }
    return result;
}

static std::uint64_t integer_power(std::uint64_t base, int exponent) {
    std::uint64_t result = 1;
    for (int i = 0; i < exponent; ++i) result *= base;
    return result;
}

static std::vector<int> prime_divisors(int value) {
    std::vector<int> result;
    for (int p = 2; p * p <= value; ++p) {
        if (value % p == 0) {
            result.push_back(p);
            while (value % p == 0) value /= p;
        }
    }
    if (value > 1) result.push_back(value);
    return result;
}

static bool irreducible(const Poly& polynomial_raw) {
    const Poly polynomial = trim(polynomial_raw);
    const int degree = static_cast<int>(polynomial.size()) - 1;
    if (degree <= 0 || polynomial.back() != 1) return false;
    const Poly x{0, 1};
    const Poly x_mod = poly_mod(x, polynomial);
    if (!poly_equal(poly_pow_mod(x, integer_power(3, degree), polynomial), x_mod)) return false;
    for (int prime : prime_divisors(degree)) {
        const Poly power = poly_pow_mod(x, integer_power(3, degree / prime), polynomial);
        if (poly_gcd(polynomial, poly_sub(power, x_mod)).size() > 1) return false;
    }
    return true;
}

static std::vector<int> codeword(const std::array<int,12>& message, const Poly& generator) {
    std::vector<int> result(35, 0);
    for (int shift = 0; shift < 12; ++shift) {
        if (message[shift] == 0) continue;
        for (std::size_t j = 0; j < generator.size(); ++j) {
            result[shift + static_cast<int>(j)] = f3_add(
                result[shift + static_cast<int>(j)],
                f3_mul(message[shift], generator[j])
            );
        }
    }
    return result;
}

static bool zero_f3_autocorrelation(const std::vector<int>& row) {
    for (int shift = 0; shift < 35; ++shift) {
        int correlation = 0;
        for (int q = 0; q < 35; ++q) {
            correlation = f3_add(correlation, f3_mul(row[q], row[(q + shift) % 35]));
        }
        if (correlation != 0) return false;
    }
    return true;
}

static Mask mask_of(const std::vector<int>& row) {
    Mask mask = 0;
    for (int q = 0; q < 35; ++q) {
        if (row[q] != 0) mask |= Mask{1} << q;
    }
    return mask;
}

static Mask mask_of_positions(const std::array<int,12>& positions) {
    Mask mask = 0;
    for (int q : positions) mask |= Mask{1} << q;
    return mask;
}

static Mask affine_transform(Mask mask, int multiplier, int translation) {
    Mask result = 0;
    for (int q = 0; q < 35; ++q) {
        if (mask & (Mask{1} << q)) {
            result |= Mask{1} << ((multiplier * q + translation) % 35);
        }
    }
    return result;
}

static std::vector<int> units_mod_35() {
    std::vector<int> units;
    for (int a = 0; a < 35; ++a) {
        if (std::gcd(a, 35) == 1) units.push_back(a);
    }
    return units;
}

static int f4_correlation(
    const std::array<int,35>& values,
    const std::array<int,12>& support,
    int shift
) {
    int result = 0;
    for (int q : support) {
        const int target = (q + shift) % 35;
        if (values[target] != 0) result = f4_add(result, f4_ratio(values[q], values[target]));
    }
    return result;
}

static void print_pair_vector(const std::vector<std::pair<int,int>>& values) {
    std::cout << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '[' << values[i].first << ',' << values[i].second << ']';
    }
    std::cout << ']';
}

int main() {
    try {
        const std::vector<Poly> fixed_factors{
            {2,1},
            {1,1,1,1,1},
            {1,1,1,1,1,1,1},
        };
        const std::array<Poly,2> reciprocal_pair{
            Poly{1,2,2,1,2,1,0,1,2,0,1,0,1},
            Poly{1,0,1,0,2,1,0,1,2,1,2,2,1},
        };
        std::vector<Poly> factors = fixed_factors;
        factors.push_back(reciprocal_pair[0]);
        factors.push_back(reciprocal_pair[1]);
        Poly product{1};
        for (const Poly& factor : factors) product = poly_mul(product, factor);
        Poly target(36, 0);
        target[0] = 2;
        target[35] = 1;
        bool factorization_valid = poly_equal(product, target);
        bool irreducibility_valid = std::all_of(
            factors.begin(), factors.end(), [](const Poly& factor) { return irreducible(factor); }
        );
        Poly reversed = reciprocal_pair[0];
        std::reverse(reversed.begin(), reversed.end());
        bool reciprocal_valid = poly_equal(reversed, reciprocal_pair[1]);
        if (!(factorization_valid && irreducibility_valid && reciprocal_valid)) {
            throw std::runtime_error("F3 factor verification failed");
        }

        Poly mandatory{1};
        for (const Poly& factor : fixed_factors) mandatory = poly_mul(mandatory, factor);
        const std::array<Poly,2> generators{
            poly_mul(mandatory, reciprocal_pair[0]),
            poly_mul(mandatory, reciprocal_pair[1]),
        };
        if (generators[0].size() != 24 || generators[1].size() != 24) {
            throw std::runtime_error("unexpected F3 generator degree");
        }

        std::array<int,2> word_counts{0,0};
        std::array<std::set<Mask>,2> support_sets;
        const std::uint64_t message_count = integer_power(3, 12);
        for (int code = 0; code < 2; ++code) {
            for (std::uint64_t encoded = 0; encoded < message_count; ++encoded) {
                std::uint64_t digits = encoded;
                std::array<int,12> message{};
                for (int i = 0; i < 12; ++i) {
                    message[i] = static_cast<int>(digits % 3U);
                    digits /= 3U;
                }
                const std::vector<int> word = codeword(message, generators[code]);
                const int weight = static_cast<int>(std::count_if(
                    word.begin(), word.end(), [](int value) { return value != 0; }
                ));
                if (weight != 12) continue;
                if (!zero_f3_autocorrelation(word)) {
                    throw std::runtime_error("weight-12 F3 codeword failed autocorrelation");
                }
                ++word_counts[code];
                support_sets[code].insert(mask_of(word));
            }
        }

        std::set<Mask> intersection;
        std::set_intersection(
            support_sets[0].begin(), support_sets[0].end(),
            support_sets[1].begin(), support_sets[1].end(),
            std::inserter(intersection, intersection.begin())
        );
        std::set<Mask> supports = support_sets[0];
        supports.insert(support_sets[1].begin(), support_sets[1].end());

        const std::array<int,12> representative{
            0,1,2,3,7,10,12,16,21,22,26,28
        };
        const Mask representative_mask = mask_of_positions(representative);
        const std::vector<int> units = units_mod_35();
        std::set<Mask> orbit;
        std::vector<std::pair<int,int>> stabilizer;
        for (int a : units) {
            for (int b = 0; b < 35; ++b) {
                const Mask transformed = affine_transform(representative_mask, a, b);
                orbit.insert(transformed);
                if (transformed == representative_mask) stabilizer.push_back({a,b});
            }
        }
        if (supports != orbit) throw std::runtime_error("F3 supports are not the claimed affine orbit");
        if (!supports.contains(representative_mask)) {
            throw std::runtime_error("representative support is absent from the F3 orbit");
        }
        for (Mask mask : orbit) {
            std::array<int,12> positions{};
            int index = 0;
            for (int q = 0; q < 35; ++q) {
                if (mask & (Mask{1} << q)) positions[index++] = q;
            }
            if (positions < representative) {
                throw std::runtime_error("stored representative is not lexicographically affine-canonical");
            }
        }

        const std::array<int,3> nonzero_values{1,2,3};
        const std::uint64_t assignments = integer_power(3, 11);
        std::uint64_t c8_zero = 0;
        std::uint64_t c2_c8_zero = 0;
        std::uint64_t c2_c8_c10_zero = 0;
        std::uint64_t all_four_zero = 0;
        std::uint64_t identity_failures = 0;
        std::array<int,35> values{};
        values.fill(0);
        values[representative[0]] = 1;
        std::array<int,12> exponents{};
        exponents.fill(0);

        for (std::uint64_t encoded = 0; encoded < assignments; ++encoded) {
            std::uint64_t digits = encoded;
            for (int i = 1; i < 12; ++i) {
                exponents[i] = static_cast<int>(digits % 3U);
                digits /= 3U;
                values[representative[i]] = nonzero_values[exponents[i]];
            }
            const int c2 = f4_correlation(values, representative, 2);
            const int c4 = f4_correlation(values, representative, 4);
            const int c8 = f4_correlation(values, representative, 8);
            const int c10 = f4_correlation(values, representative, 10);
            if (c8 == 0) {
                ++c8_zero;
                const int p = f4_ratio(values[2], values[10]);
                if (p != f4_ratio(values[28], values[1])) {
                    ++identity_failures;
                    continue;
                }
                int four_terms = 0;
                for (const auto& pair : std::array<std::pair<int,int>,4>{{
                    {0,10}, {28,3}, {2,12}, {26,1}
                }}) {
                    four_terms = f4_add(four_terms, f4_ratio(values[pair.first], values[pair.second]));
                }
                if (f4_mul(p, c2) != four_terms) {
                    ++identity_failures;
                    continue;
                }
                const int residual = f4_add(
                    f4_ratio(values[12], values[22]),
                    f4_ratio(values[16], values[26])
                );
                if (c10 != f4_add(four_terms, residual)) {
                    ++identity_failures;
                    continue;
                }
                if (c2 == 0) {
                    ++c2_c8_zero;
                    if (c10 == 0) {
                        ++c2_c8_c10_zero;
                        if (residual != 0) {
                            ++identity_failures;
                            continue;
                        }
                        if (f4_ratio(values[12], values[16]) != f4_ratio(values[22], values[26])) {
                            ++identity_failures;
                            continue;
                        }
                        const int forced_c4 = f4_ratio(values[3], values[7]);
                        if (forced_c4 == 0 || c4 != forced_c4) {
                            ++identity_failures;
                            continue;
                        }
                    }
                }
            }
            if (c2 == 0 && c4 == 0 && c8 == 0 && c10 == 0) ++all_four_zero;
        }

        if (word_counts != std::array<int,2>{420,420}) throw std::runtime_error("unexpected F3 word counts");
        if (support_sets[0].size() != 210 || support_sets[1].size() != 210) {
            throw std::runtime_error("unexpected F3 support counts");
        }
        if (!intersection.empty() || supports.size() != 420 || orbit.size() != 420) {
            throw std::runtime_error("unexpected support union or orbit size");
        }
        if (stabilizer != std::vector<std::pair<int,int>>{{1,0},{29,28}}) {
            throw std::runtime_error("unexpected affine stabilizer");
        }
        if (identity_failures != 0 || all_four_zero != 0) {
            throw std::runtime_error("F4 four-shift contradiction failed");
        }

        std::cout << '{';
        std::cout << "\"statement\":\"No 12-sparse perfect sixth-root sequence of length 35 exists.\",";
        std::cout << "\"uses_prescribed_zero_positions\":false,";
        std::cout << "\"f3_factorization_valid\":true,";
        std::cout << "\"mod3_candidates_per_code\":[" << word_counts[0] << ',' << word_counts[1] << "],";
        std::cout << "\"mod3_supports_per_code\":[" << support_sets[0].size() << ',' << support_sets[1].size() << "],";
        std::cout << "\"mod3_support_intersection_size\":" << intersection.size() << ',';
        std::cout << "\"mod3_support_count\":" << supports.size() << ',';
        std::cout << "\"affine_orbit_size\":" << orbit.size() << ',';
        std::cout << "\"affine_stabilizer\":";
        print_pair_vector(stabilizer);
        std::cout << ',';
        std::cout << "\"canonical_support\":[";
        for (std::size_t i = 0; i < representative.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << representative[i];
        }
        std::cout << "],";
        std::cout << "\"core_shifts\":[2,4,8,10],";
        std::cout << "\"normalized_f4_assignments_checked\":" << assignments << ',';
        std::cout << "\"core_counts\":{";
        std::cout << "\"c8_zero\":" << c8_zero << ',';
        std::cout << "\"c2_c8_zero\":" << c2_c8_zero << ',';
        std::cout << "\"c2_c8_c10_zero\":" << c2_c8_c10_zero << ',';
        std::cout << "\"all_four_zero\":" << all_four_zero << ',';
        std::cout << "\"algebraic_identity_failures\":" << identity_failures;
        std::cout << "},";
        std::cout << "\"support_masks\":[";
        bool first = true;
        for (Mask mask : supports) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << mask;
        }
        std::cout << "]}" << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "mu6 checker error: " << error.what() << '\n';
        return 1;
    }
}
