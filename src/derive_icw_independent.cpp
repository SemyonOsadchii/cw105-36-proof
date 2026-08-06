// Standalone exhaustive cross-check for the multiplier-4 ICW_3(35,36) reduction.
//
// This file deliberately does not share implementation code with derive_icw.py.
// It scans all 7^9 orbit coefficient assignments in base 7, retains those with
// sum 6 and squared norm 36, and then checks all 35 periodic correlations.

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {
constexpr int N = 35;
constexpr int ORBIT_COUNT = 9;
constexpr std::int64_t TOTAL_ASSIGNMENTS = 40353607; // 7^9

constexpr std::array<int, ORBIT_COUNT> orbit_sizes = {1, 6, 6, 6, 3, 6, 2, 2, 3};

const std::array<std::vector<int>, ORBIT_COUNT> orbits = {{
    {0},
    {1, 4, 16, 29, 11, 9},
    {2, 8, 32, 23, 22, 18},
    {3, 12, 13, 17, 33, 27},
    {5, 20, 10},
    {6, 24, 26, 34, 31, 19},
    {7, 28},
    {14, 21},
    {15, 25, 30},
}};

using Coefficients = std::array<int, ORBIT_COUNT>;
using Row = std::array<int, N>;

Row build_row(const Coefficients& coefficients) {
    Row row{};
    for (int orbit = 0; orbit < ORBIT_COUNT; ++orbit) {
        for (const int index : orbits[orbit]) {
            row[index] = coefficients[orbit];
        }
    }
    return row;
}

bool is_icw(const Row& row) {
    for (int shift = 0; shift < N; ++shift) {
        int correlation = 0;
        for (int i = 0; i < N; ++i) {
            correlation += row[i] * row[(i + shift) % N];
        }
        const int target = (shift == 0) ? 36 : 0;
        if (correlation != target) {
            return false;
        }
    }
    return true;
}

void print_coefficients(const Coefficients& coefficients) {
    std::cout << '[';
    for (int i = 0; i < ORBIT_COUNT; ++i) {
        if (i) std::cout << ',';
        std::cout << coefficients[i];
    }
    std::cout << ']';
}

void print_row(const Row& row) {
    std::cout << '[';
    for (int i = 0; i < N; ++i) {
        if (i) std::cout << ',';
        std::cout << row[i];
    }
    std::cout << ']';
}

bool unit_three_maps(const Row& source, const Row& target) {
    Row image{};
    for (int i = 0; i < N; ++i) {
        image[(3 * i) % N] = source[i];
    }
    return image == target;
}
} // namespace

int main() {
    std::vector<Coefficients> scalar_candidates;
    std::vector<Row> icw_solutions;
    scalar_candidates.reserve(1500);

    for (std::int64_t code = 0; code < TOTAL_ASSIGNMENTS; ++code) {
        std::int64_t value = code;
        Coefficients coefficients{};
        int linear_sum = 0;
        int squared_norm = 0;

        for (int orbit = 0; orbit < ORBIT_COUNT; ++orbit) {
            const int coefficient = static_cast<int>(value % 7) - 3;
            value /= 7;
            coefficients[orbit] = coefficient;
            linear_sum += orbit_sizes[orbit] * coefficient;
            squared_norm += orbit_sizes[orbit] * coefficient * coefficient;
        }

        if (linear_sum != 6 || squared_norm != 36) {
            continue;
        }
        scalar_candidates.push_back(coefficients);
        const Row row = build_row(coefficients);
        if (is_icw(row)) {
            icw_solutions.push_back(row);
        }
    }

    bool equivalent_by_unit_three = false;
    if (icw_solutions.size() == 2) {
        equivalent_by_unit_three = unit_three_maps(icw_solutions[0], icw_solutions[1]) ||
                                   unit_three_maps(icw_solutions[1], icw_solutions[0]);
    }

    std::cout << "{\n";
    std::cout << "  \"total_assignments\": " << TOTAL_ASSIGNMENTS << ",\n";
    std::cout << "  \"scalar_candidate_count\": " << scalar_candidates.size() << ",\n";
    std::cout << "  \"scalar_candidates\": [";
    for (std::size_t i = 0; i < scalar_candidates.size(); ++i) {
        if (i) std::cout << ',';
        print_coefficients(scalar_candidates[i]);
    }
    std::cout << "],\n";
    std::cout << "  \"invariant_solution_count\": " << icw_solutions.size() << ",\n";
    std::cout << "  \"invariant_solutions\": [";
    for (std::size_t i = 0; i < icw_solutions.size(); ++i) {
        if (i) std::cout << ',';
        print_row(icw_solutions[i]);
    }
    std::cout << "],\n";
    std::cout << "  \"equivalent_by_unit_3\": "
              << (equivalent_by_unit_three ? "true" : "false") << "\n";
    std::cout << "}\n";

    return (scalar_candidates.size() == 1434 && icw_solutions.size() == 2 &&
            equivalent_by_unit_three)
               ? 0
               : 1;
}
