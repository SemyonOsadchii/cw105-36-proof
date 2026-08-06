// Direct F4 phase enumeration on all 420 characteristic-three support candidates.
//
// This is a deliberately blunt cross-check of the shorter affine-orbit/four-shift
// proof.  It reads the independently generated support certificate, assigns every
// nonzero F4 value to each active position (modulo one common scalar), and checks
// all independent shifts 1..17 directly.

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Mask = std::uint64_t;

static int f4_mul(int a, int b) {
    const int a0 = a & 1;
    const int a1 = (a >> 1) & 1;
    const int b0 = b & 1;
    const int b1 = (b >> 1) & 1;
    const int c0 = (a0 & b0) ^ (a1 & b1);
    const int c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1);
    return c0 | (c1 << 1);
}
static int f4_conjugate(int value) { return f4_mul(value, value); }

static std::vector<Mask> read_masks(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open support-mask file: " + path);
    std::vector<Mask> masks;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::stringstream stream(line);
        Mask mask = 0;
        stream >> mask;
        if (!stream || !stream.eof()) throw std::runtime_error("invalid support mask: " + line);
        if (__builtin_popcountll(mask) != 12 || (mask >> 35) != 0) {
            throw std::runtime_error("support mask is not a 12-subset of Z_35");
        }
        masks.push_back(mask);
    }
    std::sort(masks.begin(), masks.end());
    if (std::adjacent_find(masks.begin(), masks.end()) != masks.end()) {
        throw std::runtime_error("duplicate support mask");
    }
    if (masks.size() != 420) throw std::runtime_error("expected exactly 420 F3 supports");
    return masks;
}

struct SupportData {
    std::array<int,12> positions{};
    std::array<std::vector<std::pair<int,int>>,17> pairs;
};

static SupportData prepare_support(Mask mask) {
    SupportData data;
    int count = 0;
    std::array<int,35> local_index{};
    local_index.fill(-1);
    for (int q = 0; q < 35; ++q) {
        if (!(mask & (Mask{1} << q))) continue;
        data.positions[count] = q;
        local_index[q] = count;
        ++count;
    }
    if (count != 12) throw std::runtime_error("internal support size error");
    for (int shift = 1; shift <= 17; ++shift) {
        for (int i = 0; i < 12; ++i) {
            const int target = (data.positions[i] + shift) % 35;
            const int j = local_index[target];
            if (j >= 0) data.pairs[shift - 1].push_back({i,j});
        }
    }
    return data;
}

static std::uint64_t integer_power(std::uint64_t base, int exponent) {
    std::uint64_t result = 1;
    for (int i = 0; i < exponent; ++i) result *= base;
    return result;
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("usage: mu6_35_12_direct_all SUPPORT_MASKS");
        const std::vector<Mask> masks = read_masks(argv[1]);

        const std::array<int,3> values{1,3,2};
        int term[3][3]{};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                term[a][b] = f4_mul(values[a], f4_conjugate(values[b]));
            }
        }

        const std::uint64_t assignments_per_support = integer_power(3, 11);
        std::uint64_t assignments_checked = 0;
        std::uint64_t perfect_assignments = 0;
        std::uint64_t compatible_supports = 0;
        std::array<std::uint64_t,18> first_failure_counts{};

        for (Mask mask : masks) {
            const SupportData data = prepare_support(mask);
            std::array<int,12> exponent{};
            exponent.fill(0);
            std::uint64_t perfect_on_support = 0;
            for (std::uint64_t encoded = 0; encoded < assignments_per_support; ++encoded) {
                if (encoded != 0) {
                    int digit = 1;
                    while (digit < 12) {
                        ++exponent[digit];
                        if (exponent[digit] < 3) break;
                        exponent[digit] = 0;
                        ++digit;
                    }
                }
                ++assignments_checked;
                bool perfect = true;
                for (int shift_index = 0; shift_index < 17; ++shift_index) {
                    int correlation = 0;
                    for (const auto& [i,j] : data.pairs[shift_index]) {
                        correlation ^= term[exponent[i]][exponent[j]];
                    }
                    if (correlation != 0) {
                        ++first_failure_counts[shift_index + 1];
                        perfect = false;
                        break;
                    }
                }
                if (perfect) {
                    ++perfect_assignments;
                    ++perfect_on_support;
                }
            }
            if (perfect_on_support) ++compatible_supports;
        }

        const std::uint64_t expected = masks.size() * assignments_per_support;
        if (assignments_checked != expected) throw std::runtime_error("assignment count mismatch");

        std::cout << '{';
        std::cout << "\"support_count\":" << masks.size() << ',';
        std::cout << "\"uses_prescribed_zero_positions\":false,";
        std::cout << "\"global_phase_fixed\":true,";
        std::cout << "\"assignments_per_support\":" << assignments_per_support << ',';
        std::cout << "\"assignments_checked\":" << assignments_checked << ',';
        std::cout << "\"compatible_support_count\":" << compatible_supports << ',';
        std::cout << "\"zero_correlation_assignments_fixed_phase\":" << perfect_assignments << ',';
        std::cout << "\"first_failure_shift_counts\":[";
        for (int shift = 1; shift <= 17; ++shift) {
            if (shift > 1) std::cout << ',';
            std::cout << first_failure_counts[shift];
        }
        std::cout << "]}" << '\n';
        return perfect_assignments == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "direct-all checker error: " << error.what() << '\n';
        return 1;
    }
}
