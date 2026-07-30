#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using u64 = std::uint64_t;
using i64 = std::int64_t;

static constexpr int NMAX = 200000;

struct Star {
    int p;
    int r;
    int j;
    char branch; // D, R, or C
};

static inline int addmod(int a, int b, int p) {
    int s = a + b;
    if (s >= p) s -= p;
    return s;
}

static inline int submod(int a, int b, int p) {
    int s = a - b;
    if (s < 0) s += p;
    return s;
}

static inline int mulmod(i64 a, i64 b, int p) {
    return int((u64(a) * u64(b)) % u64(p));
}

static int powmod(int a, int e, int p) {
    int out = 1;
    while (e > 0) {
        if (e & 1) out = mulmod(out, a, p);
        a = mulmod(a, a, p);
        e >>= 1;
    }
    return out;
}

static inline int Pmod(int m, int p) {
    int x = m % p;
    int x2 = mulmod(x, x, p);
    int x3 = mulmod(x2, x, p);
    int v = 0;
    v = (v + mulmod(34 % p, x3, p)) % p;
    v = (v + mulmod(51 % p, x2, p)) % p;
    v = (v + mulmod(27 % p, x, p)) % p;
    v = (v + 5) % p;
    return v;
}

static inline int cube_mod(int m, int p) {
    int x = m % p;
    return mulmod(mulmod(x, x, p), x, p);
}

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> prime(n + 1, true);
    prime[0] = prime[1] = false;
    for (int p = 2; i64(p) * p <= n; ++p) {
        if (!prime[p]) continue;
        for (int k = p * p; k <= n; k += p) prime[k] = false;
    }
    std::vector<int> out;
    for (int p = 7; p <= n; ++p) if (prime[p]) out.push_back(p);
    return out;
}

// Compute b_0,...,b_limit modulo p.  Requires limit < p.
static std::vector<int> apery_b_mod(int p, int limit) {
    assert(0 <= limit && limit < p);
    std::vector<int> b(limit + 1, 0), inv(limit + 1, 0);
    b[0] = 1 % p;
    if (limit == 0) return b;
    b[1] = 5 % p;
    inv[1] = 1;
    for (int k = 2; k <= limit; ++k) {
        inv[k] = int((p - (u64(p / k) * u64(inv[p % k])) % u64(p)) % p);
    }
    for (int m = 1; m < limit; ++m) {
        int rhs = submod(mulmod(Pmod(m, p), b[m], p),
                         mulmod(cube_mod(m, p), b[m - 1], p), p);
        int dinv = powmod(inv[m + 1], 3, p);
        b[m + 1] = mulmod(rhs, dinv, p);
    }
    return b;
}

// Compute both Apéry and companion solutions modulo p through limit < p.
static std::pair<std::vector<int>, std::vector<int>>
apery_pair_mod(int p, int limit) {
    assert(0 <= limit && limit < p);
    std::vector<int> b(limit + 1, 0), a(limit + 1, 0), inv(limit + 1, 0);
    b[0] = 1 % p;
    a[0] = 0;
    if (limit == 0) return {a, b};
    b[1] = 5 % p;
    a[1] = 6 % p;
    inv[1] = 1;
    for (int k = 2; k <= limit; ++k) {
        inv[k] = int((p - (u64(p / k) * u64(inv[p % k])) % u64(p)) % p);
    }
    for (int m = 1; m < limit; ++m) {
        int dinv = powmod(inv[m + 1], 3, p);
        int bm = submod(mulmod(Pmod(m, p), b[m], p),
                        mulmod(cube_mod(m, p), b[m - 1], p), p);
        int am = submod(mulmod(Pmod(m, p), a[m], p),
                        mulmod(cube_mod(m, p), a[m - 1], p), p);
        b[m + 1] = mulmod(bm, dinv, p);
        a[m + 1] = mulmod(am, dinv, p);
    }
    return {a, b};
}

static bool contains_zero(const std::vector<std::vector<int>>& zeros, int p, int j) {
    const auto& z = zeros[p];
    return std::binary_search(z.begin(), z.end(), j);
}

static std::string star_string(const Star& s) {
    std::ostringstream os;
    os << "(" << s.p << ",r=" << s.r << ",j=" << s.j << "," << s.branch << ")";
    return os.str();
}

struct BranchStats {
    int count = 0;
    int pspan = 0;
    int jspan = 0;
    int min4_pspan = -1;
    int min4_jspan = -1;
};

static BranchStats branch_stats(const std::vector<Star>& stars, char branch) {
    std::vector<Star> v;
    for (const auto& s : stars) if (s.branch == branch) v.push_back(s);
    BranchStats out;
    out.count = int(v.size());
    if (v.empty()) return out;
    std::sort(v.begin(), v.end(), [](const Star& x, const Star& y) { return x.p < y.p; });
    out.pspan = v.back().p - v.front().p;
    int jmin = v.front().j, jmax = v.front().j;
    for (const auto& s : v) { jmin = std::min(jmin, s.j); jmax = std::max(jmax, s.j); }
    out.jspan = jmax - jmin;
    if (v.size() >= 4) {
        out.min4_pspan = std::numeric_limits<int>::max();
        out.min4_jspan = std::numeric_limits<int>::max();
        for (std::size_t i = 0; i + 3 < v.size(); ++i) {
            int ps = v[i + 3].p - v[i].p;
            int lo = v[i].j, hi = v[i].j;
            for (std::size_t k = i; k <= i + 3; ++k) {
                lo = std::min(lo, v[k].j); hi = std::max(hi, v[k].j);
            }
            out.min4_pspan = std::min(out.min4_pspan, ps);
            out.min4_jspan = std::min(out.min4_jspan, hi - lo);
        }
    }
    return out;
}

int main() {
#ifdef _OPENMP
    std::cerr << "threads=" << omp_get_max_threads() << "\n";
#endif
    const auto primes = sieve_primes(NMAX);
    std::vector<std::vector<int>> zeros(NMAX + 1);
    std::vector<std::pair<int, Star>> all_records;

#pragma omp parallel
    {
        std::vector<std::pair<int, Star>> local;
#pragma omp for schedule(dynamic, 8)
        for (int idx = 0; idx < int(primes.size()); ++idx) {
            int p = primes[idx];
            int half = (p - 1) / 2;
            int direct_max = std::min(half, NMAX - p);
            bool reflected_possible = (3LL * p <= 2LL * NMAX + 1);
            int jmax = reflected_possible ? half : direct_max;
            if (jmax < 1) continue;
            auto b = apery_b_mod(p, jmax);
            auto& zp = zeros[p];
            for (int j = 1; j <= jmax; ++j) {
                if (b[j] != 0) continue;
                zp.push_back(j);
                int nd = p + j;
                int nr = 2 * p - 1 - j;
                bool central = (p == 2 * j + 1);
                if (nd <= NMAX) {
                    Star s{p, j, j, central ? 'C' : 'D'};
                    local.push_back({nd, s});
                }
                if (!central && nr <= NMAX) {
                    Star s{p, p - 1 - j, j, 'R'};
                    local.push_back({nr, s});
                }
            }
        }
#pragma omp critical
        all_records.insert(all_records.end(), local.begin(), local.end());
    }

    std::vector<std::vector<Star>> hits(NMAX + 1);
    for (auto& rec : all_records) hits[rec.first].push_back(rec.second);
    for (auto& row : hits) {
        std::sort(row.begin(), row.end(), [](const Star& x, const Star& y) {
            return x.p < y.p;
        });
        for (std::size_t i = 1; i < row.size(); ++i) assert(row[i - 1].p != row[i].p);
    }

    // Exhaustive reflection audit for small primes.
    for (int p : primes) {
        if (p > 2000) break;
        auto b = apery_b_mod(p, p - 1);
        for (int j = 0; j <= (p - 1) / 2; ++j) {
            assert(b[j] == b[p - 1 - j]);
        }
    }

    std::vector<int> rows_ge3;
    int maxK = 0;
    int rows_ge4 = 0;
    for (int n = 1; n <= NMAX; ++n) {
        int K = int(hits[n].size());
        maxK = std::max(maxK, K);
        if (K >= 4) ++rows_ge4;
        if (K >= 3) rows_ge3.push_back(n);
    }

    const std::vector<int> expected_50000 = {
        321, 11576, 18444, 22101, 26164, 47066, 47859
    };
    std::vector<int> actual_50000;
    for (int n : rows_ge3) if (n <= 50000) actual_50000.push_back(n);
    assert(actual_50000 == expected_50000);

    std::cout << "maximum_n=" << NMAX << "\n";
    std::cout << "prime_count=" << primes.size() << "\n";
    std::cout << "zero_pairs_folded=";
    std::size_t zcount = 0;
    for (const auto& z : zeros) zcount += z.size();
    std::cout << zcount << "\n";
    std::cout << "max_K=" << maxK << "\n";
    std::cout << "rows_K_ge_3=" << rows_ge3.size() << "\n";
    std::cout << "rows_K_ge_4=" << rows_ge4 << "\n";

    for (int n : rows_ge3) {
        auto stars = hits[n];
        std::sort(stars.begin(), stars.end(), [](const Star& x, const Star& y) {
            return x.p < y.p;
        });
        int pspan = stars.back().p - stars.front().p;
        int jmin = stars.front().j, jmax = stars.front().j;
        for (const auto& s : stars) { jmin = std::min(jmin, s.j); jmax = std::max(jmax, s.j); }

        std::vector<std::string> cross;
        std::vector<std::string> repeated;
        for (std::size_t i = 0; i < stars.size(); ++i) {
            for (std::size_t k = i + 1; k < stars.size(); ++k) {
                if (stars[i].j == stars[k].j) {
                    std::ostringstream os;
                    os << stars[i].p << "=" << stars[k].p << "@j" << stars[i].j;
                    repeated.push_back(os.str());
                }
            }
        }

        // Audit actual direct/reflected targets and every cross-hit with the companion determinant.
        for (std::size_t i = 0; i < stars.size(); ++i) {
            int p = stars[i].p;
            int limit = stars[i].r;
            for (const auto& s : stars) limit = std::max(limit, s.j);
            assert(limit < p);
            auto [a, b] = apery_pair_mod(p, limit);
            assert(b[stars[i].r] == 0); // direct evaluation, no reflection shortcut
            assert(b[stars[i].j] == 0);
            int ji = stars[i].j;
            assert(ji >= 1);
            int lhs = submod(mulmod(a[ji], b[ji - 1], p),
                             mulmod(a[ji - 1], b[ji], p), p);
            int rhs = mulmod(6 % p, powmod(cube_mod(ji, p), p - 2, p), p);
            assert(lhs == rhs);
            assert(a[ji] != 0);
            for (std::size_t k = 0; k < stars.size(); ++k) {
                if (i == k) continue;
                int jk = stars[k].j;
                bool hit = (b[jk] == 0);
                bool stored = contains_zero(zeros, p, jk);
                assert(hit == stored);
                int det = submod(mulmod(a[ji], b[jk], p),
                                 mulmod(b[ji], a[jk], p), p);
                assert((det == 0) == hit);
                if (hit) {
                    std::ostringstream os;
                    os << p << "->" << stars[k].p << "(j" << ji << "->j" << jk << ")";
                    cross.push_back(os.str());
                }
            }
        }

        std::string pattern;
        for (const auto& s : stars) pattern.push_back(s.branch);
        BranchStats ds = branch_stats(stars, 'D');
        BranchStats rs = branch_stats(stars, 'R');
        BranchStats cs = branch_stats(stars, 'C');

        std::cout << "ROW n=" << n << " K=" << stars.size() << " pattern=" << pattern
                  << " pspan=" << pspan << " jspan=" << (jmax - jmin) << "\n";
        std::cout << "  stars:";
        for (const auto& s : stars) std::cout << " " << star_string(s);
        std::cout << "\n";
        auto print_stats = [](char name, const BranchStats& x) {
            std::cout << "  branch_" << name << ":count=" << x.count
                      << ",pspan=" << x.pspan << ",jspan=" << x.jspan;
            if (x.min4_pspan >= 0) {
                std::cout << ",min4_pspan=" << x.min4_pspan
                          << ",min4_jspan=" << x.min4_jspan;
            } else {
                std::cout << ",min4=none";
            }
            std::cout << "\n";
        };
        print_stats('D', ds);
        print_stats('R', rs);
        print_stats('C', cs);
        std::cout << "  cross_hits:";
        if (cross.empty()) std::cout << " none";
        else for (const auto& x : cross) std::cout << " " << x;
        std::cout << "\n";
        std::cout << "  repeated_nodes:";
        if (repeated.empty()) std::cout << " none";
        else for (const auto& x : repeated) std::cout << " " << x;
        std::cout << "\n";
    }

    std::cout << "AUDIT_OK\n";
    return 0;
}
