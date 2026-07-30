#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using i64 = long long;

static constexpr int N = 50000;

static std::vector<int> primes_upto(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int p = 2; 1LL * p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int k = p * p; k <= n; k += p) is_prime[k] = false;
    }
    std::vector<int> ps;
    for (int p = 2; p <= n; ++p) if (is_prime[p]) ps.push_back(p);
    return ps;
}

static int mod_norm(i64 x, int p) {
    x %= p;
    if (x < 0) x += p;
    return static_cast<int>(x);
}

static i64 apery_P(i64 m) {
    const i64 m2 = m * m;
    const i64 m3 = m2 * m;
    return 34 * m3 + 51 * m2 + 27 * m + 5;
}

struct Star {
    int p;
    int r;
    int j;
    char branch; // D=direct, R=reflected, C=central
};

int main() {
    const auto primes = primes_upto(N);
    std::vector<std::vector<int>> hits(N + 1);
    long long recurrence_steps = 0;
    long long zero_pairs = 0;

    for (int p : primes) {
        const int limit = std::min(p - 1, N - p);
        if (limit < 1) continue;

        std::vector<int> inv(limit + 1, 0);
        inv[1] = 1;
        for (int k = 2; k <= limit; ++k) {
            inv[k] = mod_norm(p - (1LL * (p / k) * inv[p % k]) % p, p);
            assert(1LL * k * inv[k] % p == 1);
        }

        int b_prev = 1 % p;      // b_0
        int b_cur = 5 % p;       // b_1
        if (b_cur == 0) {
            hits[p + 1].push_back(p);
            ++zero_pairs;
        }

        for (int m = 1; m < limit; ++m) {
            const i64 m3 = 1LL * m * m * m;
            const int rhs = mod_norm(
                1LL * mod_norm(apery_P(m), p) * b_cur
                - 1LL * mod_norm(m3, p) * b_prev,
                p
            );
            const int u = inv[m + 1];
            const int den_inv = static_cast<int>(1LL * u * u % p * u % p);
            const int b_next = static_cast<int>(1LL * rhs * den_inv % p);
            ++recurrence_steps;
            const int r = m + 1;
            if (b_next == 0) {
                const int n = p + r;
                assert(n <= N);
                assert(2 * p > n && p <= n);
                hits[n].push_back(p);
                ++zero_pairs;
            }
            b_prev = b_cur;
            b_cur = b_next;
        }
    }

    std::vector<int> rows3;
    int maxK = 0;
    int rows2 = 0;
    int rows_ge3 = 0;
    int rows_ge4 = 0;

    std::cout << "SCAN N=" << N
              << " primes=" << primes.size()
              << " recurrence_steps=" << recurrence_steps
              << " target_pairs=" << zero_pairs << "\n";

    for (int n = 1; n <= N; ++n) {
        auto &ps = hits[n];
        std::sort(ps.begin(), ps.end());
        ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
        const int K = static_cast<int>(ps.size());
        maxK = std::max(maxK, K);
        if (K == 2) ++rows2;
        if (K >= 3) {
            ++rows_ge3;
            rows3.push_back(n);
            if (K >= 4) ++rows_ge4;
            std::cout << "ROW n=" << n << " K=" << K << " stars=";
            bool first = true;
            for (int p : ps) {
                const int r = n - p;
                const int reflected = p - 1 - r;
                assert(r >= 0 && reflected >= 0);
                const int j = std::min(r, reflected);
                char branch = 'C';
                if (r < reflected) branch = 'D';
                else if (reflected < r) branch = 'R';
                if (!first) std::cout << ";";
                first = false;
                std::cout << "(p=" << p
                          << ",r=" << r
                          << ",j=" << j
                          << ",b=" << branch << ")";
            }
            std::cout << "\n";
        }
    }

    const std::vector<int> expected = {
        321, 11576, 18444, 22101, 26164, 47066, 47859
    };
    assert(rows3 == expected);
    assert(maxK == 3);
    assert(rows_ge3 == 7);
    assert(rows_ge4 == 0);

    std::cout << "SUMMARY maxK=" << maxK
              << " rows_K2=" << rows2
              << " rows_Kge3=" << rows_ge3
              << " rows_Kge4=" << rows_ge4 << "\n";
    std::cout << "EXPECTED_TRIPLE_LIST_VERIFIED=yes\n";
    return 0;
}
