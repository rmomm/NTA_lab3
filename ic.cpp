#include "ic.h"

uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result = 1;
    a %= mod;

    while (b > 0) {
        if (b & 1) {
            result = (result * a) % mod;
        }

        a = (a * a) % mod;

        b >>= 1;
    }

    return result;
}

uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t temp = b;
        b = a % b;
        a = temp;
    }

    return a;
}

int64_t ext_gcd(int64_t a, int64_t b, int64_t& x,int64_t& y) {
    if (b == 0) {
        x = 1;
        y = 0;

        return a;
    }

    int64_t x1, y1;
    int64_t g = ext_gcd(b, a % b, x1, y1);

    x = y1;
    y = x1 - (a / b) * y1;

    return g;
}

uint64_t inv_mod(uint64_t a, uint64_t mod) {
    int64_t x, y;

    int64_t g =ext_gcd(a, mod, x, y);

    if (g != 1) {
        return 0;
    }

    x %= mod;

    if (x < 0) {
        x += mod;
    }

    return x;
}

vector<uint64_t> buildFactorBase(uint64_t bound) {
    vector<uint64_t> base;

    for (uint64_t p = 2; p <= bound; ++p) {
        bool prime = true;

        for (uint64_t d = 2; d * d <= p; ++d) {
            if (p % d == 0) {
                prime = false;
                break;
            }
        }

        if (prime) {
            base.push_back(p);
        }
    }

    return base;
}

bool factorOverBase(uint64_t v, const vector<uint64_t>& base, vector<int>& exp) {
    exp.assign(base.size(), 0);

    for (size_t i = 0; i < base.size(); ++i) {
        while (v % base[i] == 0) {
            v /= base[i];
            exp[i]++;
        }
    }

    return (v == 1);
}