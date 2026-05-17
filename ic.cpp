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

uint64_t computeBound(uint64_t n) {
    double c = 3.38;

    double ln_n = log((double)n);
    double ln_ln_n = log(ln_n);

    double B = c * exp(0.5 * sqrt(ln_n * ln_ln_n));

    return (uint64_t)B + 1;
}

vector<uint64_t> buildFactorBase(uint64_t B) {
    vector<uint64_t> base;

    vector<bool> isPrime(B + 1, true);

    for (uint64_t i = 2; i <= B; i++) {
        if (isPrime[i]) {
            base.push_back(i);
            for (uint64_t j = i * 2; j <= B; j += i)
                isPrime[j] = false;
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

bool tryRelation(uint64_t alpha, uint64_t k, uint64_t p,
    const vector<uint64_t>& base,
    Relation& r) {

    uint64_t val = pow_mod(alpha, k, p);

    vector<int> exp;
    if (!factorOverBase(val, base, exp))
        return false;

    r.k = k;
    r.exponents = exp;
    return true;
}

vector<Relation> collectRelations(uint64_t alpha, uint64_t p,
    const vector<uint64_t>& base,
    size_t need) {

    vector<Relation> rel;
    uint64_t n = p - 1;

    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<uint64_t> dist(0, n - 1);

    while (rel.size() < need) {
        uint64_t k = dist(gen);
        Relation r;

        if (tryRelation(alpha, k, p, base, r))
            rel.push_back(r);
    }
    return rel;
}

vector<uint64_t> factor(uint64_t n) {
    vector<uint64_t> f;

    for (uint64_t d = 2; d * d <= n; ++d) {
        if (n % d == 0) {
            f.push_back(d);
            while (n % d == 0) n /= d;
        }
    }

    if (n > 1) f.push_back(n);

    return f;
}

bool is_generator(uint64_t a, uint64_t p) {
    uint64_t phi = p - 1;
    auto fac = factor(phi);

    for (auto q : fac) {
        if (pow_mod(a, phi / q, p) == 1)
            return false;
    }

    return true;
}

uint64_t findGenerator(uint64_t p) {
    for (uint64_t g = 2; g < p; ++g) {
        if (is_generator(g, p))
            return g;
    }

    return 0;
}

uint64_t solveCongruence(uint64_t u, uint64_t v, uint64_t mod) {
    uint64_t g = gcd(u, mod);

    u /= g;
    v /= g;
    mod /= g;

    uint64_t inv = inv_mod(u, mod);

    return (v * inv) % mod;
}

vector<uint64_t> SLSMod(vector<vector<int>> A, vector<uint64_t> b, uint64_t mod) {
    int n = A.size();
    int m = A[0].size();

    vector<vector<int>> mat(n, vector<int>(m));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            mat[i][j] = ((A[i][j] % (int)mod) + mod) % mod;
        }
    }

    vector<int> where(m, -1);

    int row = 0;

    for (int col = 0; col < m && row < n; col++) {
        int sel = row;
        while (sel < n && mat[sel][col] == 0)
            sel++;

        if (sel == n)
            continue;

        swap(mat[sel], mat[row]);
        swap(b[sel], b[row]);

        where[col] = row;

        int64_t inv = inv_mod(mat[row][col], mod);
        if (inv == 0) continue;

        for (int j = col; j < m; j++) {
            mat[row][j] = (int)((1LL * mat[row][j] * inv) % mod);
        }
        b[row] = (b[row] * inv) % mod;

        for (int i = 0; i < n; i++) {
            if (i != row && mat[i][col] != 0) {
                int64_t factor = mat[i][col];

                for (int j = col; j < m; j++) {
                    mat[i][j] = (mat[i][j] - factor * mat[row][j]) % (int)mod;
                    if (mat[i][j] < 0) mat[i][j] += mod;
                }

                b[i] = (b[i] - factor * b[row]) % (int)mod;
                if ((int64_t)b[i] < 0) b[i] += mod;
            }
        }

        row++;
    }

    vector<uint64_t> x(m, 0);

    for (int i = 0; i < m; i++) {
        if (where[i] != -1) {
            x[i] = b[where[i]] % mod;
        }
    }

    return x;
}

uint64_t computeLogBeta(uint64_t alpha, uint64_t beta, uint64_t p, const vector<uint64_t>& base, const vector<uint64_t>& log_p) {
    uint64_t mod = p - 1;

    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<uint64_t> dist(0, mod - 1);

    vector<int> exp;

    while (true) {

        uint64_t l = dist(gen);

        uint64_t value = (beta * pow_mod(alpha, l, p)) % p;

        if (factorOverBase(value, base, exp)) {

            uint64_t sum = 0;

            for (size_t i = 0; i < base.size(); i++) {

                sum = (sum + (uint64_t)exp[i] * log_p[i]) % mod;
            }

            uint64_t x =
                (sum + mod - l % mod) % mod;

            if (pow_mod(alpha, x, p) == beta) {
                return x;
            }
        }
    }
}

uint64_t index_calculus( uint64_t alpha, uint64_t beta, uint64_t p){
    uint64_t mod = p - 1;

    uint64_t B = computeBound(p);
    vector<uint64_t> base = buildFactorBase(B);

    size_t need = base.size() + 10;
    vector<Relation> relations = collectRelations(alpha, p, base, need);

    vector<vector<int>> A;
    vector<uint64_t> b;

    A.reserve(relations.size());
    b.reserve(relations.size());

    for (auto& r : relations) {
        A.push_back(r.exponents);
        b.push_back(r.k % mod);
    }

    vector<uint64_t> log_p = SLSMod(A, b, mod);
    uint64_t x = computeLogBeta(alpha, beta, p, base, log_p);
    return x;

}

uint64_t index_calculus_general(uint64_t alpha, uint64_t beta, uint64_t p) {
    uint64_t mod = p - 1;

    if (is_generator(alpha, p)) {
        return index_calculus(alpha, beta, p);
    }

    uint64_t gamma = findGenerator(p);

    uint64_t u = index_calculus(gamma, alpha, p);
    uint64_t v = index_calculus(gamma, beta, p);

    return solveCongruence(u, v, mod);
}

vector<Relation> collectRelationsParallel(uint64_t alpha, uint64_t p, const vector<uint64_t>& base, size_t n) {
    vector<Relation> relations;

    mutex mtx;
    set<uint64_t> used_k;
    uint64_t mod = p - 1;

    unsigned threads_count = thread::hardware_concurrency();

    if (threads_count == 0) {
        threads_count = 4;
    }

    vector<thread> threads;

    auto worker = [&]() {
            random_device rd;
            mt19937_64 gen(rd());

            uniform_int_distribution<uint64_t>
                dist(0, mod - 1);

            while (true) {
                uint64_t k = dist(gen);

                Relation r;

                if (!tryRelation(alpha, k, p, base, r))
                {
                    continue;
                }

                lock_guard<mutex> lock(mtx);

                if (relations.size() >= n) {
                    return;
                }

                if (used_k.count(k)) {
                    continue;
                }

                used_k.insert(k);

                relations.push_back(r);
            }
        };

    for (unsigned i = 0; i < threads_count; i++) {
        threads.emplace_back(worker);
    }

    for (auto& t : threads) {
        t.join();
    }

    return relations;
}

uint64_t index_calculus_parallel(uint64_t alpha, uint64_t beta, uint64_t p){
    uint64_t mod = p - 1;
    uint64_t B = computeBound(p);
    vector<uint64_t> base = buildFactorBase(B);

    size_t need = base.size() + 10;

    vector<Relation> relations =
        collectRelationsParallel(alpha, p, base, need);

    vector<vector<int>> A;
    vector<uint64_t> b;

    A.reserve(relations.size());
    b.reserve(relations.size());

    for (auto& r : relations) {
        A.push_back(r.exponents);
        b.push_back(r.k % mod);
    }

    vector<uint64_t> log_p = SLSMod(A, b, mod);
    uint64_t x = computeLogBeta( alpha, beta, p, base,log_p);
    return x;
}

uint64_t index_calculus_general_parallel( uint64_t alpha, uint64_t beta, uint64_t p) {
    uint64_t mod = p - 1;

    if (is_generator(alpha, p)) {
        return index_calculus_parallel(alpha, beta, p);
    }

    uint64_t gamma = findGenerator(p);
    uint64_t u = index_calculus_parallel(gamma, alpha, p);
    uint64_t v = index_calculus_parallel(gamma, beta, p);
    return solveCongruence(u, v, mod);
}