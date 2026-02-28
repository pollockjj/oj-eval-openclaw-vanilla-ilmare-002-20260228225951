#include "int2048.h"

namespace sjtu {

// ======================== Helpers ========================

void int2048::trim() {
    while (d.size() > 1 && d.back() == 0) d.pop_back();
    if (d.size() == 1 && d[0] == 0) neg = false;
}

bool int2048::isZero() const {
    return d.size() == 1 && d[0] == 0;
}

int int2048::cmpAbs(const int2048 &o) const {
    if (d.size() != o.d.size()) return d.size() < o.d.size() ? -1 : 1;
    for (int i = (int)d.size() - 1; i >= 0; i--) {
        if (d[i] != o.d[i]) return d[i] < o.d[i] ? -1 : 1;
    }
    return 0;
}

void int2048::addAbs(const int2048 &o) {
    long long carry = 0;
    size_t n = std::max(d.size(), o.d.size());
    d.resize(n);
    for (size_t i = 0; i < n; i++) {
        long long sum = d[i] + (i < o.d.size() ? o.d[i] : 0) + carry;
        carry = sum / BASE;
        d[i] = sum % BASE;
    }
    if (carry) d.push_back(carry);
}

void int2048::subAbs(const int2048 &o) {
    long long borrow = 0;
    for (size_t i = 0; i < d.size(); i++) {
        long long diff = d[i] - (i < o.d.size() ? o.d[i] : 0) - borrow;
        if (diff < 0) { diff += BASE; borrow = 1; }
        else borrow = 0;
        d[i] = diff;
    }
    trim();
}

// ======================== NTT ========================

static long long power(long long a, long long b, long long mod) {
    long long res = 1; a %= mod;
    while (b > 0) {
        if (b & 1) res = (__int128)res * a % mod;
        a = (__int128)a * a % mod;
        b >>= 1;
    }
    return res;
}

static const long long MOD1 = 998244353, G1 = 3;

static void ntt(std::vector<long long> &a, bool inv) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        long long w = inv ? power(G1, MOD1 - 1 - (MOD1 - 1) / len, MOD1) : power(G1, (MOD1 - 1) / len, MOD1);
        for (int i = 0; i < n; i += len) {
            long long wn = 1;
            for (int j = 0; j < len / 2; j++) {
                long long u = a[i + j], v = (__int128)a[i + j + len / 2] * wn % MOD1;
                a[i + j] = (u + v) % MOD1;
                a[i + j + len / 2] = (u - v + MOD1) % MOD1;
                wn = (__int128)wn * w % MOD1;
            }
        }
    }
    if (inv) {
        long long inv_n = power(n, MOD1 - 2, MOD1);
        for (auto &x : a) x = (__int128)x * inv_n % MOD1;
    }
}

// NTT multiplication using 3 primes + CRT
static const long long MOD2 = 985661441, G2 = 3;  // 2^23 * 117 + 1
static const long long MOD3 = 754974721, G3 = 11; // 2^24 * 45 + 1

static void ntt2(std::vector<long long> &a, bool inv, long long mod, long long g) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        long long w = inv ? power(g, mod - 1 - (mod - 1) / len, mod) : power(g, (mod - 1) / len, mod);
        for (int i = 0; i < n; i += len) {
            long long wn = 1;
            for (int j = 0; j < len / 2; j++) {
                long long u = a[i + j], v = (__int128)a[i + j + len / 2] * wn % mod;
                a[i + j] = (u + v) % mod;
                a[i + j + len / 2] = (u - v + mod) % mod;
                wn = (__int128)wn * w % mod;
            }
        }
    }
    if (inv) {
        long long inv_n = power(n, mod - 2, mod);
        for (auto &x : a) x = (__int128)x * inv_n % mod;
    }
}

std::vector<long long> int2048::multiply_vec(const std::vector<long long> &a, const std::vector<long long> &b) {
    if (a.empty() || b.empty()) return {0};
    static const long long NB = 10000;
    int sa = 2 * a.size(), sb = 2 * b.size();
    int n = 1;
    while (n < sa + sb) n <<= 1;
    
    // Split digits into base-NB
    std::vector<long long> fa(n, 0), fb(n, 0);
    for (size_t i = 0; i < a.size(); i++) {
        fa[2*i] = a[i] % NB;
        fa[2*i+1] = a[i] / NB;
    }
    for (size_t i = 0; i < b.size(); i++) {
        fb[2*i] = b[i] % NB;
        fb[2*i+1] = b[i] / NB;
    }
    
    // Max convolution value: NB^2 * n = 10^8 * n
    // For n up to ~10M, need 3 primes.
    // M1*M2*M3 ~ 7.4*10^26, enough for up to ~7.4*10^18 per convolution value.
    // But each individual prime is ~10^9, so conv values up to ~10^9 are OK per mod.
    // The CRT reconstruction is needed when the true value > any single prime.
    
    auto fa2 = fa, fb2 = fb;
    
    // NTT with mod1
    ntt(fa, false); ntt(fb, false);
    for (int i = 0; i < n; i++) fa[i] = (__int128)fa[i] * fb[i] % MOD1;
    ntt(fa, true);
    
    // NTT with mod2
    ntt2(fa2, false, MOD2, G2); ntt2(fb2, false, MOD2, G2);
    for (int i = 0; i < n; i++) fa2[i] = (__int128)fa2[i] * fb2[i] % MOD2;
    ntt2(fa2, true, MOD2, G2);
    
    // CRT with 2 primes: reconstruct x ≡ r1 (mod M1), x ≡ r2 (mod M2)
    // x = r1 + M1 * ((r2 - r1) * M1^{-1} mod M2)
    long long inv_m1_m2 = power(MOD1, MOD2 - 2, MOD2);
    
    std::vector<long long> conv(n);
    long long carry = 0;
    for (int i = 0; i < n; i++) {
        long long r1 = fa[i], r2 = fa2[i];
        long long t = (__int128)(r2 - r1 + MOD2) % MOD2 * inv_m1_m2 % MOD2;
        __int128 x = (__int128)r1 + (__int128)t * MOD1 + carry;
        conv[i] = (long long)(x % NB);
        carry = (long long)(x / NB);
    }
    
    // Convert base-NB to base-BASE
    int rlen = (n + 1) / 2 + 2;
    std::vector<long long> result(rlen, 0);
    for (int i = 0; i < n; i++) {
        int idx = i / 2;
        if (i % 2 == 0) result[idx] += conv[i];
        else result[idx] += conv[i] * NB;
    }
    carry = 0;
    for (int i = 0; i < rlen; i++) {
        result[i] += carry;
        carry = result[i] / BASE;
        result[i] %= BASE;
    }
    while (carry) { result.push_back(carry % BASE); carry /= BASE; }
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

// ======================== Constructors ========================

int2048::int2048() : d(1, 0), neg(false) {}

int2048::int2048(long long v) : neg(false) {
    if (v < 0) {
        neg = true;
        unsigned long long uv = (unsigned long long)(-(v + 1)) + 1ULL;
        while (uv > 0) { d.push_back((long long)(uv % BASE)); uv /= BASE; }
    } else if (v == 0) {
        d.push_back(0);
    } else {
        while (v > 0) { d.push_back(v % BASE); v /= BASE; }
    }
    if (d.empty()) d.push_back(0);
}

int2048::int2048(const std::string &s) : neg(false) { read(s); }
int2048::int2048(const int2048 &o) : d(o.d), neg(o.neg) {}

// ======================== read / print ========================

void int2048::read(const std::string &s) {
    d.clear(); neg = false;
    int start = 0;
    if (!s.empty() && s[0] == '-') { neg = true; start = 1; }
    else if (!s.empty() && s[0] == '+') { start = 1; }
    while (start < (int)s.size() - 1 && s[start] == '0') start++;
    for (int i = (int)s.size() - 1; i >= start; i -= BASEDIGITS) {
        long long val = 0;
        int beg = std::max(start, i - BASEDIGITS + 1);
        for (int j = beg; j <= i; j++) val = val * 10 + (s[j] - '0');
        d.push_back(val);
    }
    if (d.empty()) d.push_back(0);
    trim();
}

void int2048::print() {
    if (neg) putchar('-');
    printf("%lld", d.back());
    for (int i = (int)d.size() - 2; i >= 0; i--) printf("%08lld", d[i]);
}

// ======================== add / minus (Integer 1) ========================

int2048 &int2048::add(const int2048 &o) {
    if (neg == o.neg) { addAbs(o); }
    else {
        int c = cmpAbs(o);
        if (c >= 0) { subAbs(o); }
        else { int2048 tmp(o); tmp.subAbs(*this); *this = tmp; }
    }
    return *this;
}

int2048 add(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::minus(const int2048 &o) {
    int2048 tmp(o);
    tmp.neg = !tmp.neg;
    if (tmp.isZero()) tmp.neg = false;
    return this->add(tmp);
}

int2048 minus(int2048 a, const int2048 &b) { return a.minus(b); }

// ======================== Unary ========================

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const {
    int2048 res(*this);
    if (!res.isZero()) res.neg = !res.neg;
    return res;
}

// ======================== Assignment ========================

int2048 &int2048::operator=(const int2048 &o) {
    if (this != &o) { d = o.d; neg = o.neg; }
    return *this;
}

// ======================== Arithmetic operators ========================

int2048 &int2048::operator+=(const int2048 &o) { return this->add(o); }
int2048 operator+(int2048 a, const int2048 &b) { return a += b; }

int2048 &int2048::operator-=(const int2048 &o) { return this->minus(o); }
int2048 operator-(int2048 a, const int2048 &b) { return a -= b; }

int2048 &int2048::operator*=(const int2048 &o) {
    if (isZero() || o.isZero()) { d.assign(1, 0); neg = false; return *this; }
    bool rn = neg != o.neg;
    if (d.size() <= 64 || o.d.size() <= 64) {
        size_t total = d.size() + o.d.size();
        std::vector<long long> res(total, 0);
        for (size_t i = 0; i < d.size(); i++) {
            long long carry = 0;
            for (size_t j = 0; j < o.d.size(); j++) {
                __int128 cur = (__int128)d[i] * o.d[j] + res[i+j] + carry;
                res[i+j] = (long long)(cur % BASE);
                carry = (long long)(cur / BASE);
            }
            if (carry) res[i + o.d.size()] += carry;
        }
        d = res;
    } else {
        d = multiply_vec(d, o.d);
    }
    neg = rn; trim(); return *this;
}

int2048 operator*(int2048 a, const int2048 &b) { return a *= b; }

// ======================== Division ========================

// Multiply big int by single cell
static void mulSingle(const std::vector<long long> &a, long long s, std::vector<long long> &res) {
    res.resize(a.size() + 1, 0);
    long long carry = 0;
    for (size_t i = 0; i < a.size(); i++) {
        __int128 val = (__int128)a[i] * s + carry;
        res[i] = (long long)(val % int2048::BASE);
        carry = (long long)(val / int2048::BASE);
    }
    res[a.size()] = carry;
    while (res.size() > 1 && res.back() == 0) res.pop_back();
}

static int cmpVec(const std::vector<long long> &a, const std::vector<long long> &b) {
    if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
    for (int i = (int)a.size() - 1; i >= 0; i--)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

static void subVecInPlace(std::vector<long long> &a, const std::vector<long long> &b) {
    long long borrow = 0;
    for (size_t i = 0; i < a.size(); i++) {
        long long diff = a[i] - (i < b.size() ? b[i] : 0) - borrow;
        if (diff < 0) { diff += int2048::BASE; borrow = 1; } else borrow = 0;
        a[i] = diff;
    }
    while (a.size() > 1 && a.back() == 0) a.pop_back();
}

// Schoolbook long division: a / b, both positive, b.d.size() >= 2
static void schoolbookDiv(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r) {
    int n = (int)a.d.size(), m = (int)b.d.size();
    q.d.assign(n - m + 1, 0); q.neg = false;
    std::vector<long long> rem;
    
    for (int i = n - 1; i >= 0; i--) {
        rem.insert(rem.begin(), a.d[i]);
        while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
        if (i > n - m) continue;
        if (cmpVec(rem, b.d) < 0) continue;
        
        long long lo = 1, hi = int2048::BASE - 1;
        while (lo < hi) {
            long long mid = lo + (hi - lo + 1) / 2;
            std::vector<long long> tmp;
            mulSingle(b.d, mid, tmp);
            if (cmpVec(tmp, rem) <= 0) lo = mid; else hi = mid - 1;
        }
        
        q.d[i] = lo;
        std::vector<long long> prod;
        mulSingle(b.d, lo, prod);
        subVecInPlace(rem, prod);
    }
    q.trim();
    r.d = rem; r.neg = false; r.trim();
}

// Newton reciprocal: compute floor(BASE^(2*k) / B) where B has m cells, k >= m
// Returns a number with about k - m + 1 cells (the inverse scaled by BASE^k)
// 
// Actually, we compute inv = floor(BASE^(m + k) / B) where k is the desired quotient precision.
// This inv has about k + 1 cells.
// Then q = floor(A * inv / BASE^(m + k))
//
// Newton iteration at a FIXED scale S = m + k:
// We want X = floor(BASE^S / B)
// Newton: X' = X * (2*BASE^S - B*X) / BASE^S
//       = 2*X - B*X^2 / BASE^S
//
// With precision doubling: at precision p cells,
// X_p has p cells, and we only compute with p cells of B.
// 
// More precisely, at precision p:
// b_p = floor(B / BASE^(m - p))  (top p cells of B, has p cells)
// scale_p = p + p = 2p  (we work at scale 2p locally)
// X_p ≈ BASE^(2p) / b_p  (X_p has about p cells)
// Newton: X'_{2p} = X_p * (2*BASE^(2p) - b_{2p} * X_p) / BASE^(2p)
//   This X'_{2p} ≈ BASE^(2*2p) / b_{2p} ??? No, that changes the scale.
//
// OK I think the cleanest formulation is:
// 
// At each precision p, we have X_p such that X_p ≈ BASE^(2p) / b_p
// where b_p = top p cells of B (= floor(B / BASE^(m-p)) if m >= p, else B).
//
// To double precision to 2p:
// 1. b_{2p} = top 2p cells of B
// 2. X_p_scaled = X_p * BASE^(2p)  (scale up to match new scale)
//    Now X_p_scaled ≈ BASE^(4p) / b_p
//    But we want ≈ BASE^(4p) / b_{2p}
//    Since b_{2p} ≈ b_p * BASE^p + b_low_p, and b_low_p < BASE^p:
//    BASE^(4p) / b_{2p} ≈ BASE^(4p) / (b_p * BASE^p) = BASE^(3p) / b_p
//    So X_p_scaled is off by roughly factor BASE^p / BASE^(2p) ... 
//
// This is getting complicated. Let me use a completely different approach.

// Compute X ≈ floor(BASE^(2*prec) / b_prec)
// where b_prec = top 'prec' cells of B (floor(B / BASE^(m - prec)))
// B has m cells total. Returns X with about prec + 1 cells.
static int2048 computeReciprocal(const int2048 &B, int m, int prec) {
    // Get top prec cells of B
    int2048 bp;
    if (prec >= m) {
        bp = B;
    } else {
        bp.d.assign(B.d.end() - prec, B.d.end());
    }
    bp.neg = false; bp.trim();
    
    if (prec <= 2) {
        // Base case: compute floor(BASE^(2*prec) / bp) using schoolbook
        // bp has at most 2 cells. BASE^(2*prec) has 2*prec+1 cells.
        // For prec=1: BASE^2 / bp[0]
        // For prec=2: BASE^4 / (bp[0] + bp[1]*BASE)
        int2048 num;
        num.d.assign(2 * prec + 1, 0);
        num.d[2 * prec] = 1;
        num.neg = false;
        
        // Simple long division since bp is small
        int2048 qq, rr;
        int2048::divmod(num, bp, qq, rr);
        return qq;
    }
    
    int half = (prec + 1) / 2;
    
    // Recurse at half precision
    int2048 x = computeReciprocal(B, m, half);
    // x ≈ BASE^(2*half) / b_half where b_half = top half cells of B
    
    // Scale up: x_scaled = x * BASE^(prec - half)
    // This approximates BASE^(prec + half) / b_half ≈ BASE^(2*prec) / b_prec
    int pad = prec - half;
    int2048 xs;
    xs.d.assign(pad, 0);
    xs.d.insert(xs.d.end(), x.d.begin(), x.d.end());
    xs.neg = false; xs.trim();
    
    // Newton refinement:
    // X' = 2*xs - bp * xs^2 / BASE^(2*prec)
    // To avoid computing xs^2 then bp*xs^2 (3 muls), use:
    // X' = xs * (2*BASE^(2*prec) - bp * xs) / BASE^(2*prec)
    // = xs * delta / BASE^(2*prec)
    // where delta = 2*BASE^(2*prec) - bp * xs
    
    // Compute bp * xs
    int2048 bx = bp;
    bx *= xs;
    
    // delta = 2*BASE^(2*prec) - bx
    int S = 2 * prec;
    int2048 twoBaseS;
    twoBaseS.d.assign(S + 1, 0);
    twoBaseS.d[S] = 2;
    twoBaseS.neg = false;
    
    int2048 delta = twoBaseS;
    delta -= bx;
    
    // X' = xs * delta / BASE^(2*prec)
    int2048 xd = xs;
    xd *= delta;
    
    // Shift right by 2*prec
    int2048 result;
    if ((int)xd.d.size() <= S) {
        result = int2048(0LL);
    } else {
        result.d.assign(xd.d.begin() + S, xd.d.end());
        result.neg = xd.neg;
        result.trim();
    }
    
    return result;
}

static void newtonDiv(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r) {
    int n = (int)a.d.size(), m = (int)b.d.size();
    int2048 abs_b(b); abs_b.neg = false;
    int2048 abs_a(a); abs_a.neg = false;
    
    if (n <= 2 * m + 2) {
        // Balanced case: use Newton reciprocal
        // Compute X ≈ BASE^(2m) / B, has ~m+1 cells
        int2048 X = computeReciprocal(abs_b, m, m);
        
        int S = 2 * m;
        int2048 AX = abs_a;
        AX *= X;
        
        if ((int)AX.d.size() <= S) {
            q = int2048(0LL);
        } else {
            q.d.assign(AX.d.begin() + S, AX.d.end());
            q.neg = false; q.trim();
        }
        
        r = abs_a;
        int2048 qb = q;
        qb *= abs_b;
        r -= qb;
        
        // At most ~2 adjustments
        while (r.neg && !r.isZero()) { q -= int2048(1LL); r += abs_b; }
        while (!r.neg && r.cmpAbs(abs_b) >= 0) { q += int2048(1LL); r -= abs_b; }
    } else {
        // Unbalanced: chunk-based long division with m-cell chunks
        q = int2048(0LL);
        int2048 rem(0LL);
        
        int nchunks = (n + m - 1) / m;
        for (int i = nchunks - 1; i >= 0; i--) {
            int lo = i * m;
            int hi = std::min((i + 1) * m, n);
            int chunk_size = hi - lo;
            
            // rem = rem * BASE^chunk_size + a[lo..hi)
            if (!rem.isZero()) {
                std::vector<long long> new_d(chunk_size, 0);
                new_d.insert(new_d.end(), rem.d.begin(), rem.d.end());
                rem.d = new_d;
            } else {
                rem.d.assign(chunk_size, 0);
            }
            for (int j = lo; j < hi; j++) {
                rem.d[j - lo] = abs_a.d[j];
            }
            rem.neg = false;
            rem.trim();
            
            // Divide rem by b (rem has at most 2m cells, balanced)
            int2048 q_digit, new_rem;
            int2048::divmod(rem, abs_b, q_digit, new_rem);
            rem = new_rem;
            
            // q = q * BASE^chunk_size + q_digit
            if (!q.isZero()) {
                std::vector<long long> new_d(chunk_size, 0);
                new_d.insert(new_d.end(), q.d.begin(), q.d.end());
                q.d = new_d;
            } else {
                q.d.assign(chunk_size, 0);
            }
            for (size_t j = 0; j < q_digit.d.size(); j++) {
                if (j < q.d.size()) q.d[j] += q_digit.d[j];
                else q.d.push_back(q_digit.d[j]);
            }
            long long carry = 0;
            for (size_t j = 0; j < q.d.size(); j++) {
                q.d[j] += carry;
                carry = q.d[j] / int2048::BASE;
                q.d[j] %= int2048::BASE;
            }
            if (carry) q.d.push_back(carry);
            q.trim();
        }
        r = rem;
    }
    
    q.neg = false; r.neg = false;
    q.trim(); r.trim();
}


void int2048::divmod(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r) {
    int c = a.cmpAbs(b);
    if (c < 0) { q = int2048(0LL); r = a; r.neg = false; return; }
    if (c == 0) { q = int2048(1LL); r = int2048(0LL); return; }
    
    if (b.d.size() == 1) {
        long long divisor = b.d[0];
        q.d.resize(a.d.size()); q.neg = false;
        long long rem = 0;
        for (int i = (int)a.d.size() - 1; i >= 0; i--) {
            __int128 cur = (__int128)rem * BASE + a.d[i];
            q.d[i] = (long long)(cur / divisor);
            rem = (long long)(cur % divisor);
        }
        q.trim();
        r = int2048(0LL); r.d = {rem}; r.trim();
        return;
    }
    
    int n = (int)a.d.size(), m = (int)b.d.size();
    long long complexity = (long long)(n - m + 1) * m;
    
    if (complexity <= 500000LL) {
        schoolbookDiv(a, b, q, r);
    } else {
        newtonDiv(a, b, q, r);
    }
}

int2048 &int2048::operator/=(const int2048 &o) {
    bool rn = neg != o.neg;
    int2048 q, r;
    divmod(*this, o, q, r);
    if (!r.isZero() && rn) q.add(int2048(1LL));
    *this = q;
    neg = rn;
    if (isZero()) neg = false;
    return *this;
}

int2048 operator/(int2048 a, const int2048 &b) { return a /= b; }

int2048 &int2048::operator%=(const int2048 &o) {
    int2048 q = *this / o;
    *this -= q * o;
    return *this;
}

int2048 operator%(int2048 a, const int2048 &b) { return a %= b; }

// ======================== Stream operators ========================

std::istream &operator>>(std::istream &is, int2048 &v) {
    std::string s; is >> s; v.read(s); return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &v) {
    if (v.neg) os << '-';
    os << v.d.back();
    char buf[16];
    for (int i = (int)v.d.size() - 2; i >= 0; i--) {
        snprintf(buf, sizeof(buf), "%08lld", v.d[i]);
        os << buf;
    }
    return os;
}

// ======================== Comparison operators ========================

bool operator==(const int2048 &a, const int2048 &b) { return a.neg == b.neg && a.d == b.d; }
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) {
    if (a.neg != b.neg) return a.neg;
    return a.neg ? a.cmpAbs(b) > 0 : a.cmpAbs(b) < 0;
}
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
