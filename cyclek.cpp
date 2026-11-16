// Simulated-annealing search for a k-tick cycle of emitter states that
// rebuilds the pattern once per cycle.
// stdin:  L n kmin kmax iters  start[6]  then per emitter t: L lines "pos mask"
//         then T, S  (masks: 7 little-endian 64-bit words)
// stdout: "none" or "k" followed by 6 lines of k state indices (0 = off, 1..L = line)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <random>
const int W = 7;  // 448 bits: enough for R = 12 (397 cells)
struct u128 {
    uint64_t w[W];
    u128() { for (int i = 0; i < W; i++) w[i] = 0; }
    u128(int) { for (int i = 0; i < W; i++) w[i] = 0; }
    u128 operator&(const u128& o) const { u128 r; for (int i = 0; i < W; i++) r.w[i] = w[i] & o.w[i]; return r; }
    u128 operator|(const u128& o) const { u128 r; for (int i = 0; i < W; i++) r.w[i] = w[i] | o.w[i]; return r; }
    u128 operator~() const { u128 r; for (int i = 0; i < W; i++) r.w[i] = ~w[i]; return r; }
    u128& operator&=(const u128& o) { for (int i = 0; i < W; i++) w[i] &= o.w[i]; return *this; }
    u128& operator|=(const u128& o) { for (int i = 0; i < W; i++) w[i] |= o.w[i]; return *this; }
    bool any() const { for (int i = 0; i < W; i++) if (w[i]) return true; return false; }
    bool operator!() const { return !any(); }
    int count() const { int c = 0; for (int i = 0; i < W; i++) c += __builtin_popcountll(w[i]); return c; }
};

int L, n, kmin, kmax;
long long iters;
int startp[6];
int pos[6][32];
u128 bm[6][32];
u128 Tm, Sm;
int K;
int seq[6][64];
std::mt19937_64 rng(12345);

// mask as W space-separated 64-bit words, least significant first
static u128 rd() {
    u128 r;
    for (int i = 0; i < W; i++) {
        unsigned long long x;
        if (scanf("%llu", &x) != 1) exit(1);
        r.w[i] = x;
    }
    return r;
}

int ringDist(int a, int b) {
    int d = ((b - a) % n + n) % n;
    return std::min(d, n - d);
}

// an emitter's cyclic sequence is executable: consecutive lines differ by <= 1,
// and during a run of r off ticks it can move r-1 steps (toggle off, move, toggle on)
bool valid(int t) {
    int firstOn = -1;
    for (int i = 0; i < K; i++) if (seq[t][i]) { firstOn = i; break; }
    if (firstOn < 0) return true;
    int i = firstOn;
    for (int c = 0; c < K; c++) {
        int a = seq[t][i];
        int j = (i + 1) % K;
        if (seq[t][j]) {
            if (abs(seq[t][j] - a) > 1) return false;
            i = j;
        } else {
            int r = 0;
            while (!seq[t][j]) { r++; j = (j + 1) % K; }
            if (ringDist(pos[t][a], pos[t][seq[t][j]]) > r - 1) return false;
            c += r;
            i = j;
        }
    }
    return true;
}

int evalErr() {
    u128 M;
    int err = 0;
    for (int pass = 0; pass < 2; pass++) {
        bool matched = false;
        for (int i = 0; i < K; i++) {
            u128 one, two, three;
            for (int t = 0; t < 6; t++) {
                u128 b = bm[t][seq[t][i]];
                three |= two & b; two |= one & b; one |= b;
            }
            M = (M | (two & ~three)) & ~three;
            if (!(Tm & ~M) && !(Sm & M)) { matched = true; M &= ~Tm; }
        }
        if (!matched) {
            err += (Tm & ~M).count() + (Sm & M).count() + 1;
        }
    }
    return err;
}

bool anneal() {
    for (int t = 0; t < 6; t++) for (int i = 0; i < K; i++) seq[t][i] = 0;
    int cur = evalErr();
    double T0 = 2.0, T1 = 0.05;
    for (long long it = 0; it < iters; it++) {
        if (cur == 0) return true;
        double temp = T0 * pow(T1 / T0, (double)it / iters);
        int t = rng() % 6, i = rng() % K;
        int old = seq[t][i];
        int nv;
        int r = rng() % 4;
        if (r == 0) nv = 0;
        else if (r == 1 && old) nv = std::max(1, std::min(L, old + ((rng() & 1) ? 1 : -1)));
        else nv = rng() % (L + 1);
        if (nv == old) continue;
        seq[t][i] = nv;
        if (!valid(t)) { seq[t][i] = old; continue; }
        int e = evalErr();
        if (e <= cur || std::uniform_real_distribution<double>(0, 1)(rng) < exp((cur - e) / temp)) cur = e;
        else seq[t][i] = old;
    }
    return cur == 0;
}

int main() {
    if (scanf("%d %d %d %d %lld", &L, &n, &kmin, &kmax, &iters) != 5) return 1;
    for (int t = 0; t < 6; t++) if (scanf("%d", &startp[t]) != 1) return 1;
    for (int t = 0; t < 6; t++) {
        bm[t][0] = u128(); pos[t][0] = startp[t];
        for (int i = 1; i <= L; i++) {
            if (scanf("%d", &pos[t][i]) != 1) return 1;
            bm[t][i] = rd();
        }
    }
    Tm = rd(); Sm = rd();
    for (K = kmin; K <= kmax; K++) {
        for (int restart = 0; restart < 3; restart++) {
            if (anneal()) {
                printf("%d\n", K);
                for (int t = 0; t < 6; t++) {
                    for (int i = 0; i < K; i++) printf("%d ", seq[t][i]);
                    printf("\n");
                }
                return 0;
            }
        }
    }
    printf("none\n");
    return 0;
}
