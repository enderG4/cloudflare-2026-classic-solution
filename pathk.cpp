// Simulated-annealing search for a straight (non-repeating) program of K ticks,
// starting from the initial emitter state, that produces the pattern `count`
// times. Per tick each emitter is either off or on a line; validity follows
// the movement rules (one ring step per tick, toggling costs a tick).
// stdin:  L n count kmax iters seed  start[6]
//         per emitter t: L lines "pos mask" ; then T, S (masks: 7 64-bit words)
// stdout: "none" or "K" followed by 6 lines of K state indices (0 = off, 1..L = line)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <random>

const int W = 7;
const int MAXK = 400;
struct Mask { uint64_t w[W]; };

int L, n, cnt, kmax;
long long iters;
int startp[6];
int pos[6][32];
Mask bm[6][32];
Mask Tm, Sm;
int K;
int seq[6][MAXK];
int best[6][MAXK];
int bigPenalty;
std::mt19937_64 rng;

static Mask rd() {
    Mask r;
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

// emitter starts off at startp; r leading off ticks allow r moves before the
// first toggle-on; between two on-states r off ticks allow r-1 moves
bool valid(int t) {
    int prev = -1;  // index of last on tick
    for (int i = 0; i < K; i++) {
        int s = seq[t][i];
        if (!s) continue;
        if (prev < 0) {
            if (ringDist(startp[t], pos[t][s]) > i) return false;
        } else if (prev == i - 1) {
            if (abs(s - seq[t][prev]) > 1) return false;
        } else {
            if (ringDist(pos[t][seq[t][prev]], pos[t][s]) > i - prev - 2) return false;
        }
        prev = i;
    }
    return true;
}

int evalErr() {
    uint64_t M[W] = {0};
    int matches = 0;
    for (int i = 0; i < K; i++) {
        uint64_t one[W] = {0}, two[W] = {0}, three[W] = {0};
        for (int t = 0; t < 6; t++) {
            int s = seq[t][i];
            if (!s) continue;
            const uint64_t* b = bm[t][s].w;
            for (int k = 0; k < W; k++) { three[k] |= two[k] & b[k]; two[k] |= one[k] & b[k]; one[k] |= b[k]; }
        }
        bool complete = true, clean = true;
        for (int k = 0; k < W; k++) {
            M[k] = (M[k] | (two[k] & ~three[k])) & ~three[k];
            if (Tm.w[k] & ~M[k]) complete = false;
            if (Sm.w[k] & M[k]) clean = false;
        }
        if (complete && clean) {
            if (++matches >= cnt) return 0;
            for (int k = 0; k < W; k++) M[k] &= ~Tm.w[k];
        }
    }
    int c = 0;
    for (int k = 0; k < W; k++) c += __builtin_popcountll(Tm.w[k] & ~M[k]) + __builtin_popcountll(Sm.w[k] & M[k]);
    return (cnt - matches) * bigPenalty + c;
}

bool anneal(bool warm) {
    if (!warm) for (int t = 0; t < 6; t++) for (int i = 0; i < K; i++) seq[t][i] = 0;
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
        else if (r == 2 && i > 0 && seq[t][i - 1]) nv = seq[t][i - 1];
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
    unsigned long long seed;
    if (scanf("%d %d %d %d %lld %llu", &L, &n, &cnt, &kmax, &iters, &seed) != 6) return 1;
    rng.seed(seed);
    kmax = std::min(kmax, MAXK);
    for (int t = 0; t < 6; t++) if (scanf("%d", &startp[t]) != 1) return 1;
    for (int t = 0; t < 6; t++) {
        memset(&bm[t][0], 0, sizeof(Mask));
        pos[t][0] = startp[t];
        for (int i = 1; i <= L; i++) {
            if (scanf("%d", &pos[t][i]) != 1) return 1;
            bm[t][i] = rd();
        }
    }
    Tm = rd(); Sm = rd();
    bigPenalty = 0;
    for (int k = 0; k < W; k++) bigPenalty += __builtin_popcountll(Tm.w[k]) + __builtin_popcountll(Sm.w[k]);
    bigPenalty++;
    // binary search on K (annealing is randomised, so this is approximate)
    int lo = 1, hi = kmax, bestK = -1;
    while (lo <= hi) {
        K = (lo + hi) / 2;
        bool ok = false;
        for (int restart = 0; restart < 2 && !ok; restart++) ok = anneal(false);
        if (ok) {
            bestK = K;
            for (int t = 0; t < 6; t++) for (int i = 0; i < K; i++) best[t][i] = seq[t][i];
            hi = K - 1;
        } else {
            lo = K + 1;
        }
    }
    if (bestK < 0) { printf("none\n"); return 0; }
    printf("%d\n", bestK);
    for (int t = 0; t < 6; t++) {
        for (int i = 0; i < bestK; i++) printf("%d ", best[t][i]);
        printf("\n");
    }
    return 0;
}
