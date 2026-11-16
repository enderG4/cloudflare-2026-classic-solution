// Shortens a working program: repeatedly deletes one tick, then repairs the
// program with local search (shifting, adding and removing moves/toggles)
// until it produces the pattern `count` times again.
// stdin:  n K count budget_ms seed
//         start[6]
//         for each emitter t, for each ring position p: beam mask (7 words)
//         T mask, S mask
//         K lines of 6 instructions (F B T -)
// stdout: the shortest working program found, one line per tick
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

const int W = 7;
struct Mask {
    uint64_t w[W];
    Mask() { memset(w, 0, sizeof w); }
};
static inline Mask band(const Mask& a, const Mask& b) { Mask r; for (int i = 0; i < W; i++) r.w[i] = a.w[i] & b.w[i]; return r; }
static inline int popc(const Mask& a) { int c = 0; for (int i = 0; i < W; i++) c += __builtin_popcountll(a.w[i]); return c; }

int n, K, cnt;
int startp[6];
std::vector<Mask> bm[6];
Mask Tm, Sm;
int bigPenalty;

enum { NOP = 0, FWD = 1, BWD = 2, TOG = 3 };
std::vector<int> prog[6];

struct State {
    int pos[6];
    bool act[6];
    Mask M;
    int matches;
};
std::vector<State> cache;  // cache[i] = state before tick i
int matchTick;             // first tick after which matches >= cnt (or -1)

static Mask readMask() {
    Mask m;
    for (int i = 0; i < W; i++) {
        unsigned long long x;
        if (scanf("%llu", &x) != 1) exit(1);
        m.w[i] = x;
    }
    return m;
}

// simulate from tick i0 (cache[i0] must be valid); returns error, 0 = works
int evalFrom(int i0) {
    State s = cache[i0];
    matchTick = -1;
    for (int i = 0; i < i0; i++) if (cache[i + 1].matches >= cnt) { matchTick = i; break; }
    if (matchTick >= 0) return 0;
    for (int i = i0; i < K; i++) {
        for (int t = 0; t < 6; t++) {
            int ins = prog[t][i];
            if (ins == FWD) s.pos[t] = (s.pos[t] + 1) % n;
            else if (ins == BWD) s.pos[t] = (s.pos[t] + n - 1) % n;
            else if (ins == TOG) s.act[t] = !s.act[t];
        }
        uint64_t one[W] = {0}, two[W] = {0}, three[W] = {0};
        for (int t = 0; t < 6; t++) {
            if (!s.act[t]) continue;
            const Mask& b = bm[t][s.pos[t]];
            for (int k = 0; k < W; k++) {
                three[k] |= two[k] & b.w[k];
                two[k] |= one[k] & b.w[k];
                one[k] |= b.w[k];
            }
        }
        bool complete = true, clean = true;
        for (int k = 0; k < W; k++) {
            s.M.w[k] = (s.M.w[k] | (two[k] & ~three[k])) & ~three[k];
            if (Tm.w[k] & ~s.M.w[k]) complete = false;
            if (Sm.w[k] & s.M.w[k]) clean = false;
        }
        if (complete && clean) {
            s.matches++;
            for (int k = 0; k < W; k++) s.M.w[k] &= ~Tm.w[k];
        }
        cache[i + 1] = s;
        if (s.matches >= cnt) { matchTick = i; return 0; }
    }
    Mask miss, junk;
    for (int k = 0; k < W; k++) { miss.w[k] = Tm.w[k] & ~s.M.w[k]; junk.w[k] = Sm.w[k] & s.M.w[k]; }
    return (cnt - s.matches) * bigPenalty + popc(miss) + popc(junk);
}

std::mt19937_64 rng;
static inline int rnd(int m) { return (int)(rng() % (uint64_t)m); }

struct Move { int t, i, j, oi, oj; };

// propose a local change; returns first affected tick or -1
int propose(Move& mv) {
    int t = rnd(6);
    int kind = rnd(4);
    mv.t = t;
    if (kind == 0) {  // swap neighbours in time
        int i = rnd(K - 1);
        if (prog[t][i] == prog[t][i + 1]) return -1;
        mv.i = i; mv.j = i + 1;
    } else if (kind == 1) {  // move one instruction to a free slot nearby
        int i = rnd(K);
        if (prog[t][i] == NOP) return -1;
        int j = i + (rnd(2) ? 1 : -1) * (1 + rnd(4));
        if (j < 0 || j >= K || prog[t][j] != NOP) return -1;
        mv.i = i; mv.j = j;
    } else if (kind == 2) {  // insert a pair (detour or on/off blink)
        int i = rnd(K);
        int j = i + 1 + rnd(6);
        if (j >= K || prog[t][i] != NOP || prog[t][j] != NOP) return -1;
        mv.i = i; mv.j = j;
        mv.oi = prog[t][i]; mv.oj = prog[t][j];
        int r = rnd(3);
        prog[t][i] = r == 0 ? FWD : r == 1 ? BWD : TOG;
        prog[t][j] = r == 0 ? BWD : r == 1 ? FWD : TOG;
        return i;
    } else {  // remove a cancelling pair
        int i = rnd(K);
        int a = prog[t][i];
        if (a == NOP) return -1;
        int want = a == FWD ? BWD : a == BWD ? FWD : TOG;
        int j = -1;
        for (int d = 1; d <= 8 && i + d < K; d++) if (prog[t][i + d] == want) { j = i + d; break; }
        if (j < 0) return -1;
        mv.i = i; mv.j = j; mv.oi = prog[t][i]; mv.oj = prog[t][j];
        prog[t][i] = NOP; prog[t][j] = NOP;
        return i;
    }
    mv.oi = prog[mv.t][mv.i]; mv.oj = prog[mv.t][mv.j];
    std::swap(prog[mv.t][mv.i], prog[mv.t][mv.j]);
    return std::min(mv.i, mv.j);
}

void undo(const Move& mv) { prog[mv.t][mv.i] = mv.oi; prog[mv.t][mv.j] = mv.oj; }

// local search at fixed length; true if it reaches a working program
bool repair(long long maxIters) {
    int cur = evalFrom(0);
    if (cur == 0) return true;
    double T0 = 1.5, T1 = 0.05;
    for (long long it = 0; it < maxIters; it++) {
        Move mv;
        int i0 = propose(mv);
        if (i0 < 0) continue;
        int e = evalFrom(i0);
        double temp = T0 * pow(T1 / T0, (double)it / maxIters);
        if (e <= cur || std::uniform_real_distribution<double>(0, 1)(rng) < exp((cur - e) / temp)) {
            cur = e;
            if (cur == 0) return true;
        } else {
            undo(mv);
            evalFrom(i0);
        }
    }
    return false;
}

// delete tick d; its non-noop instructions move to the nearest free slot
bool deleteTick(int d) {
    for (int t = 0; t < 6; t++) {
        int ins = prog[t][d];
        prog[t].erase(prog[t].begin() + d);
        if (ins == NOP) continue;
        int placed = -1;
        for (int r = 0; r < K - 1 && placed < 0; r++) {
            for (int j : {d - 1 - r, d + r}) {
                if (j >= 0 && j < K - 1 && prog[t][j] == NOP) { placed = j; break; }
            }
        }
        if (placed < 0) return false;
        prog[t][placed] = ins;
    }
    K--;
    cache.resize(K + 1);
    return true;
}

int main() {
    long long budgetMs;
    unsigned long long seed;
    if (scanf("%d %d %d %lld %llu", &n, &K, &cnt, &budgetMs, &seed) != 5) return 1;
    rng.seed(seed);
    for (int t = 0; t < 6; t++) if (scanf("%d", &startp[t]) != 1) return 1;
    for (int t = 0; t < 6; t++) { bm[t].resize(n); for (int p = 0; p < n; p++) bm[t][p] = readMask(); }
    Tm = readMask(); Sm = readMask();
    bigPenalty = popc(Tm) + popc(Sm) + 1;
    for (int t = 0; t < 6; t++) prog[t].resize(K);
    for (int i = 0; i < K; i++) {
        for (int t = 0; t < 6; t++) {
            char c[4];
            if (scanf("%3s", c) != 1) return 1;
            prog[t][i] = c[0] == 'F' ? FWD : c[0] == 'B' ? BWD : c[0] == 'T' ? TOG : NOP;
        }
    }
    cache.resize(K + 1);
    State s0;
    for (int t = 0; t < 6; t++) { s0.pos[t] = startp[t]; s0.act[t] = false; }
    s0.matches = 0;
    cache[0] = s0;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count(); };

    if (evalFrom(0) != 0) { fprintf(stderr, "initial program does not work\n"); return 1; }
    auto truncate = [&]() {
        evalFrom(0);
        if (matchTick >= 0 && matchTick + 1 < K) {
            K = matchTick + 1;
            for (int t = 0; t < 6; t++) prog[t].resize(K);
            cache.resize(K + 1);
        }
    };
    truncate();
    std::vector<int> best[6];
    for (int t = 0; t < 6; t++) best[t] = prog[t];
    int bestK = K;

    int fails = 0;
    while (elapsed() < budgetMs && K > 1 && fails < 4 * K) {
        int d = rnd(K);
        if (!deleteTick(d) || !repair(20000 + 2000LL * K)) {
            fails++;
            K = bestK;
            for (int t = 0; t < 6; t++) prog[t] = best[t];
            cache.resize(K + 1);
            continue;
        }
        fails = 0;
        truncate();
        bestK = K;
        for (int t = 0; t < 6; t++) best[t] = prog[t];
    }
    const char* sym = "-FBT";
    for (int i = 0; i < bestK; i++) {
        for (int t = 0; t < 6; t++) printf("%c%c", sym[best[t][i]], t < 5 ? ' ' : '\n');
    }
    return 0;
}
