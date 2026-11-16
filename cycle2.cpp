// Exhaustive search for a 2-tick cycle X,Y (X->Y and Y->X reachable in one
// tick) that rebuilds the pattern every cycle.
// stdin:  L  start[6]  then for each emitter t: L lines of "pos mask"
//         then T mask, S mask   (masks as decimal unsigned 128-bit halves: hi lo)
// stdout: "none" or "X0..X5 Y0..Y5" as state indices (0 = off, 1..L = line)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
typedef unsigned __int128 u128;

int L, n, cnt;
int startp[6];
int pos[6][32];
u128 bm[6][32];  // bm[t][0] = 0 (off)
u128 Tm, Sm;
int bestCost = 1 << 30;
int bestX[6], bestY[6];
int curX[6], curY[6];
long long nodes = 0;

static u128 rd() {
    unsigned long long hi, lo;
    if (scanf("%llu %llu", &hi, &lo) != 2) exit(1);
    return ((u128)hi << 64) | lo;
}

int ringDist(int a, int b) {
    int d = ((b - a) % n + n) % n;
    return std::min(d, n - d);
}

// does an X state pair with a Y state for emitter t in one tick both ways?
bool linked(int x, int y) {
    if (x == 0 || y == 0) return true;      // toggling (off emitter sits on the line)
    return abs(x - y) <= 1;                 // same line or one step
}

void dfs(int t, u128 x1, u128 x2, u128 x3, u128 y1, u128 y2, u128 y3) {
    nodes++;
    // a target hit 3+ times in Y is destroyed every cycle
    if (y3 & Tm) return;
    int rem = 6 - t;
    if (rem == 0) {
        u128 cy = y2 & ~y3, cx = x2 & ~x3;
        if (cy & Sm) return;
        if ((cx & Sm) & ~y3) return;
        if (Tm & ~(cy | (cx & ~y3))) return;
        int cost = 0;
        bool same = true;
        for (int e = 0; e < 6; e++) {
            int s = curX[e] ? curX[e] : curY[e];
            if (s) cost = std::max(cost, ringDist(startp[e], pos[e][s]));
            if (curX[e] != curY[e]) same = false;
        }
        cost += cnt * (same ? 1 : 2);
        if (cost < bestCost) {
            bestCost = cost;
            for (int e = 0; e < 6; e++) { bestX[e] = curX[e]; bestY[e] = curY[e]; }
        }
        return;
    }
    if (rem == 1) {
        // every target needs >=1 hit already in Y or X (one more beam can add at most 1)
        if (Tm & ~(y1 | x1)) return;
    }
    for (int y = 0; y <= L; y++) {
        u128 b = bm[t][y];
        u128 ny3 = y3 | (y2 & b), ny2 = y2 | (y1 & b), ny1 = y1 | b;
        if (ny3 & Tm) continue;
        for (int x = 0; x <= L; x++) {
            if (!linked(x, y)) continue;
            u128 c = bm[t][x];
            curX[t] = x; curY[t] = y;
            dfs(t + 1, x1 | c, x2 | (x1 & c), x3 | (x2 & c), ny1, ny2, ny3);
        }
    }
}

int main() {
    if (scanf("%d %d %d", &L, &n, &cnt) != 3) return 1;
    for (int t = 0; t < 6; t++) if (scanf("%d", &startp[t]) != 1) return 1;
    for (int t = 0; t < 6; t++) {
        bm[t][0] = 0; pos[t][0] = startp[t];
        for (int i = 1; i <= L; i++) {
            if (scanf("%d", &pos[t][i]) != 1) return 1;
            bm[t][i] = rd();
        }
    }
    Tm = rd(); Sm = rd();
    dfs(0, 0, 0, 0, 0, 0, 0);
    if (bestCost == (1 << 30)) { printf("none\n"); return 0; }
    for (int e = 0; e < 6; e++) printf("%d ", bestX[e]);
    for (int e = 0; e < 6; e++) printf("%d ", bestY[e]);
    printf("\n");
    return 0;
}
