import os
import subprocess
import sys
from itertools import product
from multiprocessing import Pool

from solve4 import D, build_ring, add
from solve6 import parse_patterns, board_lines_to_matter, build_program_for_pattern

DIRS = [D[k] for k in range(1, 7)]
# axis (cube index) that stays constant along emitter t's beam
AX = [next(i for i in range(3) if d[i] == 0) for d in DIRS]
PAIR = {a: [t for t in range(6) if AX[t] == a] for a in range(3)}
MAX_TRIMS = 4
LOOKAHEAD = 4


class Geo:
    def __init__(self, R):
        self.R = R
        self.n = 6 * R
        self.ring = build_ring(R)
        self.start = [t * R for t in range(6)]
        self.cells = [(q, r, -q - r) for q in range(-R + 1, R) for r in range(-R + 1, R)
                      if abs(q + r) < R]
        self.line = {}
        for a in range(3):
            for v in range(-R + 1, R):
                self.line[(a, v)] = [c for c in self.cells if c[a] == v]
        self.beam = [[() for _ in range(self.n)] for _ in range(6)]
        self.idx = {c: i for i, c in enumerate(self.cells)}
        self.bmask = [[0] * self.n for _ in range(6)]
        self.pos_for = [dict() for _ in range(6)]
        for t in range(6):
            for p in range(self.n):
                cells = []
                c = add(self.ring[p], DIRS[t])
                while max(map(abs, c)) < R:
                    cells.append(c)
                    c = add(c, DIRS[t])
                self.beam[t][p] = tuple(cells)
                self.bmask[t][p] = sum(1 << self.idx[x] for x in cells)
                if cells:
                    v = cells[0][AX[t]]
                    assert v not in self.pos_for[t]
                    self.pos_for[t][v] = p

    def mask(self, cells):
        return sum(1 << self.idx[x] for x in cells if x in self.idx)

    def dist(self, a, b):
        d = (b - a) % self.n
        return min(d, self.n - d)

    def step(self, a, b):
        return "F" if (b - a) % self.n <= (a - b) % self.n else "B"


class Sim:
    def __init__(self, geo, T, S):
        self.g = geo
        self.T = T
        self.S = S
        self.pos = list(geo.start)
        self.act = [False] * 6
        self.M = set()
        self.accepted = 0

    def copy(self):
        s = Sim.__new__(Sim)
        s.g, s.T, s.S = self.g, self.T, self.S
        s.pos, s.act, s.M, s.accepted = list(self.pos), list(self.act), set(self.M), self.accepted
        return s

    def tick(self, line, protected=None):
        """Apply one instruction line. Returns False if a write lands on an
        unprotected cell in a way that breaks the final pattern."""
        n = self.g.n
        for t, ins in enumerate(line):
            if ins == "F":
                self.pos[t] = (self.pos[t] + 1) % n
            elif ins == "B":
                self.pos[t] = (self.pos[t] - 1) % n
            elif ins == "T":
                self.act[t] = not self.act[t]
        hits = {}
        for t in range(6):
            if self.act[t]:
                for c in self.g.beam[t][self.pos[t]]:
                    hits[c] = hits.get(c, 0) + 1
        ok = True
        for c, k in hits.items():
            if k == 2:
                if c not in self.M:
                    if protected is not None and c in self.S and c not in protected:
                        ok = False
                    self.M.add(c)
            elif k >= 3 and c in self.M:
                if protected is not None and c in self.T and c not in protected:
                    ok = False
                self.M.discard(c)
        if self.T <= self.M and not (self.S & self.M):
            self.accepted += 1
            self.M -= self.T
        return ok


def neighbors(c):
    return [add(c, d) for d in DIRS]


def side_effect_cells(a, v, bset, cset, b, c):
    out = []
    for fb in bset:
        for fc in cset:
            vb, vc = fb[b], fc[c]
            x = [0, 0, 0]
            x[b], x[c], x[a] = vb, vc, -vb - vc
            x = tuple(x)
            if x[a] != v:
                out.append(x)
    return out


def trim_options(geo, a, v, need, S, later):
    """Ways to split the cells to trim between the two other axes, such that
    crossings of the trimmer beams don't create unprotected skirt matter."""
    b, c = [x for x in range(3) if x != a]
    R = geo.R
    opts = []
    k = len(need)
    for mask in range(1 << k):
        bset = [need[i] for i in range(k) if mask >> i & 1]
        cset = [need[i] for i in range(k) if not mask >> i & 1]
        if len(bset) > 2 or len(cset) > 2:
            continue
        bad = False
        for x in side_effect_cells(a, v, bset, cset, b, c):
            if max(map(abs, x)) < R and x in S and x not in later:
                bad = True
                break
        if not bad:
            opts.append((b, bset, c, cset))
    # prefer fewer side effects (all trims on one axis) first
    opts.sort(key=lambda o: min(len(o[1]), len(o[3])))
    return opts


def peel(geo, T, S, primary=None, min_new=1, strict=False):
    """Reverse-time greedy: the last fill on a line decides all of its cells.
    Returns line ops in forward order plus targets no line could cover."""
    undecided = set(T) | set(S)
    picked = []
    lines = list(geo.line.keys())
    while True:
        best = None
        for (a, v) in lines:
            if (a, v) in [p[0] for p in picked]:
                continue
            cells = geo.line[(a, v)]
            if strict and a != primary:
                continue
            newT = sum(1 for x in cells if x in T and x in undecided)
            if newT < min_new:
                continue
            need = [x for x in cells if x in S and x in undecided]
            if len(need) > MAX_TRIMS:
                continue
            later = set(geo.cells) - undecided
            if not trim_options(geo, a, v, need, S, later):
                continue
            score = (newT + (2 if a == primary else 0), -len(need))
            if best is None or score > best[0]:
                best = (score, (a, v))
        if best is None and min_new == 1:
            best = pick_cleaner(geo, T, S, undecided, picked, primary if strict else None)
        if best is None:
            break
        key = best[1]
        picked.append((key,))
        undecided -= set(geo.line[key])
    order = [p[0] for p in reversed(picked)]
    uncovered = [x for x in T if x in undecided]
    return order, uncovered


def feasible_targets(geo, T, S, undecided, skip, axis=None):
    tot = 0
    for key, cells in geo.line.items():
        if key in skip or (axis is not None and key[0] != axis):
            continue
        newT = sum(1 for x in cells if x in T and x in undecided)
        if newT and sum(1 for x in cells if x in S and x in undecided) <= MAX_TRIMS:
            tot += newT
    return tot


def pick_cleaner(geo, T, S, undecided, picked, axis=None):
    """When no target line is usable, pick a line that only clears skirt cells
    (processed later in time) so that more target lines become usable."""
    skip = {p[0] for p in picked}
    later = set(geo.cells) - undecided
    best = None
    for key, cells in geo.line.items():
        if key in skip:
            continue
        need = [x for x in cells if x in S and x in undecided]
        if len(need) > MAX_TRIMS or (not need and axis is None):
            continue
        if not trim_options(geo, key[0], key[1], need, S, later):
            continue
        gain = feasible_targets(geo, T, S, undecided - set(cells), skip | {key}, axis)
        if axis is not None:
            gain -= feasible_targets(geo, T, S, undecided, skip, axis)
        if gain and (best is None or gain > best[0]):
            best = (gain, key)
    return best


def assign_trims(geo, order, T, S):
    """For a forward order of line ops compute the trim options per op; None if invalid."""
    res = []
    later = set()
    for key in reversed(order):
        a, v = key
        cells = geo.line[key]
        need = [x for x in cells if x in S and x not in later]
        if len(need) > MAX_TRIMS:
            return None
        opts = trim_options(geo, a, v, need, S, later)
        if not opts:
            return None
        res.append((key, opts))
        later |= set(cells)
    return list(reversed(res))


def config_for_line(geo, op, cur):
    """Choose concrete emitter placement for a line op; returns (cost, desired)."""
    (a, v), opts = op
    best = None
    for (b, bset, c, cset) in opts:
        desired = [None] * 6
        for t in PAIR[a]:
            desired[t] = geo.pos_for[t][v]
        choices = []
        for ax, cs in ((b, bset), (c, cset)):
            em = PAIR[ax]
            if len(cs) == 0:
                choices.append([{}])
            elif len(cs) == 1:
                choices.append([{em[0]: geo.pos_for[em[0]][cs[0][ax]]},
                                {em[1]: geo.pos_for[em[1]][cs[0][ax]]}])
            else:
                choices.append([{em[0]: geo.pos_for[em[0]][cs[0][ax]], em[1]: geo.pos_for[em[1]][cs[1][ax]]},
                                {em[0]: geo.pos_for[em[0]][cs[1][ax]], em[1]: geo.pos_for[em[1]][cs[0][ax]]}])
        for cb, cc in product(*choices):
            d = list(desired)
            for t, p in list(cb.items()) + list(cc.items()):
                d[t] = p
            ds = [geo.dist(cur[t], d[t]) for t in range(6) if d[t] is not None]
            cost = (max(ds), sum(ds))
            if best is None or cost < best[0]:
                best = (cost, d)
    return best


def config_for_point(geo, x, cur):
    best = None
    for t1 in range(6):
        for t2 in range(t1 + 1, 6):
            if AX[t1] == AX[t2]:
                continue
            d = [None] * 6
            d[t1] = geo.pos_for[t1][x[AX[t1]]]
            d[t2] = geo.pos_for[t2][x[AX[t2]]]
            d1, d2 = geo.dist(cur[t1], d[t1]), geo.dist(cur[t2], d[t2])
            cost = (max(d1, d2), d1 + d2)
            if best is None or cost < best[0]:
                best = (cost, d)
    return best


def transition_lines(geo, sim, desired, keep, hint=None, force_hint=False):
    """Instruction lines moving from sim's state to `desired` (None = inactive).
    Emitters in `keep` that are active stay active while moving; all other
    active emitters are switched off first. Idle emitters drift toward `hint`
    without lengthening the transition."""
    seqs = []
    idle = []
    for t in range(6):
        seq = []
        p = sim.pos[t]
        on = sim.act[t]
        want_on = desired[t] is not None
        if on and not (want_on and t in keep):
            seq.append("T")
            on = False
        if want_on:
            target = desired[t]
            while p != target:
                ins = geo.step(p, target)
                seq.append(ins)
                p = (p + (1 if ins == "F" else -1)) % geo.n
        elif hint is not None and hint[t] is not None:
            idle.append((t, seq, p))
        seqs.append((seq, want_on and not on))
    total = max(1, max(len(s) + (1 if late else 0) for s, late in seqs))
    if force_hint:
        total = max([total] + [len(seq) + geo.dist(p, hint[t]) for t, seq, p in idle])
    for t, seq, p in idle:
        target = hint[t]
        while p != target and len(seq) < total:
            ins = geo.step(p, target)
            seq.append(ins)
            p = (p + (1 if ins == "F" else -1)) % geo.n
    lines = [["-"] * 6 for _ in range(total)]
    for t, (seq, late) in enumerate(seqs):
        for i, ins in enumerate(seq):
            lines[i][t] = ins
        if late:
            lines[total - 1][t] = "T"
    return lines


def run_op(geo, sim, desired, protected, fill_axis, hint=None, force_hint=False):
    """Try transitions from cheapest to safest; apply the first harmless one."""
    active_now = [t for t in range(6) if sim.act[t]]
    cands = []
    for mask in range(1 << len(active_now)):
        keep = {active_now[i] for i in range(len(active_now)) if mask >> i & 1}
        lines = transition_lines(geo, sim, desired, keep, hint, force_hint)
        cands.append((len(lines), -len(keep), lines))
    cands.sort(key=lambda c: (c[0], c[1]))
    cands = [c[2] for c in cands]
    for lines in cands:
        trial = sim.copy()
        if all(trial.tick(l, protected) for l in lines):
            best_len = len(lines)
            break
    else:
        best_len = len(cands[-1]) + 1
    lines = None if force_hint else schedule_search(geo, sim, desired, protected, best_len)
    if lines is not None:
        trial = sim.copy()
        if all(trial.tick(l, protected) for l in lines):
            cands.insert(0, lines)
    for lines in cands:
        trial = sim.copy()
        if all(trial.tick(l, protected) for l in lines):
            return trial, lines
    # the all-off transition is always harmless; reaching here means the op itself is unsafe
    trial = sim.copy()
    for l in cands[-1]:
        trial.tick(l)
    return trial, cands[-1]


def schedule_search(geo, sim, desired, protected, max_len, budget=4000):
    """Find a transition shorter than max_len by letting each emitter wait,
    step toward its target or toggle on any tick, so that no intermediate
    tick writes matter onto an unprotected skirt cell (or clears a target)."""
    n = geo.n
    want_act = [d is not None for d in desired]

    def need(t, p, a):
        k = 0
        if desired[t] is not None:
            k += geo.dist(p, desired[t])
        return k + (1 if a != want_act[t] else 0)

    lb = max(1, max(need(t, sim.pos[t], sim.act[t]) for t in range(6)))
    if lb >= max_len:
        return None
    key = (id(protected), max_len)
    cache = getattr(sim, "_mask_cache", None)
    if cache is None or cache[0] != id(protected):
        bad_s = geo.mask(sim.S) & ~geo.mask(protected)
        bad_t = geo.mask(sim.T) & ~geo.mask(protected)
        sim._mask_cache = (id(protected), bad_s, bad_t)
    _, bad_s, bad_t = sim._mask_cache
    bm = geo.bmask
    nodes = [0]

    def harmless(pos, act):
        one = two = three = 0
        for t in range(6):
            if act[t]:
                b = bm[t][pos[t]]
                three |= two & b
                two |= one & b
                one |= b
        return not ((two & ~three & bad_s) or (three & bad_t))

    def options(t, p, a, r):
        out = []
        if desired[t] is not None and p != desired[t]:
            s = geo.step(p, desired[t])
            np_ = (p + (1 if s == "F" else -1)) % n
            if need(t, np_, a) <= r:
                out.append((s, np_, a))
        if a != want_act[t] and need(t, p, not a) <= r:
            out.append(("T", p, not a))
        if need(t, p, a) <= r:
            out.append(("-", p, a))
        return out

    for L in range(lb, max_len):
        failed = set()

        def dfs(k, pos, act):
            if k == L:
                return []
            st = (k, pos, act)
            if st in failed or nodes[0] > budget:
                return None
            r = L - k - 1
            opts = [options(t, pos[t], act[t], r) for t in range(6)]
            if any(not o for o in opts):
                failed.add(st)
                return None
            for combo in product(*opts):
                nodes[0] += 1
                if nodes[0] > budget:
                    return None
                npos = tuple(c[1] for c in combo)
                nact = tuple(c[2] for c in combo)
                if k < L - 1 and not harmless(npos, nact):
                    continue
                rest = dfs(k + 1, npos, nact)
                if rest is not None:
                    return [[c[0] for c in combo]] + rest
            failed.add(st)
            return None

        res = dfs(0, tuple(sim.pos), tuple(sim.act))
        if res is not None:
            return res
        if nodes[0] > budget:
            return None
    return None


def nn_order(geo, order, opts_by_key, start_pos):
    """Reorder line ops greedily by travel cost, keeping the peel's precedence edges."""
    idx = {k: i for i, k in enumerate(order)}
    must_before = {k: set() for k in order}
    later_cells = []
    acc = set()
    for k in reversed(order):
        later_cells.append(set(acc))
        acc |= set(geo.line[k])
    later_cells.reverse()
    for i, k in enumerate(order):
        for j in range(i + 1, len(order)):
            k2 = order[j]
            if set(geo.line[k]) & set(geo.line[k2]):
                must_before[k2].add(k)
    done = []
    cur = list(start_pos)
    remaining = list(order)
    while remaining:
        avail = [k for k in remaining if must_before[k] <= set(done)]
        best = None
        for k in avail:
            cost, d = config_for_line(geo, (k, opts_by_key[k]), cur)
            if best is None or cost < best[0]:
                best = (cost, k, d)
        _, k, d = best
        done.append(k)
        remaining.remove(k)
        cur = [d[t] if d[t] is not None else cur[t] for t in range(6)]
    return done


def reverse_nn_order(geo, lines, T, S, seed_cost):
    """Build a forward order back to front: each step prepends the valid line
    whose emitter placement is closest to the one that follows it."""
    later = set()
    placed = []
    cur = None
    remaining = list(lines)
    while remaining:
        best = None
        for key in remaining:
            need = [x for x in geo.line[key] if x in S and x not in later]
            if len(need) > MAX_TRIMS:
                continue
            opts = trim_options(geo, key[0], key[1], need, S, later)
            if not opts:
                continue
            if cur is None:
                cost = seed_cost(key)
                d = config_for_line(geo, (key, opts), geo.start)[1]
            else:
                cost, d = config_for_line(geo, (key, opts), cur)
            if best is None or cost < best[0]:
                best = (cost, key, d)
        if best is None:
            return None
        _, key, d = best
        placed.append(key)
        remaining.remove(key)
        later |= set(geo.line[key])
        cur = [d[t] if d[t] is not None else (cur[t] if cur else geo.start[t]) for t in range(6)]
    placed.reverse()
    return placed


def build_rep_ops(geo, order, T, S):
    trims = assign_trims(geo, order, T, S)
    return trims


def compile_program(geo, T, S, orders, uncovered, count, pen_at="start"):
    sim = Sim(geo, T, S)
    program = []
    guard = 0
    while sim.accepted < count:
        guard += 1
        if guard > count + 5:
            return None
        start_acc = sim.accepted
        # choose the cheapest valid order for this rep from current position
        best = None
        for order in orders:
            ops = assign_trims(geo, order, T, S)
            if ops is None:
                continue
            c0, _ = config_for_line(geo, ops[0], sim.pos) if ops else (0, None)
            if best is None or c0 < best[0]:
                best = (c0, ops)
        ops = best[1] if best else []
        seq = [("L", op) for op in ops]
        if uncovered:
            seq.insert(0 if pen_at == "start" else len(seq), ("PEN", uncovered))
        suffix = [set()]
        for kind, item in reversed(seq):
            suffix.append(suffix[-1] | (set(item) if kind == "PEN" else set(geo.line[item[0]])))
        suffix.reverse()

        def config(kind, item, cur):
            if kind == "PEN":
                return None
            return config_for_line(geo, item, cur)[1]

        for i, (kind, item) in enumerate(seq):
            if kind == "PEN":
                sim, lines = run_pen(geo, sim, item, suffix[i + 1])
                if sim is None:
                    return None
            else:
                desired = config(kind, item, sim.pos)
                keep_spare_active(geo, sim, desired, suffix[i + 1])
                hint = None
                if i + 1 < len(seq):
                    # idle emitters head for where they are next needed (a few ops ahead)
                    hint = [None] * 6
                    after = [desired[t] if desired[t] is not None else sim.pos[t] for t in range(6)]
                    for j in range(i + 1, min(i + 1 + LOOKAHEAD, len(seq))):
                        nxt = config(seq[j][0], seq[j][1], after)
                        if nxt is None:
                            break
                        for t in range(6):
                            if desired[t] is None and hint[t] is None and nxt[t] is not None:
                                hint[t] = nxt[t]
                        after = [nxt[t] if nxt[t] is not None else after[t] for t in range(6)]
                sim, lines = run_op(geo, sim, desired, suffix[i], item[0][0], hint)
            program.extend(lines)
            if sim.accepted >= count:
                return program
            if sim.accepted > start_acc:
                break
        else:
            if sim.accepted == start_acc:
                return None
    return program


def keep_spare_active(geo, sim, desired, later):
    """Emitters the op doesn't need but that are on right now stay on where
    they are, if their beams don't spoil anything at the op's final tick."""
    later_m = geo.mask(later)
    bad_s = geo.mask(sim.S) & ~later_m
    bad_t = geo.mask(sim.T) & ~later_m

    def ok(d):
        one = two = three = 0
        for t in range(6):
            if d[t] is not None:
                b = geo.bmask[t][d[t]]
                three |= two & b
                two |= one & b
                one |= b
        return not ((two & ~three & bad_s) or (three & bad_t))

    for t in range(6):
        if desired[t] is None and sim.act[t]:
            desired[t] = sim.pos[t]
            if not ok(desired):
                desired[t] = None


def pen_tour(geo, a, b, targets, S, protected, start, tiebreak=True):
    """Greedy tour of the crossing point of an a-beam and a b-beam over all
    targets; each tick both beams shift by at most one line, so the crossing
    moves to a neighbouring cell (or two cells diagonally). Skirt cells that
    no later op rewrites are walls. Among equally near targets it prefers the
    one with the fewest unvisited neighbours, so dead ends get finished first."""
    R = geo.R
    c = 3 - a - b
    lo, hi = -R + 1, R - 1
    moves = [(da, db) for da in (-1, 0, 1) for db in (-1, 0, 1) if da or db]
    ok = {}
    for va in range(lo, hi + 1):
        for vb in range(lo, hi + 1):
            x = [0, 0, 0]
            x[a], x[b], x[c] = va, vb, -va - vb
            x = tuple(x)
            ok[(va, vb)] = max(map(abs, x)) >= R or x not in S or x in protected

    todo = {(t[a], t[b]) for t in targets}

    def degree(s):
        return sum(1 for da, db in moves if (s[0] + da, s[1] + db) in todo)

    cur = start
    path = [cur]
    todo.discard(cur)
    while todo:
        prev = {cur: None}
        frontier = [cur]
        found = []
        while frontier and not found:
            nxt = []
            for s0 in frontier:
                for da, db in moves:
                    s = (s0[0] + da, s0[1] + db)
                    if s in prev or not ok.get(s, False):
                        continue
                    prev[s] = s0
                    if s in todo:
                        found.append(s)
                    nxt.append(s)
            frontier = nxt
        if not found:
            return None
        tgt = min(found, key=degree) if tiebreak else found[0]
        seg = []
        s = tgt
        while s != cur:
            seg.append(s)
            s = prev[s]
        seg.reverse()
        for s in seg:
            todo.discard(s)
        path.extend(seg)
        cur = tgt
    return path


def pen_plan(geo, targets, S, protected, cur, extra_starts=3):
    best = None
    for a in range(3):
        for b in range(a + 1, 3):
            cells = {(t[a], t[b]) for t in targets}
            combos = [(ea, eb) for ea in PAIR[a] for eb in PAIR[b]]

            def entry(ea, eb, s):
                return max(geo.dist(cur[ea], geo.pos_for[ea][s[0]]),
                           geo.dist(cur[eb], geo.pos_for[eb][s[1]]))

            starts = {min(cells, key=lambda s: entry(ea, eb, s)) for ea, eb in combos}
            deg = lambda s: sum(1 for da in (-1, 0, 1) for db in (-1, 0, 1)
                                if (da or db) and (s[0] + da, s[1] + db) in cells)
            starts |= set(sorted(cells, key=deg)[:extra_starts])
            for s in starts:
                path = pen_tour(geo, a, b, targets, S, protected, s)
                if path is None:
                    continue
                for ea, eb in combos:
                    cost = entry(ea, eb, s) + len(path)
                    if best is None or cost < best[0]:
                        best = (cost, ea, eb, path)
    return best


def run_pen(geo, sim, targets, protected):
    plan = pen_plan(geo, targets, sim.S, protected, sim.pos)
    if plan is None:
        return None, None
    _, ea, eb, path = plan
    desired = [None] * 6
    desired[ea] = geo.pos_for[ea][path[0][0]]
    desired[eb] = geo.pos_for[eb][path[0][1]]
    sim, lines = run_op(geo, sim, desired, protected, None)
    lines = list(lines)
    for va, vb in path[1:]:
        line = ["-"] * 6
        for t, v in ((ea, va), (eb, vb)):
            p = geo.pos_for[t][v]
            if p != sim.pos[t]:
                line[t] = geo.step(sim.pos[t], p)
        if not sim.tick(line, protected):
            return None, None
        lines.append(line)
    return sim, lines


def verify(geo, T, S, program, count):
    sim = Sim(geo, T, S)
    for l in program:
        sim.tick(l)
    return sim.accepted >= count


CYCLE_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cycle2")
CYCLE_MAX_R = 5


def enc128(x):
    return f"{x >> 64} {x & ((1 << 64) - 1)}"


def cycle_program(geo, T, S, count):
    """Two configurations one tick apart that rebuild the pattern every second
    tick (or every tick when both are equal), found by exhaustive search."""
    if geo.R > CYCLE_MAX_R or not os.path.exists(CYCLE_BIN):
        return None
    R = geo.R
    lines = [f"{2 * R - 1} {geo.n} {count}", " ".join(map(str, geo.start))]
    for t in range(6):
        for v in range(-R + 1, R):
            p = geo.pos_for[t][v]
            lines.append(f"{p} {enc128(geo.bmask[t][p])}")
    lines.append(enc128(geo.mask(T)))
    lines.append(enc128(geo.mask(S)))
    out = subprocess.run([CYCLE_BIN], input="\n".join(lines), capture_output=True, text=True).stdout.split()
    if not out or out[0] == "none":
        return None
    st = list(map(int, out))
    vals = list(range(-R + 1, R))

    def conf(states):
        return [geo.pos_for[t][vals[s - 1]] if s else None for t, s in enumerate(states)]

    X, Y = conf(st[:6]), conf(st[6:])
    sim = Sim(geo, T, S)
    hint = [Y[t] if X[t] is None else None for t in range(6)]
    sim, program = run_op(geo, sim, X, set(), None, hint, force_hint=True)
    program = list(program)
    cur, nxt = X, Y
    while sim.accepted < count and len(program) < 4 * count + 200:
        line = ["-"] * 6
        for t in range(6):
            if (cur[t] is None) != (nxt[t] is None):
                line[t] = "T"
            elif cur[t] is not None and cur[t] != nxt[t]:
                line[t] = geo.step(cur[t], nxt[t])
        sim.tick(line)
        program.append(line)
        cur, nxt = nxt, cur
    return program if sim.accepted >= count else None


CYCLEK_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cyclek")
CYCLEK_ITERS = 300000
CYCLEK_KMAX = 60


def enc448(x):
    return " ".join(str((x >> (64 * i)) & ((1 << 64) - 1)) for i in range(7))


def cyclek_program(geo, T, S, count, kmax):
    """k-tick cycle of emitter states (found by annealing in cyclek) that
    rebuilds the pattern once per cycle; repeated until `count` matches."""
    if kmax < 1 or not os.path.exists(CYCLEK_BIN):
        return None
    R = geo.R
    lines = [f"{2 * R - 1} {geo.n} 1 {kmax} {CYCLEK_ITERS}", " ".join(map(str, geo.start))]
    for t in range(6):
        for v in range(-R + 1, R):
            p = geo.pos_for[t][v]
            lines.append(f"{p} {enc448(geo.bmask[t][p])}")
    lines.append(enc448(geo.mask(T)))
    lines.append(enc448(geo.mask(S)))
    out = subprocess.run([CYCLEK_BIN], input="\n".join(lines), capture_output=True, text=True).stdout.split()
    if not out or out[0] == "none":
        return None
    K = int(out[0])
    vals = list(range(-R + 1, R))
    st = list(map(int, out[1:]))
    seq = [[geo.pos_for[t][vals[s - 1]] if s else None for s in st[t * K:(t + 1) * K]] for t in range(6)]

    def next_on(t, i):
        for d in range(1, K + 1):
            p = seq[t][(i + d) % K]
            if p is not None:
                return p
        return None

    sim = Sim(geo, T, S)
    desired = [seq[t][0] for t in range(6)]
    hint = [next_on(t, 0) if desired[t] is None else None for t in range(6)]
    sim, program = run_op(geo, sim, desired, set(), None, hint, force_hint=True)
    program = list(program)
    i = 0
    limit = len(program) + K * (count + 2)
    while sim.accepted < count and len(program) < limit:
        i += 1
        line = ["-"] * 6
        for t in range(6):
            want = seq[t][i % K]
            p, on = sim.pos[t], sim.act[t]
            if want is not None:
                if not on:
                    if p != want:
                        return None
                    line[t] = "T"
                elif p != want:
                    if geo.dist(p, want) != 1:
                        return None
                    line[t] = geo.step(p, want)
            elif on:
                line[t] = "T"
            else:
                target = next_on(t, i % K)
                if target is not None and p != target:
                    line[t] = geo.step(p, target)
        sim.tick(line)
        program.append(line)
    return program if sim.accepted >= count else None


COMPRESS_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compress")
COMPRESS_MS = 8000


def compress_program(geo, T, S, count, program):
    """Let the compress helper delete ticks from a working program and repair it."""
    if not os.path.exists(COMPRESS_BIN):
        return None
    lines = [f"{geo.n} {len(program)} {count} {COMPRESS_MS} 7", " ".join(map(str, geo.start))]
    for t in range(6):
        for p in range(geo.n):
            lines.append(enc448(geo.bmask[t][p]))
    lines.append(enc448(geo.mask(T)))
    lines.append(enc448(geo.mask(S)))
    lines += [" ".join(l) for l in program]
    out = subprocess.run([COMPRESS_BIN], input="\n".join(lines), capture_output=True, text=True).stdout
    prog = [l.split() for l in out.strip().splitlines()]
    return prog or None


def sweep_plan(geo, T, S, a):
    """All lines of one axis in order; lines needing too many trims are left to the pen."""
    order, uncovered = [], []
    for v in range(-geo.R + 1, geo.R):
        cells = geo.line[(a, v)]
        tg = [x for x in cells if x in T]
        if not tg:
            continue
        need = [x for x in cells if x in S]
        if len(need) <= MAX_TRIMS and trim_options(geo, a, v, need, S, set()):
            order.append((a, v))
        else:
            uncovered.extend(tg)
    return order, uncovered


def solve_pattern(R, count, board_lines):
    geo = Geo(R)
    T = frozenset(board_lines_to_matter(board_lines, R))
    S = frozenset(x for t in T for x in neighbors(t) if max(map(abs, x)) < R and x not in T)
    best = build_program_for_pattern(R, set(T), count)
    if not verify(geo, T, S, best, count):
        best = None
    cyc = cycle_program(geo, T, S, count)
    if cyc is not None and verify(geo, T, S, cyc, count) and (best is None or len(cyc) < len(best)):
        best = cyc
    variants = [(p, 1) for p in (None, 0, 1, 2)] + [(None, k) for k in (3, 6, 10**9)]
    variants += [("sweep", a) for a in range(3)] + [("strict", a) for a in range(3)]
    for primary, min_new in variants:
        if primary == "sweep":
            order, uncovered = sweep_plan(geo, T, S, min_new)
        elif primary == "strict":
            order, uncovered = peel(geo, T, S, min_new, 1, strict=True)
        else:
            order, uncovered = peel(geo, T, S, primary, min_new)
        orders = [order]
        rev = list(reversed(order))
        if assign_trims(geo, rev, T, S) is not None:
            orders.append(rev)
        ops0 = assign_trims(geo, order, T, S)
        if ops0:
            opts = {k: o for k, o in ops0}
            nn = nn_order(geo, order, opts, geo.start)
            if assign_trims(geo, nn, T, S) is not None:
                orders.append(nn)
                nr = list(reversed(nn))
                if assign_trims(geo, nr, T, S) is not None:
                    orders.append(nr)
        for seed in (lambda k: (-k[1], 0), lambda k: (k[1], 0)):
            rn = reverse_nn_order(geo, order, T, S, seed)
            if rn is not None:
                orders.append(rn)
                if assign_trims(geo, rn[::-1], T, S) is not None:
                    orders.append(rn[::-1])
        seen = set()
        for order in orders:
            if tuple(order) in seen:
                continue
            seen.add(tuple(order))
            group = [order]
            if assign_trims(geo, order[::-1], T, S) is not None:
                group.append(order[::-1])
            for pen_at in (("start", "end") if uncovered else ("start",)):
                prog = compile_program(geo, T, S, group, uncovered, count, pen_at)
                if prog is None or (best is not None and len(prog) >= len(best)):
                    continue
                if verify(geo, T, S, prog, count):
                    best = prog
    if best is not None:
        kmax = min(CYCLEK_KMAX, (len(best) - 2) // count)
        cyc = cyclek_program(geo, T, S, count, kmax)
        if cyc is not None and len(cyc) < len(best) and verify(geo, T, S, cyc, count):
            best = cyc
    if best is not None:
        short = compress_program(geo, T, S, count, best)
        if short is not None and len(short) < len(best) and verify(geo, T, S, short, count):
            best = short
    return best


def solve(in_path, out_path):
    with open(in_path) as f:
        patterns = parse_patterns(f.read())
    texts = []
    total = 0
    with Pool() as pool:
        progs = pool.starmap(solve_pattern, patterns, chunksize=1)
    for i, prog in enumerate(progs):
        if prog is None:
            raise RuntimeError(f"pattern {i}: no valid program")
        total += len(prog)
        texts.append("\n".join(" ".join(l) for l in prog))
    with open(out_path, "w") as f:
        f.write("\n\n".join(texts) + "\n\n")
    return total


if __name__ == "__main__":
    total = solve(sys.argv[1], sys.argv[2])
    print("total instruction lines:", total, file=sys.stderr)
