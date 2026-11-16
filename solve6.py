import sys

from solve4 import D, build_ring, add, is_interior, simulate


def tail_position(direction_num, fixed_idx, fixed_val, R):
    # cube index layout: 0=q, 1=r_ax, 2=s
    dirv = D[direction_num]
    inc_idx = next(i for i in range(3) if dirv[i] == 1)
    dec_idx = next(i for i in range(3) if dirv[i] == -1)
    assert dirv[fixed_idx] == 0

    inc_min = max(-R, -R - fixed_val)
    cube = [0, 0, 0]
    cube[fixed_idx] = fixed_val
    cube[inc_idx] = inc_min
    cube[dec_idx] = -(fixed_val + inc_min)
    return tuple(cube)


def ring_shortest_steps(cur, target, n):
    fwd = (target - cur) % n
    bwd = (cur - target) % n
    if fwd <= bwd:
        return fwd, "F"
    else:
        return bwd, "B"


class ProgramBuilder:
    def __init__(self, R):
        self.R = R
        self.n = 6 * R
        self.ring = build_ring(R)
        self.pos = [(k - 1) * R for k in range(1, 7)]
        self.active = [False] * 6
        self.lines = []
        self.A_activated = False

    def emit(self, instr_map):
        line = ["-"] * 6
        for idx, instr in instr_map.items():
            line[idx] = instr
        self.lines.append(line)

    def move_one(self, idx, target_pos):
        while self.pos[idx] != target_pos:
            steps, direction = ring_shortest_steps(self.pos[idx], target_pos, self.n)
            self.emit({idx: direction})
            if direction == "F":
                self.pos[idx] = (self.pos[idx] + 1) % self.n
            else:
                self.pos[idx] = (self.pos[idx] - 1) % self.n

    def toggle(self, idx):
        self.emit({idx: "T"})
        self.active[idx] = not self.active[idx]


# axis 0 = q, axis 1 = r_ax, axis 2 = s; canonical direction fixing each axis
CANON_DIR_FOR_AXIS = {0: 2, 1: 1, 2: 3}


def build_groups_sequence(by_anchor, count, anchor_desc_start, sweep_desc_start):
    """
    One full pass over all (anchor_val -> sweep_vals) groups achieves the pattern
    once. Since a match immediately deletes the matched cells (see sandbox
    semantics), each of the `count` required occurrences needs its own full pass.
    Consecutive passes alternate traversal direction (ping-pong) so the anchor
    emitter's last position in pass k is adjacent to its first position in
    pass k+1, instead of paying a full traverse-back cost every repeat.
    """
    anchors_asc = sorted(by_anchor.keys())
    groups = []
    for rep in range(count):
        ad = anchor_desc_start if rep % 2 == 0 else not anchor_desc_start
        sd = sweep_desc_start if rep % 2 == 0 else not sweep_desc_start
        a_vals = list(reversed(anchors_asc)) if ad else anchors_asc
        step = -1 if sd else 1
        for a_val in a_vals:
            sweep_vals = sorted(by_anchor[a_val], reverse=sd)
            groups.append((a_val, sweep_vals, step))
    return groups


def build_program_variant(R, target_cells, count, anchor_axis, sweep_axis,
                           anchor_desc_start=False, sweep_desc_start=False):
    anchor_dir = CANON_DIR_FOR_AXIS[anchor_axis]
    sweep_dir = CANON_DIR_FOR_AXIS[sweep_axis]

    builder = ProgramBuilder(R)
    A = anchor_dir - 1
    B = sweep_dir - 1

    if not target_cells:
        return builder.lines

    by_anchor = {}
    for cell in target_cells:
        by_anchor.setdefault(cell[anchor_axis], []).append(cell[sweep_axis])

    groups = build_groups_sequence(by_anchor, count, anchor_desc_start, sweep_desc_start)

    for a_val, sweep_vals, step in groups:
        entryA = tail_position(anchor_dir, anchor_axis, a_val, R)
        entryA_idx = builder.ring.index(entryA)

        if builder.active[B]:
            builder.toggle(B)

        builder.move_one(A, entryA_idx)

        if not builder.A_activated:
            builder.toggle(A)
            builder.A_activated = True

        prev_sv = None
        for sv in sweep_vals:
            entryB = tail_position(sweep_dir, sweep_axis, sv, R)
            entryB_idx = builder.ring.index(entryB)

            if prev_sv is not None and sv == prev_sv + step and builder.active[B]:
                steps, direction = ring_shortest_steps(builder.pos[B], entryB_idx, builder.n)
                if steps == 1:
                    builder.emit({B: direction})
                    if direction == "F":
                        builder.pos[B] = (builder.pos[B] + 1) % builder.n
                    else:
                        builder.pos[B] = (builder.pos[B] - 1) % builder.n
                    prev_sv = sv
                    continue

            if builder.active[B]:
                builder.toggle(B)
            builder.move_one(B, entryB_idx)
            builder.toggle(B)
            prev_sv = sv

    return builder.lines


AXIS_PAIRS = [
    (0, 1), (0, 2),
    (1, 0), (1, 2),
    (2, 0), (2, 1),
]


def build_program_for_pattern(R, target_cells, count):
    best = None
    for anchor_axis, sweep_axis in AXIS_PAIRS:
        for anchor_desc in (False, True):
            for sweep_desc in (False, True):
                candidate = build_program_variant(
                    R, target_cells, count, anchor_axis, sweep_axis,
                    anchor_desc_start=anchor_desc, sweep_desc_start=sweep_desc,
                )
                if best is None or len(candidate) < len(best):
                    best = candidate
    return best


def format_instr_line(line):
    return " ".join(line)


def parse_patterns(text):
    lines = text.splitlines()
    pos = 0

    def next_nonblank():
        nonlocal pos
        while pos < len(lines) and lines[pos].strip() == "":
            pos += 1
        line = lines[pos]
        pos += 1
        return line

    n = int(next_nonblank())
    patterns = []
    for _ in range(n):
        header = next_nonblank().split()
        R, rows, cols, count = int(header[0]), int(header[1]), int(header[2]), int(header[3])
        board_lines = [next_nonblank() for _ in range(rows)]
        patterns.append((R, count, board_lines))
    return patterns


def board_lines_to_matter(board_lines, R):
    matter = set()
    for i, line in enumerate(board_lines):
        k = min(i, 2 * R - i)
        indent = R - k
        inner = line[indent: len(line) - indent] if indent > 0 else line
        tokens = inner.split("-")
        r_ax = i - R
        for j, tok in enumerate(tokens):
            if tok == "M":
                if r_ax <= 0:
                    q = j - i
                else:
                    q = j - R
                s = -q - r_ax
                matter.add((q, r_ax, s))
    return matter


def count_submissions(R, program, target):
    """
    Simulate program tick by tick. Each tick: apply move/toggle instructions,
    recompute matter via the create(==2)/destroy(>=3) rule, then check whether
    the target pattern is fully present. A match counts as one accepted
    submission and immediately clears the matched cells (sandbox semantics:
    matching deletes the matter), so a subsequent tick must recreate them to
    register again.
    """
    n = 6 * R
    ring = build_ring(R)
    pos = [(k - 1) * R for k in range(1, 7)]
    active = [False] * 6
    matter = set()
    submissions = 0
    for line in program:
        for idx in range(6):
            instr = line[idx]
            if instr == "F":
                pos[idx] = (pos[idx] + 1) % n
            elif instr == "B":
                pos[idx] = (pos[idx] - 1) % n
            elif instr == "T":
                active[idx] = not active[idx]
        counts = {}
        for idx in range(6):
            if not active[idx]:
                continue
            direction = D[idx + 1]
            cell = add(ring[pos[idx]], direction)
            while is_interior(cell, R):
                counts[cell] = counts.get(cell, 0) + 1
                cell = add(cell, direction)
        for cell, c in counts.items():
            if c == 2:
                matter.add(cell)
            elif c >= 3:
                matter.discard(cell)
        if target.issubset(matter):
            submissions += 1
            matter -= target
    return submissions, matter


def solve(in_path, out_path, verify=True):
    with open(in_path) as f:
        patterns = parse_patterns(f.read())

    programs_text = []
    total_lines = 0
    for pattern_idx, (R, count, board_lines) in enumerate(patterns):
        target = board_lines_to_matter(board_lines, R)
        full_program = build_program_for_pattern(R, target, count)

        if verify:
            submissions, _ = count_submissions(R, full_program, target)
            if submissions != count:
                raise RuntimeError(
                    f"pattern {pattern_idx}: expected {count} submissions, got {submissions}"
                )

        total_lines += len(full_program)
        programs_text.append("\n".join(format_instr_line(l) for l in full_program))

    with open(out_path, "w") as f:
        f.write("\n\n".join(programs_text) + "\n\n")

    return total_lines


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    total = solve(in_path, out_path)
    print("total instruction lines:", total, file=sys.stderr)
