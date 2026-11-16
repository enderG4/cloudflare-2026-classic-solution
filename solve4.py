import sys

# Cube directions (q, r, s) with q + r + s == 0, 1-indexed 1..6
D = {
    1: (1, 0, -1),
    2: (0, 1, -1),
    3: (-1, 1, 0),
    4: (-1, 0, 1),
    5: (0, -1, 1),
    6: (1, -1, 0),
}


def scale(v, k):
    return (v[0] * k, v[1] * k, v[2] * k)


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def build_ring(r):
    # corner k (1..6) = -r * D[k]
    corners = {k: scale(D[k], -r) for k in range(1, 7)}
    # edge from corner k to corner k+1 uses direction D[k-1] (1-indexed, wraps 1->6)
    edge_dir = {k: D[k - 1 if k > 1 else 6] for k in range(1, 7)}
    ring = []
    for k in range(1, 7):
        start = corners[k]
        step = edge_dir[k]
        for t in range(r):
            ring.append(add(start, scale(step, t)))
    return ring


def is_interior(cell, r):
    q, rr, s = cell
    return max(abs(q), abs(rr), abs(s)) < r


def hexagon_board_lines(r, matter):
    lines = []
    for i in range(2 * r + 1):
        k = min(i, 2 * r - i)
        r_ax = i - r
        indent = "-" * (r - k)
        if k == 0:
            tokens = ["*"] * (r + 1)
        else:
            tokens = ["*"]
            count = r + k - 1
            for c in range(count):
                j = c + 1
                if r_ax <= 0:
                    q = j - i
                else:
                    q = j - r
                s = -q - r_ax
                tokens.append("M" if (q, r_ax, s) in matter else "O")
            tokens.append("*")
        lines.append(indent + "-".join(tokens) + indent)
    return lines


def simulate(r, program):
    n = 6 * r
    ring = build_ring(r)
    pos = [(k - 1) * r for k in range(1, 7)]  # ring index for emitter 1..6
    active = [False] * 6

    matter = set()

    for line in program:
        # step 1: apply move/toggle instructions simultaneously
        for idx in range(6):
            instr = line[idx]
            if instr == "F":
                pos[idx] = (pos[idx] + 1) % n
            elif instr == "B":
                pos[idx] = (pos[idx] - 1) % n
            elif instr == "T":
                active[idx] = not active[idx]
            # '-' no-op

        # step 2: count active beams per cell, then create (==2) / destroy (>=3)
        count = {}
        for idx in range(6):
            if not active[idx]:
                continue
            emitter_num = idx + 1
            start = ring[pos[idx]]
            direction = D[emitter_num]
            cell = add(start, direction)
            while is_interior(cell, r):
                count[cell] = count.get(cell, 0) + 1
                cell = add(cell, direction)

        for cell, cnt in count.items():
            if cnt == 2:
                matter.add(cell)
            elif cnt >= 3:
                matter.discard(cell)

    return matter


def parse_programs(text):
    lines = [ln for ln in text.splitlines()]
    pos = 0

    def next_nonblank():
        nonlocal pos
        while pos < len(lines) and lines[pos].strip() == "":
            pos += 1
        line = lines[pos]
        pos += 1
        return line

    n = int(next_nonblank())
    programs = []
    for _ in range(n):
        r = int(next_nonblank())
        l = int(next_nonblank())
        prog = []
        for _ in range(l):
            prog.append(next_nonblank().split())
        programs.append((r, prog))
    return programs


def solve(in_path, out_path):
    with open(in_path) as f:
        text = f.read()

    programs = parse_programs(text)

    boards = []
    for r, program in programs:
        matter = simulate(r, program)
        boards.append("\n".join(hexagon_board_lines(r, matter)))

    with open(out_path, "w") as f:
        f.write("\n\n".join(boards) + "\n\n")


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    solve(in_path, out_path)
