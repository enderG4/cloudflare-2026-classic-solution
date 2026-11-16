import sys


def hexagon_lines(r: int) -> list[str]:
    lines = []
    for i in range(2 * r + 1):
        k = min(i, 2 * r - i)
        indent = "-" * (r - k)
        if k == 0:
            tokens = ["*"] * (r + 1)
        else:
            tokens = ["*"] + ["O"] * (r + k - 1) + ["*"]
        lines.append(indent + "-".join(tokens) + indent)
    return lines


def solve(in_path: str, out_path: str) -> None:
    with open(in_path) as f:
        tokens = f.read().split()

    n = int(tokens[0])
    radii = [int(t) for t in tokens[1:1 + n]]

    boards = ["\n".join(hexagon_lines(r)) for r in radii]

    with open(out_path, "w") as f:
        f.write("\n\n".join(boards) + "\n\n")


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    solve(in_path, out_path)
