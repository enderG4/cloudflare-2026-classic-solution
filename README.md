# Cloudflare 2026 Classic — Solutions

Solutions to the **Cloudflare 2026 "Nanofabricator"** hackathon problems (levels 1–6).

The puzzle: a hexagonal board with six emitters that move around its outer ring. Each
tick every emitter can move forward/back (`F`/`B`), toggle its beam (`T`), or do nothing
(`-`). Beams cross the board and paint only tiles marked `O`; you must write the tick
program that reproduces the target pattern `M`.

## Contents

| File | Level | Notes |
| --- | --- | --- |
| `solve.py` | 1 | Draws plain hexagonal boards |
| `solve4.py` | 4 | Shared cube-coordinate geometry + beam simulator |
| `solve6.py` | 6 | Greedy beam-by-beam program builder |
| `solve6b.py` | 6 | Stronger search: beam planning, pruning, C++ compressor assist |
| `cycle2.cpp`, `cyclek.cpp` | – | Ring/cycle search helpers |
| `pathk.cpp` | – | Path search over ring positions |
| `compress.cpp` | – | Shrinks a working program by delete-and-repair local search |
| `run.sh` | 6 | Runs `solve6b.py` over all level-6 inputs, prints timings and line-count score |
| `sandbox-level-6.html` | – | Official interactive visualizer/debugger |
| `tests/in_level-*`, `tests/out_level-*` | – | Test inputs and generated outputs |

The PDF problem statements and the original `.zip` bundles are intentionally not
committed (see `.gitignore`).

## Usage

Level 1 and level 6 solvers take an input and output path:

```bash
python3 solve.py    tests/in_level-1_0-example.txt tests/out_level-1_0-example.txt
python3 solve6b.py  tests/in_level-6_2-medium.txt tests/out_level-6_2-medium.txt
```

Run the full level-6 suite with scoring:

```bash
./run.sh
```

## Building the C++ helpers

```bash
cmake -S . -B build
cmake --build build -j
```

Binaries are written to `build/bin/`; `solve6b.py` finds them automatically
(it also still accepts prebuilt copies in the repo root or `bin/`).

## Debugging

Open `sandbox-level-6.html` in a browser and paste in the input/output pair to replay
the emitters and see the beams.
