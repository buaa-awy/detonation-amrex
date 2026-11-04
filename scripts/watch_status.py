#!/usr/bin/env python3
"""Print a compact, non-graphical live status line from status.json."""

import argparse
import json
import sys
import time
from pathlib import Path


def fmt_seconds(value: float) -> str:
    if value < 0.0:
        return "--:--:--"
    value = int(value + 0.5)
    return f"{value // 3600:02d}:{value % 3600 // 60:02d}:{value % 60:02d}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("status", nargs="?", default="status.json")
    parser.add_argument("--interval", type=float, default=1.0)
    args = parser.parse_args()
    path = Path(args.status)

    last_revision = None
    while True:
        try:
            revision = path.stat().st_mtime_ns
            if revision != last_revision:
                data = json.loads(path.read_text(encoding="utf-8"))
                width = 24
                filled = min(width, int(width * data["percent"] / 100.0))
                bar = "#" * filled + "." * (width - filled)
                levels = "/".join(str(x) for x in data["level_cells"])
                line = (
                    f"[{bar}] {data['percent']:6.2f}%"
                    f"  t={data['t_us']:9.3f}/{data['target_us']:.3f} us"
                    f"  step={data['step']:8d}"
                    f"  dt={data['dt_s']:.3e} s"
                    f"  L={data['finest_level']}"
                    f"  cells={data['cells']} ({levels})"
                    f"  front={data['front_x_m'] * 100:.3f} cm"
                    f"  react={data.get('reaction_front_x_m', 0.0) * 100:.3f} cm"
                    f"  gap={data.get('front_separation_mm', 0.0):.3f} mm"
                    f"  D={data.get('front_speed_m_s', 0.0):.0f} m/s"
                    f"  wall={fmt_seconds(data['wall_s'])}"
                    f"  ETA={fmt_seconds(data['eta_s'])}"
                )
                if sys.stdout.isatty():
                    print("\r" + line + " " * 8, end="", flush=True)
                else:
                    print(line, flush=True)
                last_revision = revision
                if data.get("finished"):
                    print()
                    return 0
        except (FileNotFoundError, json.JSONDecodeError, KeyError):
            pass
        time.sleep(args.interval)


if __name__ == "__main__":
    raise SystemExit(main())
