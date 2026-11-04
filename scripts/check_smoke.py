#!/usr/bin/env python3
"""Validate a completed deterministic smoke-test status file."""

import argparse
import json
import math
from pathlib import Path


EXPECTED_STEP = 100
EXPECTED_LEVEL = 3
EXPECTED_MASS = 1260.95732991
EXPECTED_ENERGY = 193540164.2116359


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("status", type=Path)
    args = parser.parse_args()
    data = json.loads(args.status.read_text(encoding="utf-8"))
    errors = []
    if data.get("step") != EXPECTED_STEP:
        errors.append("step=%s (expected %d)" % (
            data.get("step"), EXPECTED_STEP))
    if data.get("finest_level") != EXPECTED_LEVEL:
        errors.append("finest_level=%s (expected %d)" % (
            data.get("finest_level"), EXPECTED_LEVEL))
    if not data.get("finished"):
        errors.append("run is not marked finished")
    mass_error = abs(data["mass"] / EXPECTED_MASS - 1.0)
    energy_error = abs(data["energy"] / EXPECTED_ENERGY - 1.0)
    if mass_error > 1.0e-10:
        errors.append("mass differs by %.3e" % mass_error)
    if energy_error > 1.0e-9:
        errors.append("energy differs by %.3e" % energy_error)
    if not all(math.isfinite(float(data[key]))
               for key in ("dt_s", "mass", "energy")):
        errors.append("non-finite final state")
    if errors:
        raise SystemExit("Smoke test FAILED: " + "; ".join(errors))
    print(
        "Smoke test PASSED: step=%d t=%.6f us finest=%d mass=%.8f energy=%.3f"
        % (
            data["step"], data["t_us"], data["finest_level"],
            data["mass"], data["energy"]))


if __name__ == "__main__":
    main()
