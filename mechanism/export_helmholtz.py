"""Export full pure-fluid EOS parameters for native Julia Clapeyron.SingleFluid.

Example: python mechanism/export_helmholtz.py CO2 carbon-dioxide.json
Requires CoolProp for parameter preprocessing only.
"""
import argparse
import json
from pathlib import Path


def export(fluid: str, destination: Path) -> None:
    from CoolProp.CoolProp import get_fluid_param_string

    records = json.loads(get_fluid_param_string(fluid, "JSON"))
    if len(records) != 1 or not records[0].get("EOS"):
        raise ValueError("expected one complete pure-fluid EOS record")
    destination.write_text(json.dumps(records[0], indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fluid")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    export(args.fluid, args.destination)
