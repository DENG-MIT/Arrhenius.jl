#!/usr/bin/env python3
"""Prepare native/Cantera YAML and LoKI classic input from one LXCat TXT set.

Only text formatting changes. Every input energy/cross-section token is retained
verbatim in both generated tables; source and table hashes go to source-audit.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

KINDS = {"EFFECTIVE": "effective", "ELASTIC": "elastic", "EXCITATION": "excitation",
         "IONIZATION": "ionization", "ATTACHMENT": "attachment"}


def parse(source: Path):
    lines = source.read_text(encoding="utf-8-sig").splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        kind = lines[i].strip()
        if kind not in KINDS or i + 2 >= len(lines):
            i += 1
            continue
        target_line = lines[i + 1].strip()
        target = target_line.split("->", 1)[0].strip()
        third = lines[i + 2].strip()
        if target != "Ar":
            i += 1
            continue
        j = i + 3
        while j < len(lines) and not lines[j].startswith("COLUMNS:"):
            j += 1
        if j == len(lines):
            raise ValueError(f"{kind}: missing COLUMNS marker")
        j += 1
        while j < len(lines) and not lines[j].lstrip().startswith("---"):
            j += 1
        j += 1
        table = []
        while j < len(lines) and not lines[j].lstrip().startswith("---"):
            parts = lines[j].split()
            if len(parts) >= 2:
                try:
                    float(parts[0]); float(parts[1])
                except ValueError:
                    pass
                else:
                    table.append((parts[0], parts[1]))
            j += 1
        if len(table) < 2:
            raise ValueError(f"{kind}: missing or short cross-section table")
        threshold = 0.0 if kind in ("EFFECTIVE", "ELASTIC") else float(third)
        blocks.append({"raw_kind": kind, "kind": KINDS[kind], "target": target,
                       "third": third, "threshold": threshold, "table": table})
        i = j + 1
    if not blocks:
        raise ValueError("No Ar collisions found")
    expected = ["effective", "excitation", "ionization"]
    if [b["kind"] for b in blocks] != expected:
        raise ValueError(f"Expected collision order {expected}, got {[b['kind'] for b in blocks]}")
    return blocks


def yaml_sequence(key, values, indent=2, per_line=8):
    pad = " " * indent
    lines = [pad + key + ": ["]
    for start in range(0, len(values), per_line):
        chunk = values[start:start+per_line]
        suffix = "," if start + per_line < len(values) else ""
        lines.append(pad + "  " + ", ".join(chunk) + suffix)
    lines.append(pad + "]")
    return lines


def make_yaml(blocks, cells: int, emax: float) -> str:
    grid = [f"{emax*k/cells:.17g}" for k in range(cells + 1)]
    out = [
        "description: LXCat Phelps Ar complete set; cross-section tables preserved from source.",
        "phases:",
        "- name: ar-plasma",
        "  thermo: plasma",
        "  kinetics: gas",
        "  elements: [Ar, E]",
        "  species: [Electron, Ar]",
        "  state: {T: 300.0, P: 133.32, X: {Ar: 1.0, Electron: 1.0e-30}}",
        "  electron-energy-distribution:",
        "    type: Boltzmann-two-term",
    ]
    out.extend(yaml_sequence("energy-levels", grid, indent=4, per_line=8))
    out.append("  reactions: []")
    out.extend([
        "species:",
        "- name: Electron",
        "  composition: {E: 1}",
        "  thermo:",
        "    model: NASA7",
        "    temperature-ranges: [200.0, 1000.0, 6000.0]",
        "    data:",
        "    - [2.5, 0, 0, 0, 0, -745.375, -11.7246902]",
        "    - [2.5, 0, 0, 0, 0, -745.375, -11.7246902]",
        "- name: Ar",
        "  composition: {Ar: 1}",
        "  thermo:",
        "    model: NASA7",
        "    temperature-ranges: [200.0, 1000.0, 6000.0]",
        "    data:",
        "    - [2.5, 0, 0, 0, 0, 0, 0]",
        "    - [2.5, 0, 0, 0, 0, 0, 0]",
        "electron-collisions:",
    ])
    for b in blocks:
        out.extend([f"- target: {b['target']}", f"  kind: {b['kind']}",
                    f"  threshold: {b['threshold']:.17g}"])
        out.extend(yaml_sequence("energy-levels", [x for x, _ in b["table"]], indent=2))
        out.extend(yaml_sequence("cross-sections", [y for _, y in b["table"]], indent=2))
    return "\n".join(out) + "\n"

def loki_classic(blocks) -> str:
    out = ["% Derived from the accompanying LXCat Phelps Ar download.",
           "% Table knot values are copied without numerical changes."]
    names = {"effective": "Ar(1S0) + e -> Ar(1S0) + e, Effective",
             "excitation": "Ar(1S0) + e -> Ar(11.5eV) + e, Excitation",
             "ionization": "Ar(1S0) + e -> Ar(+,gnd) + e + e, Ionization"}
    for b in blocks:
        if b["kind"] == "effective":
            param = f"m/M = {b['third'].split()[0]}, complete set"
        else:
            param = f"E = {b['third']} eV, complete set"
        out.extend([f"PARAM.: {param}", f"[{names[b['kind']]}]",
                    "COLUMNS: Energy (eV) | Cross section (m2)", "-----------------------------"])
        out.extend(f"{e}\t{s}" for e, s in b["table"])
        out.append("-----------------------------")
        out.append("")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--cells", type=int, default=2000)
    ap.add_argument("--emax", type=float, default=100.0)
    args = ap.parse_args()
    if args.cells < 2 or args.emax <= 0:
        raise SystemExit("cells must be >=2 and emax >0")
    raw = args.source.read_bytes()
    blocks = parse(args.source)
    args.output.mkdir(parents=True, exist_ok=True)
    yaml_path = args.output / "Phelps_Ar_probe.yaml"
    loki_path = args.output / "Phelps_Ar_LoKI.txt"
    yaml_path.write_text(make_yaml(blocks, args.cells, args.emax), encoding="utf-8")
    loki_path.write_text(loki_classic(blocks), encoding="utf-8")
    tables = []
    for b in blocks:
        canonical = "\n".join(f"{e}\t{s}" for e, s in b["table"]).encode()
        tables.append({"kind": b["kind"], "threshold_eV": b["threshold"],
                       "knots": len(b["table"]), "table_sha256": hashlib.sha256(canonical).hexdigest()})
    audit = {"source": str(args.source), "source_sha256": hashlib.sha256(raw).hexdigest(),
             "source_bytes": len(raw), "retrieved": "2026-10-17", "database": "Phelps database",
             "species": "Ar", "data_format": "LXCat plain text", "citation":
             "Phelps database, www.lxcat.net, retrieved on October 17, 2026; see source header for contributor references.",
             "grid": {"min_eV": 0.0, "max_eV": args.emax, "cells": args.cells,
                      "points": args.cells + 1}, "tables": tables,
             "generated": {yaml_path.name: hashlib.sha256(yaml_path.read_bytes()).hexdigest(),
                           loki_path.name: hashlib.sha256(loki_path.read_bytes()).hexdigest()}}
    table_path = args.output / "probe-tables.json"
    table_payload = [{"kind": b["kind"], "threshold_eV": b["threshold"],
                      "energy_eV": [float(e) for e, _ in b["table"]],
                      "cross_section_m2": [float(s) for _, s in b["table"]]} for b in blocks]
    table_path.write_text(json.dumps(table_payload, separators=(",", ":")) + "\n", encoding="utf-8")
    audit["generated"][table_path.name] = hashlib.sha256(table_path.read_bytes()).hexdigest()
    (args.output / "source-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()








