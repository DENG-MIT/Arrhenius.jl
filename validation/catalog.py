"""Freeze the linked Cantera Python example inventory without executing examples."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen

INDEX = "https://cantera.org/dev/examples/python/index.html"


class ExampleLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.main = False
        self.links = set()

    def handle_starttag(self, tag, attrs):
        if tag == "main":
            self.main = True
        if self.main and tag == "a":
            href = dict(attrs).get("href", "")
            url = urljoin(INDEX, href).split("#")[0]
            path = urlparse(url).path
            if (url.startswith(INDEX.rsplit("/", 1)[0] + "/")
                    and path.endswith(".html") and not path.endswith("/index.html")):
                self.links.add(url)

    def handle_endtag(self, tag):
        if tag == "main":
            self.main = False


def snapshot(html: bytes):
    parser = ExampleLinks()
    parser.feed(html.decode("utf-8"))
    if not parser.links:
        raise ValueError("no Python examples found in the main document")
    return {
        "index_url": INDEX,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "index_sha256": hashlib.sha256(html).hexdigest(),
        "minimum_speed_ratio": 0.95,
        "speed_ratio_definition": "Cantera elapsed / Arrhenius elapsed",
        "repeated_run_timing_contract": "Warm repeated calculations reuse the loaded mechanism in both Arrhenius and Cantera. Each reactor repetition starts from fresh state and includes reactor and solver setup, calculation, and required numerical output extraction. Imports, compilation, one-time mechanism/sidecar loading, file I/O, and validation checks are separate. Historical per-case evidence retains its recorded boundary and is not relabeled.",
        "performance_target": "WSL",
        "performance_targets": ["wsl"],
        "examples": [
            {"id": url.split("/examples/python/")[1][:-5], "url": url,
             "status": "not_validated", "julia_example": None,
             "correctness_evidence": None, "wsl_timing_evidence": None,
             "m4_timing_evidence": None,
             "performance_status": "not_validated"}
            for url in sorted(parser.links)
        ],
    }


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("output", type=Path)
    opts = args.parse_args()
    with urlopen(INDEX, timeout=30) as response:
        data = snapshot(response.read())
    opts.output.parent.mkdir(parents=True, exist_ok=True)
    opts.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"{len(data['examples'])} examples: {opts.output}")
