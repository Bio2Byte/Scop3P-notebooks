#!/usr/bin/env python3
"""Download the pinned browser libraries into a vendor directory.

Run at image build time so the running app never depends on a public CDN. The manifest
lives in ``apps/common/vendor.py`` -- one place, so a version bump cannot leave the
download and the page disagreeing.

    python scripts/fetch-vendor-assets.py [--dest DIR] [--print-hashes]

Note the trade this makes: the *build* now depends on the CDNs, even though the *runtime*
no longer does. That is the right way round -- a build failure is visible to whoever is
building, while a runtime failure lands on a user in the middle of an experiment.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_manifest(path: Path):
    """Load the manifest module directly from its file.

    Deliberately not ``from common.vendor import ...``: that would execute
    ``apps/common/__init__.py``, which imports pandas and the rest of the app's dependency
    graph. This script runs during an image build, where the point is to fetch static files
    -- it should not need the science stack to be importable, and should not fail if it is
    not yet.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("_scop3p_vendor_manifest", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load the vendor manifest from {path}")
    module = importlib.util.module_from_spec(spec)
    # Registered before executing: the manifest uses dataclasses under
    # ``from __future__ import annotations``, and dataclasses resolves those string
    # annotations through ``sys.modules[cls.__module__]``. A module that is not registered
    # fails there with a confusing AttributeError inside dataclasses itself.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_MANIFEST_PATH = Path(
    os.getenv("SCOP3P_VENDOR_MANIFEST", REPO_ROOT / "apps" / "common" / "vendor.py")
)
_manifest = _load_manifest(_MANIFEST_PATH)
ASSETS = _manifest.ASSETS
DEFAULT_VENDOR_DIR = _manifest.DEFAULT_VENDOR_DIR
VendoredAsset = _manifest.VendoredAsset

TIMEOUT_SECONDS = 60


def fetch(asset: VendoredAsset, destination: Path, *, force: bool = False) -> tuple[bool, str]:
    """Download one asset. Returns (downloaded, sha256)."""
    target = destination / asset.filename
    if target.is_file() and target.stat().st_size > 0 and not force:
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        return False, digest

    request = urllib.request.Request(
        asset.url,
        headers={"User-Agent": "scop3p-toolkit-build"},
    )
    with urllib.request.urlopen(request, timeout=TIMEOUT_SECONDS) as response:
        payload = response.read()

    if not payload:
        raise RuntimeError(f"{asset.url} returned an empty body")

    digest = hashlib.sha256(payload).hexdigest()
    if asset.sha256 and digest != asset.sha256:
        # A pinned version whose bytes changed. Refusing is the point: silently shipping
        # different code than was tested is exactly what pinning is meant to prevent.
        raise RuntimeError(
            f"{asset.filename}: expected sha256 {asset.sha256}, got {digest}. "
            "The CDN served different bytes for a pinned version; verify before updating "
            "the manifest."
        )

    destination.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    return True, digest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, default=DEFAULT_VENDOR_DIR)
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    parser.add_argument(
        "--print-hashes",
        action="store_true",
        help="print a manifest snippet with the observed hashes, to pin them",
    )
    args = parser.parse_args()

    failures: list[str] = []
    hashes: dict[str, str] = {}
    for asset in ASSETS:
        try:
            downloaded, digest = fetch(asset, args.dest, force=args.force)
        except Exception as error:  # noqa: BLE001 - report every failure, not just the first
            print(f"  FAILED  {asset.filename}: {error}", file=sys.stderr)
            failures.append(asset.filename)
            continue
        hashes[asset.key] = digest
        size_kb = (args.dest / asset.filename).stat().st_size / 1024
        state = "downloaded" if downloaded else "cached"
        print(f"  {state:10s} {asset.filename:24s} {size_kb:8.1f} KB  {digest[:12]}")

    if args.print_hashes:
        print("\n# sha256 values observed, for pinning in apps/common/vendor.py:")
        for key, digest in hashes.items():
            print(f'#   {key}: "{digest}"')

    if failures:
        print(f"\n{len(failures)} asset(s) could not be fetched: {failures}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
