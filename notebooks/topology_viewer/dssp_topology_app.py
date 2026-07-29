"""Backwards-compatible entry point.

The original module carried the whole application in one file. It now lives in
the ``topology`` package; this shim keeps ``from dssp_topology_app import
make_app`` working so existing notebooks do not need editing.
"""

from topology import build_view, load_structure, make_app, render

__all__ = ["make_app", "build_view", "load_structure", "render"]

if __name__ == "__main__":
    from IPython.display import display
    display(make_app())
