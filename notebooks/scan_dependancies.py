#!/usr/bin/env python3

import ast
import json
import sys
import importlib.metadata
from pathlib import Path

# Map common import names -> PyPI package names
IMPORT_TO_PACKAGE = {
    "Bio": "biopython",
    "PIL": "Pillow",
    "cv2": "opencv-python",
    "sklearn": "scikit-learn",
    "yaml": "PyYAML",
    "bs4": "beautifulsoup4",
    "skimage": "scikit-image",
    "mpl_toolkits": "matplotlib",
}

# Python standard-library modules
STDLIB = getattr(sys, "stdlib_module_names", set())

# Collect imports and where they occur
imports = {}

notebooks = list(Path(".").rglob("*.ipynb"))

print(f"\nFound {len(notebooks)} notebooks.\n")

for nb_path in notebooks:

    try:
        with open(nb_path, encoding="utf-8") as f:
            nb = json.load(f)
    except Exception as e:
        print(f"Could not read {nb_path}: {e}")
        continue

    for cell in nb.get("cells", []):

        if cell.get("cell_type") != "code":
            continue

        source = "".join(cell.get("source", []))

        # Remove notebook magics / shell commands
        source = "\n".join(
            line for line in source.splitlines()
            if not line.lstrip().startswith(("%", "!"))
        )

        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue

        for node in ast.walk(tree):

            if isinstance(node, ast.Import):
                for alias in node.names:
                    module = alias.name.split(".")[0]
                    imports.setdefault(module, set()).add(str(nb_path))

            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    module = node.module.split(".")[0]
                    imports.setdefault(module, set()).add(str(nb_path))


print("=" * 75)
print("THIRD-PARTY IMPORTS FOUND")
print("=" * 75)

requirements = []

for module in sorted(imports):

    # Skip Python standard library
    if module in STDLIB:
        continue

    package = IMPORT_TO_PACKAGE.get(module, module)

    try:
        version = importlib.metadata.version(package)
        status = f"{package}=={version}"
        requirements.append(status)

    except importlib.metadata.PackageNotFoundError:
        version = "NOT FOUND"
        status = f"{package}   [version not found]"

    print(f"\n{module}")
    print(f"  Package : {status}")
    print("  Used in :")

    for nb in sorted(imports[module]):
        print(f"            {nb}")


print("\n")
print("=" * 75)
print("POTENTIAL REQUIREMENTS")
print("=" * 75)

for req in sorted(set(requirements), key=str.lower):
    print(req)


# Write result to a file
output = Path("detected_requirements.txt")

with open(output, "w") as f:
    for req in sorted(set(requirements), key=str.lower):
        f.write(req + "\n")

print(f"\nSaved to: {output}")