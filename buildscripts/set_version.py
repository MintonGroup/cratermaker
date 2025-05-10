import re
import subprocess
from pathlib import Path

version = (
    subprocess.check_output(["git", "describe", "--tags"])
    .decode()
    .replace("v", "")
    .replace(".0", ".")
    .strip()
)
# Because Rust and Python have different versioning schemes, we need to scrub the raw version string to make it compatible with both.
# Example (develop): git tag: "v2025.05.3-alpha-41-g7a2bba8" => version "2025.5.3-a41+g7a2bba8"
# Example (release): git tag: "v2025.05.3-alpha" => version "2025.5.3-a0"
# First we try to split it along "alpha-", which is what the version has if it is a pre-release development version. .
version = version.split("alpha-")
if len(version) > 1:
    version[1] = version[1].replace(
        "-", "+"
    )  # replace "-" with "+" to make it compatible with Python/Rust versioning
    version = "a".join(version)
else:  # This is a release version, so there is no dash after alpha, but it needs to be changed to a0
    version = version[0].replace("alpha", "a0")

root_path = Path(__file__).resolve().parents[1]
cargo_file = root_path / "Cargo.toml"
version_file = root_path / "cratermaker" / "_version.py"

# Read Cargo.toml as text
with open(cargo_file, "r", encoding="utf-8") as f:
    cargo_contents = f.read()

# Replace the version line
new_cargo_contents = re.sub(
    r'version\s*=\s*".*?"',
    f'version = "{version}"',
    cargo_contents,
    count=1,
)

# Write back
with open(cargo_file, "w", encoding="utf-8") as f:
    f.write(new_cargo_contents)
with open(version_file, "w") as f:
    f.write(f'__version__ = version = "{version}"\n')

print(f"Updated Cargo.toml to version: {version}")
