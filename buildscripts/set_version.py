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

version = version.split("alpha-")
if len(version) > 1:
    version[1] = version[1].replace("-", "+")
    version = "a".join(version)
else:
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
