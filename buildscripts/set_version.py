import subprocess
import re

# Get version string from setuptools_scm
raw_version = subprocess.check_output(["python", "-m", "setuptools_scm"]).decode().strip()

match = re.match(r'^(\d+\.\d+\.\d+)([^\d].*)$', raw_version)
if match:
    version = f"{match.group(1)}-{match.group(2)}"
else:
    version = raw_version

cargo_path = "Cargo.toml"

# Read Cargo.toml as text
with open(cargo_path, "r", encoding="utf-8") as f:
    cargo_contents = f.read()

# Replace the version line
new_cargo_contents = re.sub(
    r'version\s*=\s*".*?"',
    f'version = "{version}"',
    cargo_contents,
    count=1,
)

# Write back
with open(cargo_path, "w", encoding="utf-8") as f:
    f.write(new_cargo_contents)

print(f"Updated Cargo.toml to version: {version}")