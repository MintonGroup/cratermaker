# set_version.py
import toml
import subprocess
import re

raw_version = subprocess.check_output(["python", "-m", "setuptools_scm"]).decode().strip()

match = re.match(r'^(\d+\.\d+\.\d+)([^\d].*)$', raw_version)
if match:
    version = f"{match.group(1)}-{match.group(2)}"
else:
    version = raw_version

cargo_path = "Cargo.toml"
data = toml.load(cargo_path)
data["package"]["version"] = version

with open(cargo_path, "w") as f:
    toml.dump(data, f)

print(f"Updated Cargo.toml to version: {version}")