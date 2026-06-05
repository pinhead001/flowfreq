"""Build script for f2py compilation of peakfqr Fortran sources."""
import os
import sys
import subprocess

# Add required paths to environment before invoking f2py
scripts_dir = os.path.join(
    os.environ["LOCALAPPDATA"],
    r"Packages\PythonSoftwareFoundation.Python.3.12_qbz5n2kfra8p0"
    r"\LocalCache\local-packages\Python312\Scripts",
)
mingw_bin = r"C:\msys64\mingw64\bin"

os.environ["PATH"] = mingw_bin + ";" + scripts_dir + ";" + os.environ["PATH"]

# Verify tools
for tool in ["gfortran", "meson"]:
    import shutil
    loc = shutil.which(tool)
    if loc:
        print(f"Found {tool}: {loc}")
    else:
        print(f"ERROR: {tool} not found on PATH")
        sys.exit(1)

src = r"C:\a\hal\_shared\peakfqr\src"
sources = [
    os.path.join(src, "emafit.f"),
    os.path.join(src, "dcdflib1.f90"),
    os.path.join(src, "imslfake.f"),
    os.path.join(src, "probfun.f"),
]

build_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(build_dir)

cmd = [
    sys.executable, "-m", "numpy.f2py",
    "-c", *sources,
    "-m", "_emafort",
    "--backend", "meson",
    "--build-dir", os.path.join(build_dir, "mbuild"),
]

print(f"Running: {' '.join(cmd)}")
result = subprocess.run(cmd, env=os.environ)
sys.exit(result.returncode)
