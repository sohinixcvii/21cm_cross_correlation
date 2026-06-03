# Installing 21cmFAST v4.1.1 on CSD3 (or Similar HPC Systems)

This document records the installation procedure that successfully worked after encountering cache, quota, and dependency issues.

---

# 1. Create a Clean Conda Environment

```bash
CONDA_NO_PLUGINS=true conda create -n 21cmfast python=3.11 -y
conda activate 21cmfast
```

---

# 2. Configure Cache Directories

On CSD3, home quotas are small and pip may fail with:

```text
OSError: [Errno 122] Disk quota exceeded
```

Create cache directories on scratch storage:

```bash
mkdir -p /nvme2/scratch/$USER/.cache/pip
mkdir -p /nvme2/scratch/$USER/.cache
```

Set environment variables:

```bash
export PIP_CACHE_DIR=/nvme2/scratch/$USER/.cache/pip
export XDG_CACHE_HOME=/nvme2/scratch/$USER/.cache
```

(Optional) Add these to `~/.bashrc`.

---

# 3. Install Core Dependencies

Install FFTW and GSL from conda-forge:

```bash
CONDA_NO_PLUGINS=true conda install -c conda-forge fftw gsl
```

Verify FFTW is present:

```bash
ls $CONDA_PREFIX/lib/libfftw3f*
```

Expected output includes:

```text
libfftw3f.so
libfftw3f_omp.so
```

---

# 4. Install 21cmFAST

```bash
pip install --no-cache-dir 21cmFAST==4.1.1
```

The installation may take several minutes because C extensions are compiled.

---

# 5. Verify Installation

Check that the package imports correctly:

```bash
python -c "import py21cmfast; print(py21cmfast.__version__)"
```

Expected output:

```text
4.1.1
```

---

# Common Problems

## Problem 1: Disk Quota Exceeded

Error:

```text
OSError: [Errno 122] Disk quota exceeded
```

Solution:

Move pip cache to scratch space:

```bash
export PIP_CACHE_DIR=/nvme2/scratch/$USER/.cache/pip
export XDG_CACHE_HOME=/nvme2/scratch/$USER/.cache
```

---

## Problem 2: Conda Plugin Error

Error:

```text
Error while loading conda entry point
conda-libmamba-solver
```

Solution:

Run conda with plugins disabled:

```bash
CONDA_NO_PLUGINS=true conda install ...
```

or

```bash
CONDA_NO_PLUGINS=true conda create ...
```

---

## Problem 3: FFTW Linking Failure

Error:

```text
cannot find -lfftw3f
cannot find -lfftw3f_omp
```

Cause:

FFTW libraries are missing.

Fix:

```bash
CONDA_NO_PLUGINS=true conda install -c conda-forge fftw gsl
```

Verify:

```bash
ls $CONDA_PREFIX/lib/libfftw3f*
```

---

## Problem 4: Build Fails During Compilation

First check:

```bash
which gcc
which g++
```

and

```bash
echo $CONDA_PREFIX
```

Make sure the conda environment is activated before installation.

---

# Recreating Environment from Scratch

```bash
CONDA_NO_PLUGINS=true conda create -n 21cmfast python=3.11 -y

conda activate 21cmfast

export PIP_CACHE_DIR=/nvme2/scratch/$USER/.cache/pip
export XDG_CACHE_HOME=/nvme2/scratch/$USER/.cache

CONDA_NO_PLUGINS=true conda install -c conda-forge fftw gsl

pip install --no-cache-dir 21cmFAST==4.1.1

python -c "import py21cmfast; print(py21cmfast.__version__)"
```

If the final command prints:

```text
4.1.1
```

the installation has succeeded.

