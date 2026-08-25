# Project Instructions for Claude Code

## Environment

Always use the `21cmfast` conda environment for all Python commands.
Activate it before running anything:

```bash
conda run -n 21cmfast <command>
```

Never use the system Python or pip directly. If a package is missing,
install it into the conda environment:

```bash
conda run -n 21cmfast pip install <package>
# or
conda install -n 21cmfast <package>
```

## Running Code

Always do a test run before considering any task complete.

- For scripts: `conda run -n 21cmfast python <script.py>`
- For tests: `conda run -n 21cmfast pytest tests/ -v`

If a test run fails, diagnose and fix before reporting back. Do not
hand back broken code.

## Code Style

- Follow PEP 8.
- Use type hints for all function signatures.
- Write docstrings for all functions and classes (NumPy docstring format).
- Keep functions small and single-purpose.
- Prefer explicit over implicit — no clever one-liners that obscure intent.

## Testing

- Every new function should have at least one corresponding test in `tests/`.
- Tests go in files named `test_<module>.py`.
- Run the full test suite after any non-trivial change.
- Do not mark a task complete if any tests are failing.

## General Rules

- Always read existing code before editing — do not assume structure.
- Do not delete files. Move to `_archive/` if something needs to be removed.
- Figures are generated programmatically. Edit the script that produces them,
  not the output files directly.
- When in doubt about a destructive action (overwriting, bulk reformatting,
  deleting), ask before proceeding.
- Do not commit anything. 
- Always update README files with latest information.
- Always update CHANGELOG with latest changes.
- Whenever any numbers or formulae are updated, update NUMBERS_AND_RESOURCES.md to reflect the latest changes