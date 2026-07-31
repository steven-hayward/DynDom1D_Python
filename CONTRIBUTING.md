# Contributing to DynDom1D_Python

Thanks for your interest in contributing!

> **Looking for the version described in the paper?** Please use the [official release](../../releases) rather than `main`, which may have changed since publication.

## Reporting bugs or suggesting features

Feel free to [open an issue](../../issues) with:
- Steps to reproduce (including the command or input used)
- Expected vs. actual behavior
- Your Python version / OS

## Contributing code

We use the standard fork-and-pull-request workflow:

1. Fork the repository and clone it locally.
2. Create a branch for your change:
   ```bash
   git checkout -b fix/short-description
   ```
3. Install dependencies:
   ```bash
   cd src/DynDom1D_Python
   pip install -r requirements.txt
   ```
4. Make your changes, and test them against the example datasets in `data/` to confirm expected output in `output/`.
5. Commit with a clear message and push to your fork.
6. Open a pull request against `main`, describing what changed and why (referencing any related issue).

Pull requests are automatically checked by our CI (GitHub Actions), which runs the test suite — please make sure these checks pass before requesting review.

If you're planning a larger change, it's worth forking and trying things out first to get a feel for the codebase (`Clusterer.py`, `Domain.py`, `Protein.py`, `HierarchySystem.py`, `DomainBuilder.py`, `FileMngr.py`) before opening a PR.

## Code style

Please follow the existing conventions in the file you're editing, and add docstrings/comments for new functions and classes.

## Getting help

If you have questions about using DynDom1D_Python, or run into trouble that isn't a bug report, please [open an issue](../../issues) as well — this is the best way to reach us and helps others who run into the same question.

## License

By contributing, you agree that your contributions will be licensed under this project's license (see `LICENSE`).
