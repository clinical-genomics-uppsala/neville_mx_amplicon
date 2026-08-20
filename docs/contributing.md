# Contributing

We welcome contributions to improve the `neville_mx_amplicon` pipeline. Please follow these guidelines to submit issues, add features, and maintain code quality.

---

## 1. Development Standards

To keep the pipeline maintainable, we adhere to the following development practices:
- **Snakemake Rule Design**: Follow the modular rules structure defined by the [hydra-genetics](https://github.com/hydra-genetics) guidelines.
- **Containerization**: All rules must run within suitable Singularity/Apptainer containers to ensure reproducibility.
- **Python Style**: Python scripts under `workflow/scripts/` should follow PEP 8 standards.

---

## 2. Formatting

Before pushing code, format rules and scripts:
- **Snakefiles/Rules**: Format using `snakefmt`.
- **Python code**: Format using `black` or check using `pycodestyle`.

---

## 3. Pull Requests & Code Review

1. Fork the repository and create a new feature branch.
2. Ensure the code is linted and all integration tests complete successfully.
3. Submit a Pull Request describing your changes and link any related issues.
