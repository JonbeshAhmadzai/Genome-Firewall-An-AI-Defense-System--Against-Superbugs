from __future__ import annotations

import importlib.util
import shutil


PYTHON_MODULES = [
    "pandas",
    "numpy",
    "scipy",
    "sklearn",
    "matplotlib",
    "seaborn",
    "plotly",
    "Bio",
    "pyarrow",
    "joblib",
    "tqdm",
    "yaml",
    "requests",
    "streamlit",
    "gradio",
    "shap",
    "reportlab",
    "pypdf",
    "fitz",
    "pdfplumber",
]

CLI_TOOLS = ["amrfinder", "mmseqs", "mash", "fastANI", "fastani"]


def main() -> None:
    print("Python modules")
    for module in PYTHON_MODULES:
        status = "OK" if importlib.util.find_spec(module) else "MISSING"
        print(f"{status:8} {module}")

    print("\nCLI tools")
    for command in CLI_TOOLS:
        path = shutil.which(command)
        status = "OK" if path else "MISSING"
        print(f"{status:8} {command}" + (f" -> {path}" if path else ""))


if __name__ == "__main__":
    main()

