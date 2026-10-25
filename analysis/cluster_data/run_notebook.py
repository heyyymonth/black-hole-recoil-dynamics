"""Execute the data audit with this interpreter; no cluster-dynamics calculation."""
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = ROOT / "analysis/cluster_data/cluster_data_audit.ipynb"


def main():
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    manager = KernelManager(kernel_name="python3")
    manager.kernel_spec.argv = [sys.executable, "-m", "ipykernel_launcher",
                                "-f", "{connection_file}"]
    client = NotebookClient(notebook, km=manager, timeout=180,
                            resources={"metadata": {"path": str(ROOT)}},
                            record_timing=False)
    for cell in notebook.cells:
        cell.metadata.pop("execution", None)
    client.execute()
    nbformat.validate(notebook)
    nbformat.write(notebook, NOTEBOOK)
    print("Executed data-audit notebook successfully. No physical experiment run.")


if __name__ == "__main__":
    main()
