# Project notes

## GWF pipeline

This project uses [GWF](https://gwf.app/) to orchestrate analyses on Slurm and locally. The pipeline lives in `workflow.py`.

For a compact local reference to the GWF API, CLI, and common patterns used here, see [`claude-gwf-ref.md`](./claude-gwf-ref.md). Consult it before writing or modifying targets, templates, `gwf.map`/`collect` chains, or backend/executor configuration.
