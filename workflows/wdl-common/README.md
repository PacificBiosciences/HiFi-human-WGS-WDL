<h1 align="center"><img width="300px" src="./images/logo_wdl_workflows.svg" alt="PacBio Pipelines - Common Tasks and Workflows"/></h1>

<h1 align="center">PacBio Pipelines<br/>Common Tasks and Workflows</h1>

Workflows and tasks reused across PacBio workflows.

---

**The WDL files here are under active development and are currently provided in an unsupported format.**

---

## Developer setup

### Prerequisites

1. **uv** (Python package manager):

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. **Rust toolchain** (for sprocket):

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

3. **Sprocket** (WDL checker, linter, and formatter):

   ```bash
   cargo install sprocket
   ```

4. **Python virtual environment and dependencies**:

   ```bash
   uv venv --managed-python --python 3.13
   source .venv/bin/activate
   uv pip install -r requirements.txt
   ```

### Install pre-commit hooks

```bash
make setup
```

This installs git pre-commit hooks that run `sprocket check`, `sprocket lint`, and
`sprocket format check` on staged `.wdl` files before each commit.

To run the checks manually against all WDL files:

```bash
make pre-commit
```

If the format check fails, fix formatting with:

```bash
make format-overwrite
```

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
