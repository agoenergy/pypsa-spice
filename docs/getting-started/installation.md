<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Installation

## Clone the Repository

First, clone the [PyPSA-SPICE repository](https://github.com/agoenergy/pypsa-spice){:target="_blank"} using the **Git** version control system. **Important:** the path to the directory where the repository is cloned **must not contain any spaces**.

If Git is not installed on your system, please follow the [Git installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git){:target="_blank"}.

```shell title="Cloning the repository"
git clone https://github.com/agoenergy/pypsa-spice.git
cd pypsa-spice
```

## Install Python Dependencies

PyPSA-SPICE needs a set of Python packages to function. We recommend using **Conda**, a package and environment management system, to handle these dependencies.

Start by installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html){:target="_blank"}, a lightweight version of [Anaconda](https://www.anaconda.com/){:target="_blank"} that includes only Conda and its core dependencies. If you already have Conda installed, you can skip this step. For installation instructions tailored to your operating system, refer to the [official Conda installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/){:target="_blank"}.

Once Conda is installed, we recommend installing **Mamba**, a fast drop-in replacement for Conda that significantly speeds up environment installation. According to the [official mamba documentation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html){:target="_blank"}, the recommended way to install it is via **Miniforge3**, a minimal Conda installer bundled with Mamba.

Installation steps for Miniforge3 vary depending on your operating system. You can find platform-specific instructions in the [official Miniforge guide](https://github.com/conda-forge/miniforge){:target="_blank"}.

For example, if you’re using Linux system, you can install Miniforge3 by running the following command outside the PyPSA-SPICE folder and following the prompts:

```shell title="Installing Miniforge3 on Linux"
cd ..
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-$ Miniforge3-Linux-x86_64.sh
```

You can switch to any other directory outside the folder you cloned PyPSA-SPICE.

**Important Notes:**

1. **Mambaforge** is deprecated as of July 2024 and was officially retired after January 2025. For this reason, we recommend using Miniforge3 instead.
2. Install Miniforge3 or mamba packages only in the base Conda environment. Installing them in other environments may lead to compatibility issues or unexpected errors.
3. Miniforge3 and Mambaforge use different environment paths. If you switch from one to the other, you will need to recreate all Conda environments, as they are not shared between the two setups.

The required Python packages for PyPSA-SPICE are listed in the `environment.yaml` (via conda or mamba) and `requirements.txt` (via pip). You can create and activate the environment (which is called `pypsa-spice`) using the following commands:

```bash title="Installing and activating the virtual environment"
mamba env create -f envs/environment.yaml
conda activate pypsa-spice
```

Note that the environment activation is local to the currently open terminal session. If you open a new terminal window, you will need to re-run the activation command.

## Install a Solver

Th network model created in PyPSA-SPICE model builder is passed to an external solver to perform total annual system cost minimization and obtain optimal power flows. PyPSA-SPICE is compatible with several solvers that can be installed via Python. In our default environment, Gurobi and HiGHs solvers are installed (packages installed not licenses). Below is a list of supported solvers along with links to their official installation guides for different operating systems:

| Solver              | License type         | Installation guide |
| ------------------- | -------------------- | ------------------ |
| Ipopt               | Free & open-source   | [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html){:target="_blank"} |
| Cbc                 | Free & open-source   | [Cbc](https://github.com/coin-or/Cbc?tab=readme-ov-file#download){:target="_blank"} |
| GLPK                | Free & open-source   | [GLPK](https://www.gnu.org/software/glpk/){:target="_blank"} /[WinGLPK](http://winglpk.sourceforge.net/){:target="_blank"} |
| HiGHs               | Free & open-source   | [highspy](https://ergo-code.github.io/HiGHS/dev/interfaces/python/#Install){:target="_blank"} |
| Gurobi              | commercial | [gurobipy](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python){:target="_blank"} |
| CPLEX               | commercial | [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio){:target="_blank"} |
| FICO® Xpress Solver | commercial | [FICO® Xpress Solver](https://www.fico.com/de/products/fico-xpress-solver){:target="_blank"} |

Commercial solvers such as Gurobi currently significantly outperform open-source solvers for large-scale problems. In some cases, using a commercial solver may be necessary to obtain a feasible solution. However, open-source solvers are being improved recently, and some of them already perform very well for small to medium-sized problems. Please check the [solver benchmarking results](https://openenergybenchmark.org/dashboard/main-results){:target="_blank"} to compare their performance.

=== "on Mac or Linux"

    ```bash
    conda activate pypsa-spice
    conda install -c conda-forge ipopt coincbc
    ```

=== "on Windows"

    ```bash
    conda activate pypsa-spice
    conda install -c conda-forge ipopt glpk
    ```

On Windows, new versions of ``ipopt`` have caused problems. Consider downgrading to version 3.11.1.

## Set up the Default Configuration

PyPSA-SPICE requires several configuration options to be defined in a `config.yaml` file located in the root directory of the project. An example file, `config.default.yaml`, is provided in the repository as a reference.

Before running PyPSA-SPICE for the first time, create your own ``config.yaml`` by copying the example file:

```bash title="Creating the config.yaml file"
cp config.default.yaml config.yaml
```

Users are encouraged to regularly compare their own ``config.yaml`` with ``config.default.yaml`` when pulling updates from the remote repository, to ensure compatibility with any new or changed configuration options.

## Quick Execution of the Model Builder using Template Data

To have a first glance of how the model builder works, template data in `/root/data/global_csv_templates` folder can be used. After completing the installation, run the entire workflow with the following command:

```bash title="Running the entire workflow using this single command"
snakemake -j1 -c4 solve_all_networks # (1)!
```

1. `-j1` means running only 1 job in parallel at a time and `-c4` means allowing each job to use up to 4 CPU threads.

or

```bash title="Running the entire workflow using this call command"
snakemake -call
```

For more information, please refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html){:target="_blank"} to adjust the cores and threads to use.
