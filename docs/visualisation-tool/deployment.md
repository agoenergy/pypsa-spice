<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Deployment of PyPSA-SPICE-Vis app

The easiest method to deploy is using [git-lfs](https://git-lfs.com/){target="_blank"} to manage the data and simply using [streamlit's community cloud](https://streamlit.io/cloud){target="_blank"}.

0. Download git-lfs:

=== "on Windows"

    ```bash
    download git-lfs directly from https://git-lfs.com/.
    ```

=== "on Linux"

    ```bash
    sudo apt-get install git-lfs
    ```

=== "On Mac"

    ```bash
    brew install git-lfs
    ```

1. To deploy the app, create a branch of the repo. This will also allow you to make changes like default values for country and scenario, graph setting, etc.

2. Install [git-lfs](https://git-lfs.com/){target="_blank"} via the following command:

```shell title="installing git-lfs"
git lfs install
```

3. Add the results data files inside this branch. These data files will be tracked using git-lfs. Note: add only the CSV files but do keep the same folder structure as that of `pypsa-spice` results.

4. Add path to this results folder inside `pypsa-spice-vis/setting/initial_project_01_deploy.yaml`. Note: the `pypsa-spice-vis/setting/initial.yaml`is ignored by git and used for local deployment.

5. In `pypsa-spice-vis/main.py`, make `DEPLOY == True`. By default, `False` for local run.

6. Initialise and track the result files using [git-lfs](https://git-lfs.com/){target="_blank"}. You will have to add all csvs to be tracked with lfs tracking system with following commands based on your operation system:

=== "on Mac or Linux"

    ```bash
    git lfs install
    git lfs track "*.csv"
    git add .gitattributes
    git commit -m "Track all CSV files with Git LFS"
    find . -name "*.csv" -exec git add {} \; 
    ```

=== "on Windows"

    ```bash
    git lfs install
    git lfs track "*.csv"
    git add .gitattributes
    git commit -m "Track all CSV files with Git LFS"
    for /r %i in (*.csv) do git add "%i"
    ```

After adding all csv files, you can initialise and commit tracked files using git-lfs:

```shell title="committing all tracked files"
git commit -m "Add all CSV files"
git push
```

Now the repository can be linked simply to Streamlit's deployment system.
