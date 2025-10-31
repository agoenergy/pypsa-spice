<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# How to Contribute

PyPSA-SPICE is an open-source repository and thus your contributions are welcome and appreciated!

If you are considering using this model builder, please reach out to us at [modelling@agora-thinktanks.org](mailto:modelling@agora-thinktanks.org){:target="_blank"}. We would be happy to help you get started.

If you encounter a bug, please create a [new issue](https://github.com/agoenergy/pypsa-spice/issues){:target="_blank"}. For new ideas or feature requests, you can start a conversation in the [discussions](https://github.com/agoenergy/pypsa-spice/discussions){:target="_blank"} section of the repository. For troubleshooting, please check the [troubleshooting](../user-guide/troubleshooting.md) guide for more information.

## Code style

If your contributions involve code changes, please follow the steps below to format your code before creating a pull request. This will help us review your changes more efficiently!

1. Fork the repository on GitHub
2. Clone your fork: `git clone https://github.com/<your-username>/pypsa-spice.git`
3. Fetch the upstream tags `git fetch --tags https://github.com/agoenergy/pypsa-spice/pypsa-spice.git`
4. Install with dependencies in editable mode: `pip install -e .[dev]`
5. Setup linter and formatter, e.g `pre-commit install` (see [Pre-Commit](#pre-commit))
6. Open a new branch and write your code
7. Push your changes to your fork and create a pull request on GitHub

## Pre-Commit

We use [pre-commit](https://pre-commit.com/) to maintain consistent code style and catch common errors before commits. The pre-commit package is included in the `pypsa-spice` environment. To enable automatic formatting before each commit (recommended), run this command once:

```bash title="Initialize pre-commit"
pre-commit install
```

This will automatically check the changes which are staged before you commit them.
