<!--
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Repository agent guidance

## PyPSA-SPICE web visualisation

- Treat `pypsa-spice-web/` as an independent web application built with React,
  TypeScript, CSS, and FastAPI.
- For work on `pypsa-spice-vis-ui` or `pypsa-spice-web/`, do not use the
  Streamlit application in `pypsa-spice-vis/` as UI or UX guidance. Do not copy
  or port its page structure, layouts, controls, component choices, styling, or
  interaction patterns into the web interface.
- Use the current React frontend, the web visualisation documentation, and the
  user's explicit requirements as the sources of truth for web UI decisions.
- The Streamlit application may be inspected only when necessary to understand
  data semantics, calculations, or feature behaviour. Reimplement those needs
  in web-native patterns rather than reproducing the Streamlit interface.
- Keep Streamlit-specific packages and implementation details out of the web
  frontend and its UI-development guidance.
- The web application does not need mobile optimization. Prioritize desktop
  layouts and interactions; mobile-specific responsive work is not required
  unless the user explicitly requests it.
- Don't apply the Agora design skill. Instead, keep the app's own design consistent.
