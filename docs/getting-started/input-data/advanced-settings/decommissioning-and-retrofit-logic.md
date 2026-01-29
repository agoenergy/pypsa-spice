# PyPSA-SPICE Decommissioning & Retrofit Logic

This note explains how PyPSA-SPICE treats plant capacity that reaches the end of its technical lifetime **when you still want to allow the model to keep operating it by paying for a retrofit**. To represent this, each plant can be split into two **paired** assets:

- **Normal asset** (e.g., `XY_SO_COAL`) — the existing fleet/capacity  
- **Retrofit asset** (e.g., `XY_SO_COAL_RETRO`) — the portion that can be continued through retrofit investment  

Key assumption: the retrofit asset represents an **upgrade of existing capacity**, so it pays **only the retrofit CAPEX** (not the full cost of building a new plant).

## Logic overview

- **Base year (2025):** capacity that is already over-lifetime is **shifted** from the normal asset to the retrofit asset (treated as retrofitted, not decommissioned).
- **Future year (e.g., 2030):** capacity that becomes over-lifetime is written to `decomission_capacity.csv` so substations can **decommission** it from the normal fleet; the same amount becomes the **maximum retrofit investment option** via `retro.p_nom_max_2030`.

```mermaid
flowchart TB
  A["Inputs:
- normal p_nom
- commissioning year + lifetime
- base year 2025
- future year(s) e.g. 2030"] --> B{"End of lifetime in 2025?"}

  B -- Yes --> C["MW_over_2025"]
  C --> D["normal.p_nom -= MW_over_2025"]
  C --> E["retro.p_nom  += MW_over_2025
(assume retrofit, not decommission)"]

  B -- No --> F{"End of lifetime in 2030?"}
  D --> F
  E --> F

  F -- Yes --> G["MW_over_2030"]
  G --> H["Write decomission_capacity.csv:
(node, normal_asset, 2030, MW_over_2030)
→ substation can decommission from normal fleet"]
  G --> I["Set retrofit option:
retro.p_nom_max_2030 = MW_over_2030
→ model may invest (pay retrofit CAPEX)"]

  F -- No --> J["No decommision/retrofit done"]

  H --> K["Outputs:
- updated p_nom (2025 normal & retro)
- decomission_capacity.csv (future retirements)
- retrofit cap p_nom_max_year"]

  I --> K
  J --> K

  N["Assumptions:
- normal vs retrofit assets
- retrofit pays retrofit CAPEX only
- 2025: shift MW to retrofit p_nom
- 2030: decommissionable, but optional retrofit"] -.-> A
```
