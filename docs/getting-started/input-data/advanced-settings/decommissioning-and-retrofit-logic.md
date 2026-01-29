# Advanced settings in PyPSA-SPICE

## Decommissioning & retrofit logic

This note explains how PyPSA‑SPICE handles plant capacity that reaches the end of its technical lifetime **while you still want to allow the model to keep operating the plant by paying for a retrofit**. To represent this, each plant can be split into two **paired** assets:

- **Normal asset** (e.g., `XY_SO_COAL`) — the existing fleet/capacity  
- **Retrofit asset** (e.g., `XY_SO_COAL_RETRO`) — the portion that can be continued through retrofit investment  

!!! Note
    Key assumption: retrofit asset is treated as an **upgrade of the existing capacity**, so it pays **only the retrofit CAPEX** rather than the full greenfield cost of constructing a new plant.

### Logic overview

- **Base year (2025):** capacity that is already over-lifetime is **shifted** from the normal asset to the retrofit asset (treated as retrofitted, not decommissioned).

- **Future year (e.g., 2030):** capacity that becomes over-lifetime is written to `decomission_capacity.csv`. Hence, substations can **decommission** it from the normal fleet; the same amount becomes the **maximum retrofit investment option** via `retro.p_nom_max_2030`.

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontSize":"11px", "htmlLabels": true}}}%%
flowchart TB
A[["<div style='text-align:left; margin:0'><u>Inputs</u>
- normal p_nom
- commissioning year + lifetime
- base year 2025
- future year(s) e.g. 2030</div>"]] --> B{"End of lifetime in 2025?"}

B -- Yes --> C["MW_over_2025"]
C --> D["normal.p_nom -= MW_over_2025"]
C --> E["retro.p_nom  += MW_over_2025
(assume retrofit, not decommission)"]

B -- No --> F{"End of lifetime in 2030?"}
D --> F
E --> F

F -- Yes --> G["MW_over_2030"]
G --> H["<div style='text-align:left; margin:0'><u>Write decomission_capacity.csv</u>
(node, normal_asset, 2030, MW_over_2030)
→ substation can decommission from normal fleet</div>"]
G --> I["<div style='text-align:left; margin:0'><u>Set retrofit option</u>
retro.p_nom_max_2030 = MW_over_2030
→ model may invest (pay retrofit CAPEX)</div>"]

F -- No --> J["No decommision/retrofit done"]

H --> K[["<div style='text-align:left; margin:0'><u>Outputs</u>
- updated p_nom (2025 normal & retro)
- decomission_capacity.csv (future retirements)
- retrofit cap p_nom_max_year</div>"]]

I --> K
J --> K

N{{"<div style='text-align:left; margin:0'><u>Assumptions</u>
- normal vs retrofit assets
- retrofit assets pay retrofit CAPEX only
- 2025: shift MW to p_nom in retrofit assets
- 2030: decommissionable, but optional retrofit</div>"}} -.-> A

```
