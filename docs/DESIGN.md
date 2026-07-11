# PoreSim Visual Design

## Color Palette

The logo and all visual assets use a strict 3-color amber palette mirroring the
[Material Design Amber scale](https://m2.material.io/design/color/the-color-system.html),
following the same structural logic as PoreMS (blue) and PoreAna (red).

| Role | Hex | RGB | Material token |
|---|---|---|---|
| Bond lines | `#ECD078` | `rgb(236, 208, 120)` | Amber 200 |
| Small/medium nodes, accent text ("SIM") | `#CA8B32` | `rgb(202, 139, 50)` | Amber 700 |
| Large anchor nodes, body text ("PORE") | `#92400E` | `rgb(146, 64, 14)` | Amber 900 |

### Rationale

- **3 colors maximum** keeps the graphic legible at small sizes (favicon, sidebar logo).
- **Perceptual hierarchy** matches PoreMS exactly: bonds are the lightest element (connectors),
  large anchor nodes are the darkest (structural weight).
- **Material Amber scale** provides perceptually uniform steps and is visually distinct from
  PoreMS (blue) and PoreAna (red) while sharing the same design grammar.

## Logo Files

| File | Usage |
|---|---|
| `docs/pics/logo.svg` | Square icon (favicon source, app icon) |
| `docs/pics/logo_text.svg` | Horizontal logo with "PoreSim" wordmark |
| `docs/pics/logo_text_sub.svg` | Logo with wordmark + subtitle line |

## Typography

Wordmark uses **Arial / Arial MT** (system sans-serif fallback).
- "**PORE**": Amber 900 (`#92400E`) — darkest, matches large anchor nodes
- "**SIM**": Amber 700 (`#CA8B32`) — medium, matches small/medium nodes

## Favicon

The favicon is derived from `logo.svg`. At 32 × 32 px the three-level amber scheme remains
distinguishable because the large anchor nodes (Amber 900) contrast against the bond lines (Amber 200).
