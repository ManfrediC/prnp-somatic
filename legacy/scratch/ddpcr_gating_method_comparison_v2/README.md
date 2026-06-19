# ddPCR gating method comparison v2

This scratch experiment implements the control-anchored method-comparison plan
from `scratch/ddpcr_gating_method_comparison/PLAN_v2_control_anchored_methods.md`.

The official workflow in `src/ddpcr` and immutable raw files in `raw/ddpcr`
remain unchanged. Files under `inputs/` are derived package-specific exports
from the raw JSONs so that each package can be run in its expected format.

Numerical outputs must be computed from full droplet classifications. Plots may
use deterministic downsampling only after classification and counting.

