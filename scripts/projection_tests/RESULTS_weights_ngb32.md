# Weighting of the projection: Tu et al. versus two-point (TPFA)

Exact solves. `recovery` (gradient fields): rms|B_face.n - B0.n| / rms|grad(phi).n|, 0 = gradient part fully removed. `face err` (tg3d): rms|B_face.n - B_exact.n| / rms|B|. `rhs = D`: full face-sum divergence; `rhs = R`: closure term removed (closed-face proxy). Field widths are fixed in units of the mean spacing (grad_short: 1.5 dx, grad_long: 4 dx), so the gradient part is equally resolved at every N.

| dist | N | field | no projection | Tu, rhs=D | TPFA, rhs=D | Tu, rhs=R | TPFA, rhs=R |
|---|---|---|---|---|---|---|---|
| lattice | 16^3 | grad_short | 1.019 | 0.495 | 0.088 | 0.495 | 0.088 |
| lattice | 16^3 | grad_long | 1.011 | 0.495 | 0.023 | 0.495 | 0.023 |
| lattice | 16^3 | tg3d | 0.023 | 0.022 | 0.022 | 0.022 | 0.022 |
| lattice | 24^3 | grad_short | 1.019 | 0.495 | 0.088 | 0.495 | 0.088 |
| lattice | 24^3 | grad_long | 1.015 | 0.497 | 0.023 | 0.497 | 0.023 |
| lattice | 24^3 | tg3d | 0.010 | 0.009 | 0.009 | 0.009 | 0.009 |
| lattice | 32^3 | grad_short | 1.019 | 0.495 | 0.088 | 0.495 | 0.088 |
| lattice | 32^3 | grad_long | 1.015 | 0.497 | 0.023 | 0.497 | 0.023 |
| lattice | 32^3 | tg3d | 0.005 | 0.005 | 0.005 | 0.005 | 0.005 |
| glass | 16^3 | grad_short | 1.013 | 0.839 | 0.630 | 0.664 | 0.122 |
| glass | 16^3 | grad_long | 1.011 | 0.685 | 0.161 | 0.675 | 0.078 |
| glass | 16^3 | tg3d | 0.023 | 0.031 | 0.034 | 0.023 | 0.023 |
| glass | 24^3 | grad_short | 1.010 | 1.170 | 1.164 | 0.673 | 0.124 |
| glass | 24^3 | grad_long | 1.017 | 0.715 | 0.273 | 0.680 | 0.062 |
| glass | 24^3 | tg3d | 0.011 | 0.024 | 0.027 | 0.011 | 0.011 |
| glass | 32^3 | grad_short | 1.015 | 1.630 | 1.792 | 0.681 | 0.129 |
| glass | 32^3 | grad_long | 1.017 | 0.760 | 0.415 | 0.680 | 0.062 |
| glass | 32^3 | tg3d | 0.006 | 0.022 | 0.026 | 0.007 | 0.007 |
