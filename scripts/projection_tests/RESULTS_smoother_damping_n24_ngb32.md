# Smoother damping at equal cost: N = 24^3, N_ngb = 32, GIZMO limiters on, gradient std

Entries: div(variant) / div(no projection) (smaller is better); cost in pair-pass equivalents (PPE).

## glass / weights tu / tg3d   (no projection: div = 1.133e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.475 | 0.252 | 0.147 | 0.063 | 0.016 | 0.007 | 0.002 |
| Jacobi omega=0.67 | 0.318 | 0.144 | 0.079 | 0.032 | 0.008 | 0.003 | 0.001 |
| Jacobi omega=0.80 | 0.223 | 0.097 | 0.052 | 0.020 | 0.005 | 0.002 | 0.001 |
| Jacobi omega=1.00 | 0.230 | 0.078 | 0.034 | 0.012 | 0.003 | 0.001 | 0.000 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.124 | 2 passes (29.3 PPE): 0.021 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.224 | 2 passes (29.3 PPE): 0.063 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.407 | 2 passes (29.3 PPE): 0.176 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.504 | 2 passes (29.3 PPE): 0.261 | | | | | |

## glass / weights tu / grad_long   (no projection: div = 1.295e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.637 | 0.521 | 0.479 | 0.444 | 0.398 | 0.361 | 0.278 |
| Jacobi omega=0.67 | 0.554 | 0.481 | 0.456 | 0.426 | 0.372 | 0.329 | 0.238 |
| Jacobi omega=0.80 | 0.515 | 0.466 | 0.445 | 0.414 | 0.354 | 0.307 | 0.212 |
| Jacobi omega=1.00 | 0.511 | 0.455 | 0.432 | 0.397 | 0.328 | 0.277 | 0.180 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.429 | 2 passes (29.3 PPE): 0.352 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.470 | 2 passes (29.3 PPE): 0.380 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.566 | 2 passes (29.3 PPE): 0.430 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.628 | 2 passes (29.3 PPE): 0.473 | | | | | |

## glass / weights tu / bump   (no projection: div = 4.092e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.887 | 0.819 | 0.765 | 0.678 | 0.522 | 0.419 | 0.252 |
| Jacobi omega=0.67 | 0.857 | 0.779 | 0.716 | 0.615 | 0.446 | 0.344 | 0.195 |
| Jacobi omega=0.80 | 0.837 | 0.750 | 0.680 | 0.571 | 0.398 | 0.299 | 0.164 |
| Jacobi omega=1.00 | 0.811 | 0.711 | 0.632 | 0.514 | 0.340 | 0.248 | 0.130 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.553 | 2 passes (29.3 PPE): 0.377 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.618 | 2 passes (29.3 PPE): 0.448 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.700 | 2 passes (29.3 PPE): 0.536 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.745 | 2 passes (29.3 PPE): 0.590 | | | | | |

## glass / weights tpfa / tg3d   (no projection: div = 1.133e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.474 | 0.249 | 0.144 | 0.061 | 0.015 | 0.006 | 0.002 |
| Jacobi omega=0.67 | 0.315 | 0.141 | 0.077 | 0.031 | 0.007 | 0.003 | 0.001 |
| Jacobi omega=0.80 | 0.217 | 0.094 | 0.050 | 0.019 | 0.005 | 0.002 | 0.001 |
| Jacobi omega=1.00 | 0.221 | 0.072 | 0.031 | 0.011 | 0.003 | 0.001 | 0.000 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.121 | 2 passes (29.3 PPE): 0.021 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.223 | 2 passes (29.3 PPE): 0.062 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.407 | 2 passes (29.3 PPE): 0.175 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.504 | 2 passes (29.3 PPE): 0.261 | | | | | |

## glass / weights tpfa / grad_long   (no projection: div = 1.295e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.636 | 0.519 | 0.477 | 0.443 | 0.396 | 0.359 | 0.275 |
| Jacobi omega=0.67 | 0.553 | 0.480 | 0.454 | 0.425 | 0.370 | 0.327 | 0.235 |
| Jacobi omega=0.80 | 0.512 | 0.465 | 0.444 | 0.413 | 0.351 | 0.304 | 0.209 |
| Jacobi omega=1.00 | 0.508 | 0.453 | 0.431 | 0.395 | 0.326 | 0.274 | 0.177 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.428 | 2 passes (29.3 PPE): 0.351 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.469 | 2 passes (29.3 PPE): 0.380 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.566 | 2 passes (29.3 PPE): 0.429 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.628 | 2 passes (29.3 PPE): 0.473 | | | | | |

## glass / weights tpfa / bump   (no projection: div = 4.092e-01; stencil pass = 14.7 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.886 | 0.816 | 0.762 | 0.674 | 0.517 | 0.413 | 0.247 |
| Jacobi omega=0.67 | 0.856 | 0.776 | 0.712 | 0.610 | 0.441 | 0.338 | 0.190 |
| Jacobi omega=0.80 | 0.835 | 0.747 | 0.677 | 0.566 | 0.392 | 0.294 | 0.159 |
| Jacobi omega=1.00 | 0.809 | 0.707 | 0.628 | 0.508 | 0.334 | 0.243 | 0.127 |
| stencil omega_S=1.0 | 1 pass (14.7 PPE): 0.551 | 2 passes (29.3 PPE): 0.375 | | | | | |
| stencil omega_S=0.8 | 1 pass (14.7 PPE): 0.617 | 2 passes (29.3 PPE): 0.446 | | | | | |
| stencil omega_S=0.6 | 1 pass (14.7 PPE): 0.699 | 2 passes (29.3 PPE): 0.535 | | | | | |
| stencil omega_S=0.5 | 1 pass (14.7 PPE): 0.744 | 2 passes (29.3 PPE): 0.589 | | | | | |

## jitter / weights tu / tg3d   (no projection: div = 6.129e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.481 | 0.260 | 0.158 | 0.077 | 0.027 | 0.014 | 0.005 |
| Jacobi omega=0.67 | 0.326 | 0.155 | 0.093 | 0.046 | 0.016 | 0.008 | 0.003 |
| Jacobi omega=0.80 | 0.232 | 0.111 | 0.068 | 0.034 | 0.012 | 0.006 | 0.002 |
| Jacobi omega=1.00 | 0.231 | 0.090 | 0.049 | 0.023 | 0.008 | 0.004 | 0.002 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.123 | 2 passes (31.0 PPE): 0.023 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.237 | 2 passes (31.0 PPE): 0.068 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.418 | 2 passes (31.0 PPE): 0.184 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.514 | 2 passes (31.0 PPE): 0.270 | | | | | |

## jitter / weights tu / grad_long   (no projection: div = 6.167e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.493 | 0.280 | 0.185 | 0.117 | 0.082 | 0.072 | 0.056 |
| Jacobi omega=0.67 | 0.345 | 0.184 | 0.130 | 0.096 | 0.075 | 0.066 | 0.048 |
| Jacobi omega=0.80 | 0.256 | 0.146 | 0.111 | 0.088 | 0.071 | 0.061 | 0.043 |
| Jacobi omega=1.00 | 0.250 | 0.127 | 0.098 | 0.081 | 0.066 | 0.056 | 0.037 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.146 | 2 passes (31.0 PPE): 0.070 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.249 | 2 passes (31.0 PPE): 0.098 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.425 | 2 passes (31.0 PPE): 0.198 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.518 | 2 passes (31.0 PPE): 0.280 | | | | | |

## jitter / weights tu / bump   (no projection: div = 7.128e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.642 | 0.509 | 0.449 | 0.388 | 0.301 | 0.244 | 0.150 |
| Jacobi omega=0.67 | 0.555 | 0.455 | 0.410 | 0.352 | 0.260 | 0.203 | 0.117 |
| Jacobi omega=0.80 | 0.507 | 0.431 | 0.388 | 0.328 | 0.233 | 0.177 | 0.098 |
| Jacobi omega=1.00 | 0.485 | 0.407 | 0.361 | 0.297 | 0.200 | 0.148 | 0.079 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.313 | 2 passes (31.0 PPE): 0.200 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.385 | 2 passes (31.0 PPE): 0.244 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.514 | 2 passes (31.0 PPE): 0.324 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.588 | 2 passes (31.0 PPE): 0.388 | | | | | |

## jitter / weights tpfa / tg3d   (no projection: div = 6.129e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.484 | 0.269 | 0.169 | 0.086 | 0.032 | 0.017 | 0.005 |
| Jacobi omega=0.67 | 0.335 | 0.168 | 0.104 | 0.054 | 0.019 | 0.010 | 0.003 |
| Jacobi omega=0.80 | 0.249 | 0.124 | 0.077 | 0.040 | 0.014 | 0.007 | 0.003 |
| Jacobi omega=1.00 | 0.256 | 0.106 | 0.059 | 0.028 | 0.010 | 0.005 | 0.002 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.128 | 2 passes (31.0 PPE): 0.023 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.239 | 2 passes (31.0 PPE): 0.070 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.419 | 2 passes (31.0 PPE): 0.186 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.514 | 2 passes (31.0 PPE): 0.272 | | | | | |

## jitter / weights tpfa / grad_long   (no projection: div = 6.167e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.496 | 0.289 | 0.196 | 0.125 | 0.085 | 0.074 | 0.058 |
| Jacobi omega=0.67 | 0.353 | 0.196 | 0.140 | 0.101 | 0.077 | 0.067 | 0.050 |
| Jacobi omega=0.80 | 0.272 | 0.157 | 0.118 | 0.091 | 0.073 | 0.063 | 0.045 |
| Jacobi omega=1.00 | 0.276 | 0.140 | 0.104 | 0.084 | 0.067 | 0.057 | 0.039 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.150 | 2 passes (31.0 PPE): 0.070 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.251 | 2 passes (31.0 PPE): 0.099 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.425 | 2 passes (31.0 PPE): 0.199 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.519 | 2 passes (31.0 PPE): 0.281 | | | | | |

## jitter / weights tpfa / bump   (no projection: div = 7.128e-01; stencil pass = 15.5 PPE)

| smoother | 1 PPE | 2 PPE | 3 PPE | 5 PPE | 10 PPE | 15 PPE | 30 PPE |
|---|---|---|---|---|---|---|---|
| Jacobi omega=0.50 | 0.634 | 0.507 | 0.451 | 0.393 | 0.309 | 0.253 | 0.158 |
| Jacobi omega=0.67 | 0.549 | 0.458 | 0.415 | 0.359 | 0.268 | 0.211 | 0.124 |
| Jacobi omega=0.80 | 0.507 | 0.435 | 0.394 | 0.336 | 0.242 | 0.186 | 0.105 |
| Jacobi omega=1.00 | 0.499 | 0.415 | 0.368 | 0.305 | 0.209 | 0.156 | 0.084 |
| stencil omega_S=1.0 | 1 pass (15.5 PPE): 0.314 | 2 passes (31.0 PPE): 0.199 | | | | | |
| stencil omega_S=0.8 | 1 pass (15.5 PPE): 0.387 | 2 passes (31.0 PPE): 0.244 | | | | | |
| stencil omega_S=0.6 | 1 pass (15.5 PPE): 0.515 | 2 passes (31.0 PPE): 0.325 | | | | | |
| stencil omega_S=0.5 | 1 pass (15.5 PPE): 0.589 | 2 passes (31.0 PPE): 0.389 | | | | | |

