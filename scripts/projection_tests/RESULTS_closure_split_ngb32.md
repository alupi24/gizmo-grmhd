# Face non-closure versus reconstruction in GIZMO's face-sum divergence

`div(D)`: rms dx|D_i|/V_i/rms|B| of the full face-sum divergence (gradient std, limiters on). `div(B.C)`: same for the closure term B_i . sum_j A_ij. `div(R)`: same for the remainder R = D - B.C. `dface(exact)`: rms change of the face normal field B_face.n produced by the exact projection of D, relative to rms|B|; `dface(exact on R)`: the same if only R had to be removed (i.e. with closed faces).

| dist | N | field | div(D) | div(B.C) | div(R) | dface(exact) | dface(exact on R) |
|---|---|---|---|---|---|---|---|
| lattice | 16^3 | uniform | 5.505e-17 | 4.944e-17 | 8.016e-17 | 2.748e-17 | 3.251e-17 |
| lattice | 16^3 | tg3d | 2.134e-02 | 4.249e-17 | 2.134e-02 | 5.330e-03 | 5.330e-03 |
| lattice | 16^3 | abc | 1.363e-16 | 7.828e-17 | 1.430e-16 | 3.535e-17 | 6.713e-17 |
| lattice | 16^3 | grad_long | 1.344e-01 | 5.387e-17 | 1.344e-01 | 1.308e-01 | 1.308e-01 |
| lattice | 16^3 | bump | 3.862e-01 | 0.000e+00 | 3.862e-01 | 2.657e-01 | 2.657e-01 |
| lattice | 24^3 | uniform | 2.943e-15 | 2.945e-15 | 1.007e-16 | 7.932e-16 | 2.988e-17 |
| lattice | 24^3 | tg3d | 7.417e-03 | 3.026e-15 | 7.417e-03 | 2.048e-03 | 2.048e-03 |
| lattice | 24^3 | abc | 2.945e-15 | 2.947e-15 | 1.680e-16 | 7.858e-16 | 4.051e-17 |
| lattice | 24^3 | grad_long | 6.401e-02 | 2.875e-15 | 6.401e-02 | 7.568e-02 | 7.568e-02 |
| lattice | 24^3 | bump | 3.862e-01 | 1.913e-15 | 3.862e-01 | 2.713e-01 | 2.713e-01 |
| lattice | 32^3 | uniform | 5.209e-17 | 4.901e-17 | 8.746e-17 | 2.999e-17 | 3.971e-17 |
| lattice | 32^3 | tg3d | 3.532e-03 | 4.195e-17 | 3.532e-03 | 1.041e-03 | 1.041e-03 |
| lattice | 32^3 | abc | 1.323e-16 | 7.232e-17 | 1.421e-16 | 3.463e-17 | 1.370e-16 |
| lattice | 32^3 | grad_long | 4.169e-02 | 4.920e-17 | 4.169e-02 | 4.953e-02 | 4.953e-02 |
| lattice | 32^3 | bump | 3.862e-01 | 0.000e+00 | 3.862e-01 | 2.727e-01 | 2.727e-01 |
| jitter | 16^3 | uniform | 6.159e-01 | 6.159e-01 | 1.803e-16 | 1.082e-01 | 0.000e+00 |
| jitter | 16^3 | tg3d | 6.134e-01 | 6.185e-01 | 4.988e-02 | 1.070e-01 | 9.213e-03 |
| jitter | 16^3 | abc | 6.034e-01 | 6.045e-01 | 2.929e-02 | 1.063e-01 | 5.082e-03 |
| jitter | 16^3 | grad_long | 6.280e-01 | 6.151e-01 | 1.331e-01 | 1.542e-01 | 1.199e-01 |
| jitter | 16^3 | bump | 8.068e-01 | 7.197e-01 | 3.930e-01 | 2.600e-01 | 2.525e-01 |
| jitter | 24^3 | uniform | 6.060e-01 | 6.060e-01 | 1.767e-16 | 1.070e-01 | 0.000e+00 |
| jitter | 24^3 | tg3d | 6.025e-01 | 6.050e-01 | 3.035e-02 | 1.066e-01 | 5.513e-03 |
| jitter | 24^3 | abc | 6.204e-01 | 6.213e-01 | 1.938e-02 | 1.088e-01 | 3.389e-03 |
| jitter | 24^3 | grad_long | 6.071e-01 | 6.040e-01 | 6.202e-02 | 1.243e-01 | 6.828e-02 |
| jitter | 24^3 | bump | 6.865e-01 | 5.738e-01 | 4.031e-01 | 2.586e-01 | 2.613e-01 |
| jitter | 32^3 | uniform | 6.098e-01 | 6.098e-01 | 1.807e-16 | 1.071e-01 | 0.000e+00 |
| jitter | 32^3 | tg3d | 6.107e-01 | 6.123e-01 | 2.214e-02 | 1.072e-01 | 3.986e-03 |
| jitter | 32^3 | abc | 6.088e-01 | 6.092e-01 | 1.449e-02 | 1.073e-01 | 2.527e-03 |
| jitter | 32^3 | grad_long | 6.113e-01 | 6.100e-01 | 4.057e-02 | 1.149e-01 | 4.466e-02 |
| jitter | 32^3 | bump | 7.418e-01 | 6.414e-01 | 3.880e-01 | 2.600e-01 | 2.527e-01 |
| glass | 16^3 | uniform | 1.079e-01 | 1.079e-01 | 1.292e-16 | 2.000e-02 | 0.000e+00 |
| glass | 16^3 | tg3d | 1.139e-01 | 1.122e-01 | 2.210e-02 | 2.103e-02 | 4.130e-03 |
| glass | 16^3 | abc | 1.101e-01 | 1.102e-01 | 1.007e-02 | 2.037e-02 | 1.749e-03 |
| glass | 16^3 | grad_long | 1.729e-01 | 1.078e-01 | 1.347e-01 | 1.161e-01 | 1.148e-01 |
| glass | 16^3 | bump | 4.023e-01 | 1.048e-01 | 3.895e-01 | 2.356e-01 | 2.358e-01 |
| glass | 24^3 | uniform | 1.121e-01 | 1.121e-01 | 1.275e-16 | 2.063e-02 | 0.000e+00 |
| glass | 24^3 | tg3d | 1.123e-01 | 1.121e-01 | 1.267e-02 | 2.071e-02 | 2.303e-03 |
| glass | 24^3 | abc | 1.118e-01 | 1.118e-01 | 6.713e-03 | 2.054e-02 | 1.161e-03 |
| glass | 24^3 | grad_long | 1.285e-01 | 1.115e-01 | 6.400e-02 | 6.946e-02 | 6.659e-02 |
| glass | 24^3 | bump | 4.011e-01 | 1.000e-01 | 3.883e-01 | 2.415e-01 | 2.407e-01 |
| glass | 32^3 | uniform | 1.111e-01 | 1.111e-01 | 1.295e-16 | 2.054e-02 | 0.000e+00 |
| glass | 32^3 | tg3d | 1.095e-01 | 1.095e-01 | 8.965e-03 | 2.026e-02 | 1.594e-03 |
| glass | 32^3 | abc | 1.109e-01 | 1.109e-01 | 4.949e-03 | 2.050e-02 | 8.558e-04 |
| glass | 32^3 | grad_long | 1.186e-01 | 1.111e-01 | 4.157e-02 | 4.808e-02 | 4.355e-02 |
| glass | 32^3 | bump | 4.027e-01 | 1.034e-01 | 3.891e-01 | 2.425e-01 | 2.421e-01 |
