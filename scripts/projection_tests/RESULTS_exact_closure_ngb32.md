# Exact face closure by an antisymmetric correction A~_ij = A_ij - k_ij (lambda_i - lambda_j)

`closure` = rms |sum_j A_ij| / sum_j |A_ij|; `dA` = rms |A~ - A| / rms |A|; `flipped` = faces whose orientation reverses; `lin. cons.` = rms || sum_j 1/2 (x_j - x_i) (x) A_ij - V_i I || / V_i; `div tg3d` = face-sum divergence of the divergence-free tg3d field, rms dx|D|/V/rms|B| (no B projection); `Jacobi 10x/100x` = sweeps to reduce the closure 10x / 100x from a cold start (k = |A|).

| set | N | faces | closure | dA | flipped | lin. cons. | div tg3d | Jacobi 10x / 100x |
|---|---|---|---|---|---|---|---|---|
| glass | 16^3 | standard | 3.66e-02 | 0.00e+00 | 0 | 8.215e-02 | 1.105e-01 |  |
| glass | 16^3 | closed, k=|A| | 3.31e-15 | 5.09e-02 | 0 | 8.538e-02 | 2.343e-02 | 2 / 5 |
| glass | 16^3 | closed, k=1 | 3.15e-15 | 3.41e-02 | 9049 | 8.069e-02 | 2.253e-02 |  |
| glass | 24^3 | standard | 3.72e-02 | 0.00e+00 | 0 | 8.411e-02 | 1.122e-01 |  |
| glass | 24^3 | closed, k=|A| | 3.01e-15 | 5.18e-02 | 0 | 8.763e-02 | 1.361e-02 | 2 / 5 |
| glass | 24^3 | closed, k=1 | 2.30e-15 | 3.47e-02 | 30894 | 8.267e-02 | 1.297e-02 |  |
| glass | 32^3 | standard | 3.71e-02 | 0.00e+00 | 0 | 8.368e-02 | 1.118e-01 |  |
| glass | 32^3 | closed, k=|A| | 3.59e-15 | 5.16e-02 | 0 | 8.727e-02 | 9.780e-03 | 2 / 5 |
| glass | 32^3 | closed, k=1 | 3.08e-15 | 3.46e-02 | 72799 | 8.223e-02 | 9.223e-03 |  |
| jitter | 16^3 | standard | 1.89e-01 | 0.00e+00 | 0 | 2.586e-01 | 6.219e-01 |  |
| jitter | 16^3 | closed, k=|A| | 1.67e-14 | 2.66e-01 | 0 | 2.964e-01 | 5.967e-02 | 2 / 9 |
| jitter | 16^3 | closed, k=1 | 1.11e-14 | 1.59e-01 | 11050 | 2.575e-01 | 5.284e-02 |  |
| jitter | 24^3 | standard | 1.91e-01 | 0.00e+00 | 0 | 2.593e-01 | 6.214e-01 |  |
| jitter | 24^3 | closed, k=|A| | 1.39e-14 | 2.68e-01 | 0 | 2.988e-01 | 3.801e-02 | 2 / 8 |
| jitter | 24^3 | closed, k=1 | 1.23e-14 | 1.60e-01 | 37401 | 2.592e-01 | 3.276e-02 |  |
| jitter | 32^3 | standard | 1.90e-01 | 0.00e+00 | 0 | 2.560e-01 | 6.118e-01 |  |
| jitter | 32^3 | closed, k=|A| | 1.98e-14 | 2.66e-01 | 0 | 2.957e-01 | 2.801e-02 | 2 / 8 |
| jitter | 32^3 | closed, k=1 | 1.62e-14 | 1.59e-01 | 88034 | 2.562e-01 | 2.412e-02 |  |
