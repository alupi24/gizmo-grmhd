# Face closure: standard MFM faces versus HYDRO_TENSOR_FACE_CORRECTIONS

closure = rms_i |sum_j A_ij| / sum_j |A_ij|. `own-side` = the part of the standard closure coming from particle i's own kernel moments, -V_i E_i^{-1} S1_i (the term the tensor correction is designed to remove).

| set | N | standard | own-side part | rest | tensor (as coded, S2 - S3.B) | tensor (S2 + S3.B) | first-moment check | faces with A.dp<0 (coded / flipped) |
|---|---|---|---|---|---|---|---|---|
| lattice | 16^3 | 2.607e-17 | 1.053e-17 | 2.015e-17 | 2.605e-17 | 2.605e-17 | 1.0e+00 | 0 / 0 |
| lattice | 24^3 | 9.836e-16 | 6.616e-16 | 3.547e-16 | 5.828e-17 | 5.828e-17 | 1.0e+00 | 0 / 0 |
| glass | 16^3 | 3.656e-02 | 1.861e-02 | 1.890e-02 | 1.172e-02 | 1.189e-02 | 5.8e-14 | 0 / 0 |
| glass | 24^3 | 3.722e-02 | 1.891e-02 | 1.929e-02 | 1.203e-02 | 1.224e-02 | 3.6e-14 | 0 / 0 |
| jitter | 16^3 | 1.924e-01 | 9.393e-02 | 1.024e-01 | 6.243e-02 | 6.176e-02 | 1.1e-14 | 379 / 453 |
| jitter | 24^3 | 1.881e-01 | 9.176e-02 | 1.005e-01 | 6.097e-02 | 6.145e-02 | 1.4e-14 | 1068 / 1298 |
