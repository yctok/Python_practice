# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 23:33:53 2023

@author: user
"""

import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

psi_solps = [0.56591402, 0.59635553, 0.65365526, 0.70507622, 0.75083496,
       0.79132874, 0.82698182, 0.85806206, 0.88490369, 0.90789432,
       0.92738632, 0.94367313, 0.95706941, 0.96795829, 0.97677538,
       0.9838775 , 0.98955578, 0.99415907, 0.99803803, 1.002408  ,
       1.00753157, 1.01263476, 1.01772166, 1.02279374, 1.02785249,
       1.03288158, 1.03794617, 1.04306613, 1.04817989, 1.05328886,
       1.05838546, 1.06347049, 1.06855367, 1.07363646, 1.07872032,
       1.08380671, 1.08889011, 1.09145489]

dsa = [-0.10817856, -0.09990042, -0.08461333, -0.07136643, -0.05987293,
       -0.04989893, -0.04124898, -0.03377955, -0.02737025, -0.02191068,
       -0.01729219, -0.01343244, -0.01025504, -0.00766848, -0.00557757,
       -0.00389362, -0.00254028, -0.00144284, -0.00052025,  0.00052025,
        0.00174129,  0.00295871,  0.00417348,  0.00538548,  0.0065947 ,
        0.00779984,  0.00900968,  0.01022689,  0.01144393,  0.01266685,
        0.01389157,  0.01511209,  0.01633297,  0.01755502,  0.01877859,
        0.02000397,  0.02123034,  0.02184967]

neexppsi = [0.270082, 0.275369, 0.280657, 0.285944, 0.291232, 0.296519,
       0.301807, 0.307094, 0.312381, 0.317669, 0.322956, 0.328244,
       0.333531, 0.338819, 0.344106, 0.349394, 0.354681, 0.359968,
       0.365256, 0.370543, 0.375831, 0.381118, 0.386406, 0.391693,
       0.39698 , 0.402268, 0.407555, 0.412843, 0.41813 , 0.423418,
       0.428705, 0.433993, 0.43928 , 0.444567, 0.449855, 0.455142,
       0.46043 , 0.465717, 0.471005, 0.476292, 0.481579, 0.486867,
       0.492154, 0.497442, 0.502729, 0.508017, 0.513304, 0.518592,
       0.523879, 0.529166, 0.534454, 0.539741, 0.545029, 0.550316,
       0.555604, 0.560891, 0.566178, 0.571466, 0.576753, 0.582041,
       0.587328, 0.592616, 0.597903, 0.60319 , 0.608478, 0.613765,
       0.619053, 0.62434 , 0.629628, 0.634915, 0.640203, 0.64549 ,
       0.650777, 0.656065, 0.661352, 0.66664 , 0.671927, 0.677215,
       0.682502, 0.687789, 0.693077, 0.698364, 0.703652, 0.708939,
       0.714227, 0.719514, 0.724802, 0.730089, 0.735376, 0.740664,
       0.745951, 0.751239, 0.756526, 0.761814, 0.767101, 0.772388,
       0.777676, 0.782963, 0.788251, 0.793538, 0.798826, 0.804113,
       0.8094  , 0.814688, 0.819975, 0.825263, 0.83055 , 0.835838,
       0.841125, 0.846413, 0.8517  , 0.856987, 0.862275, 0.867562,
       0.87285 , 0.878137, 0.883425, 0.888712, 0.893999, 0.899287,
       0.904574, 0.909862, 0.915149, 0.920437, 0.925724, 0.931012,
       0.936299, 0.941586, 0.946874, 0.952161, 0.957449, 0.962736,
       0.968024, 0.973311, 0.978598, 0.983886, 0.989173, 0.994461,
       0.999748, 1.005036, 1.010323, 1.015611, 1.020898, 1.026185,
       1.031473, 1.03676 , 1.042048, 1.047335, 1.052623, 1.05791 ]

neexp = [4.42268e+19, 4.42268e+19, 4.42268e+19, 4.42268e+19, 4.42267e+19,
       4.42267e+19, 4.42267e+19, 4.42267e+19, 4.42267e+19, 4.42266e+19,
       4.42266e+19, 4.42266e+19, 4.42265e+19, 4.42265e+19, 4.42264e+19,
       4.42263e+19, 4.42263e+19, 4.42262e+19, 4.42260e+19, 4.42259e+19,
       4.42258e+19, 4.42256e+19, 4.42254e+19, 4.42251e+19, 4.42248e+19,
       4.42245e+19, 4.42241e+19, 4.42236e+19, 4.42231e+19, 4.42224e+19,
       4.42217e+19, 4.42208e+19, 4.42198e+19, 4.42186e+19, 4.42172e+19,
       4.42156e+19, 4.42137e+19, 4.42114e+19, 4.42088e+19, 4.42058e+19,
       4.42022e+19, 4.41980e+19, 4.41932e+19, 4.41875e+19, 4.41808e+19,
       4.41731e+19, 4.41640e+19, 4.41534e+19, 4.41410e+19, 4.41266e+19,
       4.41099e+19, 4.40903e+19, 4.40676e+19, 4.40411e+19, 4.40104e+19,
       4.39748e+19, 4.39336e+19, 4.38859e+19, 4.38309e+19, 4.37675e+19,
       4.36947e+19, 4.36114e+19, 4.35163e+19, 4.34081e+19, 4.32857e+19,
       4.31478e+19, 4.29933e+19, 4.28214e+19, 4.26312e+19, 4.24227e+19,
       4.21958e+19, 4.19512e+19, 4.16901e+19, 4.14142e+19, 4.11261e+19,
       4.08286e+19, 4.05250e+19, 4.02190e+19, 3.99144e+19, 3.96147e+19,
       3.93235e+19, 3.90439e+19, 3.87783e+19, 3.85287e+19, 3.82966e+19,
       3.80826e+19, 3.78872e+19, 3.77100e+19, 3.75505e+19, 3.74079e+19,
       3.72811e+19, 3.71689e+19, 3.70701e+19, 3.69835e+19, 3.69078e+19,
       3.68418e+19, 3.67844e+19, 3.67347e+19, 3.66917e+19, 3.66545e+19,
       3.66224e+19, 3.65948e+19, 3.65711e+19, 3.65506e+19, 3.65331e+19,
       3.65180e+19, 3.65051e+19, 3.64941e+19, 3.64846e+19, 3.64764e+19,
       3.64695e+19, 3.64635e+19, 3.64584e+19, 3.64541e+19, 3.64503e+19,
       3.64471e+19, 3.64444e+19, 3.64421e+19, 3.64401e+19, 3.64384e+19,
       3.64369e+19, 3.64356e+19, 3.64346e+19, 3.64337e+19, 3.64329e+19,
       3.64322e+19, 3.64316e+19, 3.64311e+19, 3.64306e+19, 3.64300e+19,
       3.64286e+19, 3.64240e+19, 3.64070e+19, 3.63417e+19, 3.60913e+19,
       3.51609e+19, 3.20949e+19, 2.49862e+19, 1.66486e+19, 1.20915e+19,
       1.05598e+19, 1.01322e+19, 9.97870e+18, 9.26270e+18, 5.41840e+18,
       3.60150e+18, 3.46050e+18, 3.45260e+18, 3.45210e+18, 3.45210e+18]


plt.figure(1)
plt.plot(psi_solps, dsa, color='r')
plt.show()

psi_to_dsa_func = interpolate.interp1d(psi_solps, dsa, fill_value = 'extrapolate')
dsa_neprofile = psi_to_dsa_func(neexppsi)

plt.figure(2)
plt.plot(neexppsi, dsa_neprofile, color='r')
plt.show()

gnexp = np.gradient(neexp) / np.gradient(dsa_neprofile)

gnexp_dsafunc = interpolate.interp1d(dsa_neprofile, gnexp, kind='cubic', fill_value = 'extrapolate')
gnexp_solpslocs = gnexp_dsafunc(dsa)

plt.figure(3)
plt.plot(dsa, gnexp_solpslocs, color='r')
plt.show()

expden_dsa_func = interpolate.interp1d(dsa_neprofile, neexp, fill_value = 'extrapolate')

ne_decay_len_end = (expden_dsa_func(dsa[-12]) - expden_dsa_func(dsa[-11])) / \
    np.mean([expden_dsa_func(dsa[-12]), expden_dsa_func(dsa[-11])])

print(ne_decay_len_end)