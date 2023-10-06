import numpy as np
import matplotlib.pyplot as plt

from function import simulate

H_0 = 1000
N_0 = 0
P_0 = 0
A_0 = 0
T_0 = 6000
a_0 = 0.00
b_0 = 1
e_0 = 0.00


INIT_VALUES = [H_0, N_0, T_0, a_0, b_0, e_0, P_0, A_0]

g_N = 0.5
K_PS = 0.7
K_AQ = 0.7
K_AS = 0.7
K_PU = 0.7
k_nq = 0.85
k_ns = 0.2
k_tn = 0.33
w = 0.0005
P_crit = 2500
A_crit = 1000
t_half_leukocytes = 7
t_half_P = 4.1375
t_half_A = 7
t_double_P = 15
gamma = 0.0000001
S_a = 1
S_n = 1
t = np.linspace(0, 200, 200)    #timeseries set from 0 hrs -> 100 hrs w/ delta = 1 hr
N_inf = 20000000
S_PQ = 0.33
S_PH = 0.01
S_PS = 0.02
S_AS = 0.04
S_AH = 0.01
S_AU = 0.33
theta_ps = 10_000_000
theta_ar = 10_000_000
theta_AS = 10_000_000
theta_UP = 10_000_000
Immune_start = 500
Active_start = 500
Immune_crit = 2500
Active_crit = 2500
y = 0.0008

RATES = [g_N, K_PS, K_AQ, K_AS, K_PU, k_nq, k_ns, k_tn, w, P_crit, A_crit, S_a, S_n, N_inf, S_PQ, S_PH, S_PS, S_AS, S_AH, S_AU, 
         theta_ps, theta_ar, theta_AS, theta_UP, Immune_start, Active_start, Immune_crit, Active_crit, y]

output = simulate(INIT_VALUES, RATES, 200)

# ---------------- plots: -------------------

plt.figure(1)
plt.plot(t, output[0])
plt.plot(t, [1000]*200, "k--")
plt.title('HSPC\'s')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.ylim(0, 2000)

plt.figure(2)
plt.plot(t, output[1])
plt.plot(t, [N_0]*200, "k--")
plt.title('Pathogens')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(3)
plt.plot(t, output[2])
plt.plot(t, np.ones(200), "k--")
plt.title('Pro-inflammatory cytokines')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(4)
plt.plot(t, output[3])
plt.title('Anti-inflammatory cytokines')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

'''plt.figure(4)
plt.plot(t, output[5], "r-")
plt.plot(t, output[7], "b-")
plt.plot(t, output[4], "g-")
plt.plot(t, output[6], "k-")
plt.title("Active / Total / Stable / Immuno Leukocytes")
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
plt.legend(["Active", "Total", "Stable", "Immuno-"], loc ="lower right")'''

plt.figure(5)
plt.plot(t, output[4])
plt.title('Stable Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(6)
plt.plot(t, output[5])
plt.title('Active Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(7)
plt.plot(t, output[6])
plt.title('Immuno-suppressive Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")

plt.figure(8)
plt.plot(t, output[7])
plt.title('Total Leukocytes')
plt.xlabel("Time (hrs)")
plt.ylabel("Density (per cubic centimeter)")
