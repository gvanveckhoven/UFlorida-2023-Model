# Allison Vollmer 5/22/23
# iGEM code for determining interactions between cytokines/cells

from classes import Macrophage, GMCSF, IFNg, IL6, Neutrophil, TGFb
import math as mt
import random as rd
import numpy as np
import matplotlib.pyplot as plt

def main():

    #macrophageLevel = input(print("Enter macrophage level: "))

    # create objects of all cytokine signaling attributes
    # I've arbitrarily set all levels to start at 1 but can be changed to different values or allow for user to input
    d_p = 100
    d_m = 100
    g_p = 8
    k_p = 0.1
    k_a = 0.1
    k_np = 100
    t_half_hspc = 100
    w = 0.01
    p_crit = 100
    a_crit = 100
    t_half_leukocytes = 100
    t_half_n = 100
    t_half_m = 100
    gamma = 0.001
    s_a = 100
    t = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    A = 100
    P = 100
    E = 100
    R = 100
    H = 100
    D = 100
    N = 100
    N_inf = 100
    M = 100
    S = 100
    a = 100
    L_t = 100
    L_a = 100
    L_b = 100
    A_l = 100
    m_pk = 100
    m_ph = 100
    m_mi = 100
    m_mh = 100
    mp_bind = 100
    ma_bind = 100
    dLdt = 100
    dadt = 100
    dHdt = 100
    dNdt = 100
    dN2dt = 100
    dPdt = 100
    
    

    "For HSPC's"
    
    "For E(t)"
   
    def functions (t):
        
        "Function for E"
        if (P - (s_a*A)<p_crit):
            E = 1 + (P/p_crit) - s_a* (A/a_crit)
        elif (P - (s_a*A)>=p_crit):
                E = 2 / (w* (P-p_crit) + 1)
        
        "Function for R"
        if (P-A < p_crit):
            R = (0.05 + 0.2 * ((P/p_crit) - s_a * (A/a_crit))) * H
        elif ((P-A) >= p_crit):
            R = 0.25 * H
            
        "Function for D"
        if (P-A < p_crit):
            D = (0.1 + 0.5 * ((P/p_crit) - s_a * (A/a_crit))) * H
        elif ((P-A) >= p_crit):
            D = 0.6 * H
        
        "Function for dHdt"
        dHdt = (E*R) - D
        
        "Function for S"
        if (N <= 0):
            S = 0
        elif(N > 0):
            S = 1
        
        "Function for dNdt"
        dNdt = (g_p * N) * (N/N_inf) - (k_np * L_a) * S
            
        "Function for L_t"
        L_t = L_a + L_b
        
        "Function for A_l"
        A_l = (0.1 * N) * (1.5/((1+mt.exp(gamma * L_t)) - 0.75))
        
        "Function for dLdt"
        dLdt = D - A_l - (t_half_leukocytes * mt.log(0.5) * L_t)
        
        "Function for alpha (a)"
        a = L_a / L_t
        
        "Function for mp_bind"
        total_sum_p = 0
        for d in range (1, L_b + 1):
            d = rd.randint(0, 100)
            total_sum_p += d
        mp_bind = total_sum_p / L_b
        
        "Function for ma_bind"
        total_sum_a = 0
        for e in range (1, L_a + 1):
            e = rd.randint(0, 100)
            total_sum_a += e
        ma_bind = total_sum_a / L_a
        
        "Function for dadt"
        dadt = ((L_t * (((P + N)/mp_bind) - (A/ma_bind))) * ((1 - a) - (L_a) * (D - A_l - t_half_leukocytes))) / 2
        
        "Function for m_pk"
        total_sum_pk = 0
        for k in range (1, L_a + 1):
            k = rd.randint(0, 100)
            total_sum_pk += k
        m_pk = total_sum_pk / L_a
        
        "Function for m_ph"
        total_sum_ph = 0
        for h in range (1, H + 1):
            h = rd.randint(0, 100)
            total_sum_ph += h
        m_ph = total_sum_ph / H
        
        "Function for d_p"
        d_p = t_half_n * mt.log(0.5) * N
        
        "Function for d_m"
        d_m = t_half_m * mt.log(0.5) * M
        
        "Function for dPdt"
        dPdt = (m_pk * L_a) + (m_ph * H) - d_p - (k_p * mp_bind * L_b)
        
        "Function for dN2dt"
        dN2dt = t_half_n * mt.exp(0.5) * N
        
        "Function for m_mi"
        total_sum_mi = 0
        for s in range (1, L_b + 1):
            s = rd.randint(0, 100)
            total_sum_mi += s
        m_mi = total_sum_pk / L_b
        
        "Function for m_mh"
        total_sum_mh = 0
        for r in range (1, H + 1):
            r = rd.randint(0, 100)
            total_sum_mi += r
        m_mh = total_sum_mh / H
        
        "Function for dAdt"
        dAdt = (m_mi * L_b) + (m_mh * H) - d_m - (k_a * ma_bind * L_a)
        
        return ([dHdt, dNdt, dLdt, dadt, dPdt, dN2dt, dAdt])
   
    arr = functions(t)
   
    plt.plot(t, arr)
  

if __name__ == "__main__":
    main()

