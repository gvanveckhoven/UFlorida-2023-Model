# Allison Vollmer 5/22/23
# George VanVeckhoven 9/10/23
# iGEM code for determining interactions between cytokines/cells

import math as mt
import random as rd
import numpy as np
import matplotlib.pyplot as plt

def main():

    # All initial values for parameters and fucntions
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
    t = []
    ittr = 10
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
    
    for i in range(ittr):
        t.append(i)
    
   
    def functions (t):
        
        H = 10
        A = 10
        P = 10
        N = 10
        N_2 = 10
        L_t = 10
        a = 10
        
        H_arr = [H]
        A_arr = [A]
        P_arr = [P]
        N_arr = [N]
        N_2_arr = [N_2]
        L_t_arr = [L_t]
        a_arr = [a]
        output = []
        
        for x in t:
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
            for d in range (1, 5):
                d_n = rd.randint(0, 100)
                total_sum_p += d_n
            mp_bind = total_sum_p / L_b
            
            "Function for ma_bind"
            total_sum_a = 0
            for e in range (1, 5):
                e_n = rd.randint(0, 100)
                total_sum_a += e_n
            ma_bind = total_sum_a / L_a
            
            "Function for dadt"
            dadt = ((L_t * (((P + N)/mp_bind) - (A/ma_bind))) * ((1 - a) - (L_a) * (D - A_l - t_half_leukocytes))) / 2
            
            "Function for m_pk"
            total_sum_pk = 0
            for k in range (1, 5):
                k_n = rd.randint(0, 100)
                total_sum_pk += k_n
            m_pk = total_sum_pk / L_a
            
            "Function for m_ph"
            total_sum_ph = 0
            for h in range (1, 5):
                h_n = rd.randint(0, 100)
                total_sum_ph += h_n
            m_ph = total_sum_ph / H
            
            "Function for d_p"
            d_p = t_half_n * mt.log(0.5) * N_2
            
            "Function for d_m"
            d_m = t_half_m * mt.log(0.5) * M
            
            "Function for dPdt"
            dPdt = (m_pk * L_a) + (m_ph * H) - d_p - (k_p * mp_bind * L_b)
            
            "Function for dN2dt"
            dN2dt = t_half_n * mt.exp(0.5) * N_2
            
            "Function for m_mi"
            total_sum_mi = 0
            for s in range (1, 5):
                s_n = rd.randint(0, 100)
                total_sum_mi += s_n
            m_mi = total_sum_pk / L_b
            
            "Function for m_mh"
            total_sum_mh = 0
            for r in range (1, 5):
                r_n = rd.randint(0, 100)
                total_sum_mi += r_n
            m_mh = total_sum_mh / H
            
            "Function for dAdt"
            dAdt = (m_mi * L_b) + (m_mh * H) - d_m - (k_a * ma_bind * L_a)
            
            "Linear appx's using previously calculated rates"
            H = (dHdt * x) + H
            
            N = (dNdt * x) + N
            
            L_t = (dLdt * x) + L_t
            
            a = (dadt * x) + a
            
            P = (dPdt * x) + P
            
            N_2 = (dN2dt * x) + N_2
            
            A = (dAdt * x) + A
            
            H_arr.append(H)
            N_arr.append(N)
            L_t_arr.append(L_t)
            a_arr.append(a)
            P_arr.append(P)
            N_2_arr.append(N_2)
            A_arr.append(A)
            
            output = [H_arr, N_arr, L_t_arr, a_arr, P_arr, N_2_arr, A_arr]
        return (output)
    
   
    def normalize (fun_arr):
        norm_arr = []
        for x in fun_arr:
            x = x/fun_arr[0]
            norm_arr.append(x)
        return norm_arr
   
    
    output_arr = []
    output_arr_nonnorm = functions(t)
    
    for arr in output_arr_nonnorm:
        output_arr.append(normalize(arr))
    print(output_arr)
    
    t.insert(0,0)
    
    plt.figure(1)
    plt.plot(t, output_arr[0])
    plt.title('HSPC\'s')
    
    plt.figure(2)
    plt.plot(t, output_arr[1])
    plt.title('Instigator Units')
    
    plt.figure(3)
    plt.plot(t, output_arr[2])
    plt.title('Leukocytes')
    
    plt.figure(4)
    plt.plot(t, output_arr[3])
    plt.title('alpha')
    
    plt.figure(5)
    plt.plot(t, output_arr[4])
    plt.title('Pathogen Units')
    
    plt.figure(6)
    plt.plot(t, output_arr[5])
    plt.title('Proinflammatory Units')
    
    plt.figure(7)
    plt.plot(t, output_arr[6])
    plt.title('Antiinflammatory Units')
    
    fig, axs = plt.subplots(7)
   
    axs[0].plot(t, output_arr[0])
    axs[0].set_title('HSPC\'s')
    axs[1].plot(t, output_arr[1])
    axs[1].set_title('Instigator Units')
    axs[2].plot(t, output_arr[2])
    axs[2].set_title('Luekocytes')
    axs[3].plot(t, output_arr[3])
    axs[3].set_title('alpha')
    axs[4].plot(t, output_arr[4])
    axs[4].set_title('Pathogen Units')
    axs[5].plot(t, output_arr[5])
    axs[5].set_title('Proinfalammatory Units')
    axs[6].plot(t, output_arr[6])
    axs[6].set_title('Antiinflammatory Units')
    
    plt.subplots_adjust(hspace = 5)

if __name__ == "__main__":
    main()

