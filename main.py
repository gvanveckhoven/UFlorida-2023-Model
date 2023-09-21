# Allison Vollmer 5/22/23
# George VanVeckhoven 9/10/23
# iGEM code for determining interactions between cytokines/cells

import math as mt
import random as rd
import numpy as np
import matplotlib.pyplot as plt

def main():

    # All initial values for parameters and fucntions
   # S_n = 100
   # A = 100
   # P = 100
   # E = 100
   # H = 100
   # D = 100
   # N = 100
   # a = 100
   # T = 100
   # R = 100
   # S = 100
   # A_l = 100
   # mp_bind = 100
   # ma_bind = 100
   # dLdt = 100
   # dadt = 100
   # dHdt = 100
   # dNdt = 100
   # dN2dt = 100
   # dPdt = 100
    
    k_p = 0.5
    k_a = 0.1
    k_np = 0.6
    w = 0.01
    p_crit = 100
    a_crit = 100
    t_half_leukocytes = 7
    t_half_n = 4.1375
    t_half_m = 1
    t_double_P = 15
    gamma = 0.001
    S_a = 1
    t = []
    ittr = 100
    N_inf = 20000000
    M = 1
    S_pr = 1
    S_ph = 1
    S_as = 5
    S_ah = 5
    theta_ps = 1
    theta_ar = 1

    
    for i in range(ittr):
        i += 1
        t.append(float(i) / 10)
    
    def functions (t):
        
        H = 1000
        A = 250
        P = 500
        N = 5000
        T = 7000
        a = 0.1
       
        S = 10
        N_2 = 10
        
        
        H_arr = [H]
        A_arr = [A]
        P_arr = [P]
        N_arr = [N]
        N_2_arr = [N_2]
        T_arr = [T]
        a_arr = [a]
        output = []
        
        for x in t:
            "Test for overflow at the start"
            try:
                A_l = (0.1 * N) * (1.5/((1 + mt.exp(gamma * T)) - 0.75))
            except OverflowError:
                break
            
            "Function for S"
            if (N <= 0):
                S_n = 0
            elif(N > 0):
                S_n = 1
            
            "Function for E"
            if ((P + S_n * N - S_a * A) < p_crit):
                E = 1
            elif ((P + S_n * N - S_a * a) >= p_crit):
                E = 2 - (2/(1+mt.exp(-1 * w * (P + S_n * N - S_a * A - p_crit))))
            
            if E < 0:
                E = 0
                
            
            "Function for R"
            if ((P + S_n * N - S_a * A) < p_crit):
                R = (0.05 + 0.25 * (P + S_n * N - S_a * A) * H)
            elif ((P + S_n * N - S_a * A * -1) >= p_crit):
                R = 0.25 * H
                
            if R < 0:
                R = 0
                
            "Function for D"
            if ((P + S_n * N - S_a * A) < p_crit):
                D = (0.1 + 0.5 * ((P/p_crit) - S_a * (A/a_crit))) * H
            elif ((P + S_n * N - S_a * A * -1) >= p_crit):
                D = 0.6 * H
                
            if D < 0:
                D = 0
            
            "Function for dHdt"
            dHdt = (E*R) - D
            
            "Function for dNdt"
            dNdt = ((mt.log(2)/t_double_P)*N) * (1-(N/N_inf)) - (k_np * R * S)
                
            "Function for T"
            T = R + S
            
            if T < 0:
                T = 0
            
           # "Function for A_l"
           # A_l = (0.1 * N) * (1.5/((1 + mt.exp(gamma * T)) - 0.75))
            
            
            "Function for dLdt"
            dTdt = D - A_l - (t_half_leukocytes * mt.log(0.5) * T)
            
            "Function for alpha (a)"
            a = R / T
            
            if a > 1:
                a = 1
            elif a < 0:
                a = 0
            
            #"Function for S"
            #S = L_a - (P + N)
            
            "Function for mp_bind""This needs to be updated"
            mp_bind = ((P + S_n * N) * S) / (P + theta_ps * S)
            
            if mp_bind < 0:
                mp_bind = 0
            
            "Function for ma_bind"
            ma_bind = (A * R) / (theta_ar + A * R)
            
            if ma_bind < 0:
                ma_bind = 0
            
            "Function for dadt"
            try:
                dadt = ((T * ((mp_bind * S) - (ma_bind * R)) * ((1 - a) - (R) * (D - A_l - t_half_leukocytes))) / mt.pow(T, 2))
            except OverflowError:
                dadt = float('inf')
            
            "Function for d_p"
            d_p = t_half_n * mt.log(0.5) * N_2
            
            "Function for d_m"
            d_m = t_half_m * mt.log(0.5) * M
            
            "Function for dPdt"
            dPdt = (S_pr * R) + (S_ph * H) - d_p - (k_p * mp_bind)
            
            "Function for dN2dt"
            dN2dt = t_half_n * mt.exp(0.5) * N_2
            
            "Function for dAdt"
            dAdt = (S_as * S) + (S_ah * H) - d_m - (k_a * ma_bind)
            
            "Linear appx's using previously calculated rates"
            H = (dHdt * x) + H
            
            if H < 0:
                H = 0
            
            N = (dNdt * x) + N
            
            if N < 0:
                N = 0
            
            T = (dTdt * x) + T
            
            if T < 0:
                T = 0
            
            a = (dadt * x) + a
            
            if a > 1:
                a = 1
            elif a < 0:
                a = 0
            
            P = (dPdt * x) + P
            
            if P < 0:
                P = 0
            
            N_2 = (dN2dt * x) + N_2
            
            if N_2 < 0:
                N_2 = 0
            
            A = (dAdt * x) + A
            
            if A < 0:
                A = 0
            
            H_arr.append(H)
            N_arr.append(N)
            T_arr.append(T)
            a_arr.append(a)
            P_arr.append(P)
            N_2_arr.append(N_2)
            A_arr.append(A)
            
            output = [H_arr, N_arr, T_arr, a_arr, P_arr, N_2_arr, A_arr]
        return (output)
    
   
    def normalize (fun_arr):
        norm_arr = []
        for x in fun_arr:
            x = x/fun_arr[0]
            norm_arr.append(x)
        return norm_arr
   
    
    output_arr = []
    output_arr_nonnorm = functions(t)
    
    "To fit the t list to the size of the values inserted incase it breaks"
    n = len(t)
    k = len(output_arr_nonnorm[0])
    for i in range(0, n - k + 1):
        t.pop()
    
    
    for arr in output_arr_nonnorm:
        output_arr.append(normalize(arr))

    
    t.insert(0,0)
    
    plt.figure(1)
    plt.plot(t, output_arr[0])
    plt.title('HSPC\'s')
    plt.savefig("HSPC's_ittr2")
    
    plt.figure(2)
    plt.plot(t, output_arr[1])
    plt.title('Instigator Units')
    plt.savefig("Instagators_ittr2")
    
    plt.figure(3)
    plt.plot(t, output_arr[2])
    plt.title('Leukocytes')
    plt.savefig("Leukocytes_ittr2")
    
    plt.figure(4)
    plt.plot(t, output_arr[3])
    plt.title('alpha')
    plt.savefig("alpha_ittr2")
    
    plt.figure(5)
    plt.plot(t, output_arr[4])
    plt.title('Proinflammatory units')
    plt.savefig("Proinflammatory_ittr2")
    
    plt.figure(6)
    plt.plot(t, output_arr[5])
    plt.title('? Units')
    plt.savefig("?_ittr2")
    
    plt.figure(7)
    plt.plot(t, output_arr[6])
    plt.title('Antiinflammatory Units')
    plt.savefig("Antiinflammatory_ittr2")
    
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
    
    plt.savefig("Graphs_of_ittr_2.jpg")

if __name__ == "__main__":
    main()

