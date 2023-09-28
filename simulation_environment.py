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
   # mu_SP = 100
   # mu_QA = 100
   # dLdt = 100
   # dAdt = 100
   # dHdt = 100
   # dndt = 100
   # dN2dt = 100
   # dPdt = 100
    
    g_N = 0.3
    k_p = 0.5
    k_a = 1
    k_np = 0.6
    k_tn = 0.5
    w = 0.00001
    p_crit = 2500
    A_crit = 1000
    t_half_leukocytes = 7
    t_half_P = 4.1375
    t_half_A = 4.1375
    t_double_P = 15
    gamma = 0.00000001
    S_a = 1
    S_n = 1
    t = np.linspace(0, 100, 100)    #timeseries set from 0 hrs -> 100 hrs w/ delta = 1 hr
    #ittr = 100
    N_inf = 20000000
    M = 1
    S_PQ = 5
    S_PH = 3
    S_PS = 1
    S_AS = 2
    S_AH = 1
    S_AU = 5
    theta_ps = 10_000_000
    theta_ar = 10_000_000
    theta_AS = 10_000_000
    theta_UP = 10_000_000
    delta_t = 1

    
    '''for i in range(ittr):
        i += 1
        t.append(float(i) / 10)'''
        
    def functions (t):
        #initial values:
            
        H = 1000    #HSPCs
        A = 500    #Anti-inflammatory cytokines
        P = 500    #Pro-inflammatory cytokines
        N = 0    #Pathogens
        T = 7000    #Total Leukocytes
        a = 0.05      #alpha
        b = 0.9       #beta
        e = 0.05       #epsilon
        Q = a*T    #Active leukocytes
        S = e*T    #Stable leukocytes
        U = b*T    #Immuno-suppressive leukocytes
        N_2 = 10    #???
        
        
        H_arr = [H]
        A_arr = [A]
        P_arr = [P]
        N_arr = [N]
        N_2_arr = [N_2]
        T_arr = [T]
        Q_arr = [Q]
        S_arr = [S]
        U_arr = [U]
        output = []
        
        H_0 = H
        P_0 = P
        A_0 = A
        
        # following arrays are for testing
        E_arr = []
        R_arr = []
        D_arr = []
        A_l_arr = []
        T_decay = []
        d_a_arr = []
        
        for x in t[1:]:

            
            '''"Function for S (Switch - dN/dt)"
            if (N <= 0):
                S_n = 0
            elif(N > 0):
                S_n = 1'''
            
            "Function for E (exhaustion - dH/dt)"
            if ((P + S_n * N - S_a * A) < p_crit):
                E = 1
            elif ((P + S_n * N - S_a * a) >= p_crit):
                E = 2 - (2/(1+mt.exp(-1 * w * (P + S_n * N - S_a * A - p_crit))))
            
            if E < 0:
                E = 0
                
            E_arr.append(E)
                
            
            "Function for R (self renewal - dH/dt)"
            if ((P + S_n * N - S_a * A) < p_crit):
                R = (0.05 + 0.25 * ((P + S_n * N - S_a * A)/p_crit)) * H
            elif ((P + S_n * N - S_a * A) >= p_crit):
                R = 0.25 * H
                
            if R < 0:
                R = 0
            
            R_arr.append(R)
            
            "Function for D (differentiation - dH/dt & dT/dt"
            if (P + S_n * N - S_a * A <= 0):
                D = 0.1*H
            elif ((P + S_n * N - S_a * A) <= p_crit):
                D = (0.1 + 0.5*((P + S_n * N - S_a * A)/p_crit))*H
            elif ((P + S_n * N - S_a * a) > p_crit):
                D = 0.6*H
            
            "Function for R_plus"
            
            if S_a*A - S_n*N - P < 0:
                R_plus = 0
            elif S_a*A - S_n*N - P < A_crit:
                R_plus = 0.05*((S_a*A - S_n*N - P)/A_crit)*H
            else:
                R_plus = 0.05*H
            
            
            D_arr.append(D)    
            
            "Function for dHdt"
            dhdt = (E*R + E*R_plus) - D
            print("D: " + str(D))
            print("E: " + str(E))
            print("R: " + str(R))
            print("R_plus: " + str(R_plus))
            
            "Function for dndt"
            dndt = (g_N*N) * (1-(N/N_inf)) - (k_np * Q)
            
            #print((mt.log(1/2)/t_double_P)*N)
                
            '''"Function for T"
            T = Q + S
            
            if T < 0:
                T = 0'''
            
            "Function for A_l"
            try:
                A_l = (k_tn * N) * ((-2/(1 + mt.exp(gamma * T))) + 1)
            except OverflowError:               
                A_l = k_tn * N
                
            A_l_arr.append(A_l)
            

            "Function for mu_SP"
            
            mu_SP = 0.45*(np.power(P*S, 0.75)) / (np.power(P*S, 0.75) + np.power(theta_ps, 0.75))
            
            if mu_SP < 0:
                mu_SP = 0
            
            "Function for mu_QA"
            #mu_QA = (A * Q) / (theta_ar + (A * Q))
            mu_QA = (np.power(A*Q, 0.75)) / (np.power(A*Q, 0.75) + np.power(theta_ar, 0.75))
            if mu_QA < 0:
                mu_QA = 0
                
            "Function for mu_SA"
            mu_SA = 0.45*(np.power(A*S, 0.75)) / (np.power(A*S, 0.75) + np.power(theta_AS, 0.75))
            
            "Function for mu_UP"
            #mu_UP = (U*(P+N))/(U*(P+N) + theta_UP)
            mu_UP = (np.power(U*(P+N), 0.75)) / (np.power(U*(P+N), 0.75) + np.power(theta_UP, 0.75))

            "Function for d_p"
            #d_p = np.log(1/2)*P_0*np.power((1/2), (x/t_half_P))*t_half_P
            d_p = 0.01 * P
            
            "Function for d_a"
            #d_a = A_0 * np.log(0.5) * (1 / t_half_A) * np.power(0.5, (x / t_half_A)) * x
            
            d_a = 0.01 * A
            d_a_arr.append(d_a)
            
            "Function for dN2dt"
            dN2dt = t_half_P * mt.exp(0.5) * N_2
            
            "Function for dPdt"
            dPdt = (S_PS * S) + (S_PQ * Q) + (S_PH * H) + d_p - (k_p * mu_SP)
            
            "Function for dAdt"
            dAdt = (S_AU * U) + (S_AS * S) + (S_AH * H) + d_a - (k_a * mu_QA)
            #print(dAdt)
            
            "Function for dS/dt"    # Stable Leukocytes
            
            dsdt = D - A_l - 0.01*S + mu_QA*Q + mu_UP*U - mu_SA*S - mu_SP*S
            
            
            T_decay.append(mt.log(1/2)*H_0*np.power((1/2), (x/t_half_leukocytes))*t_half_leukocytes)
            
            "Function for dq/dt"    # Active Leukocytes
            dqdt = (mu_SP*S - mu_QA*Q)
            
            "Function for du/dt"    # Immuno-suppressive Leukocytes
            dudt = mu_SA*S - mu_UP*U
            
            
            
            "Linear appx's using previously calculated rates"
            H = dhdt + H
            
            if H < 0:
                H = 0
            
            N = dndt + N
            
            if N < 0:
                N = 0
            
            S = dsdt + S
            
            if S < 0:
                S = 0
            
            Q = dqdt + Q
            U = dudt + U
            
            T = U + S + Q
            a = Q/T
            b = S/T
            e = U/T
            
            
            if a > 1:
                a = 1
            elif a < 0:
                a = 0
            
            P = dPdt + P
            
            if P < 0:
                P = 0
            
            N_2 = (dN2dt * x) + N_2
            
            if N_2 < 0:
                N_2 = 0
            
            A = dAdt + A
            
            if A < 0:
                A = 0
                
            H_arr.append(H)
            N_arr.append(N)
            T_arr.append(T)
            Q_arr.append(Q)
            S_arr.append(S)
            U_arr.append(U)
            P_arr.append(P)
            N_2_arr.append(N_2)
            A_arr.append(A)
            
            output = [H_arr, N_arr, T_arr, Q_arr, S_arr, P_arr, N_2_arr, A_arr, R_arr, D_arr, E_arr, A_l_arr, T_decay, d_a_arr, U_arr]
        return (output)
    
   
    def normalize (fun_arr):
        norm_arr = []
        for x in fun_arr:
            if fun_arr[0] == 0:
                y = 0
                norm_arr.append(y)
            else:
                y = x/fun_arr[0]
                norm_arr.append(y)
        return norm_arr
   
    
    output_arr = []
    output_arr_nonnorm = functions(t)
    
    "To fit the t list to the size of the values inserted incase it breaks"
    '''n = len(t)
    k = len(output_arr_nonnorm[0])
    for i in range(0, n - k + 1):
        t.pop()'''
    
    
    for arr in output_arr_nonnorm:
        output_arr.append(normalize(arr))

    
    #t.insert(0,0)
    
    plt.figure(1)
    plt.plot(t, output_arr_nonnorm[0])
    plt.title('HSPC\'s')
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    plt.ylim(0, 2000)
    #plt.savefig("HSPC's_ittr2")
    
    plt.figure(2)
    plt.plot(t, output_arr_nonnorm[1])
    plt.title('Pathogens')
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    #plt.savefig("Instagators_ittr2")
    
    plt.figure(3)
    plt.plot(t, output_arr[2])
    plt.title('Leukocytes')
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    #plt.savefig("Leukocytes_ittr2")
    
    plt.figure(4)
    plt.plot(t, output_arr_nonnorm[3], "r-")      # don't normalize these
    plt.plot(t, output_arr_nonnorm[2], "g-")
    plt.plot(t, output_arr_nonnorm[4], "b-")
    plt.plot(t, output_arr_nonnorm[14], "k-")
    plt.title("Active / Total / Stable / Immuno Leukocytes")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    plt.legend(["Active", "Total", "Stable", "Immuno-"], loc ="lower right")
    #plt.savefig("alpha_ittr2")
    
    plt.figure(5)
    plt.plot(t, output_arr[4])
    plt.title('Pro-inflammatory Cytokines')
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    #plt.savefig("Proinflammatory_ittr2")
    
    plt.figure(6)
    plt.plot(t, output_arr[5])
    plt.title('? Units')
    #plt.savefig("?_ittr2")
    
    plt.figure(7)
    plt.plot(t, output_arr[7])
    plt.title('Anti-inflammatory cytokines')
    plt.xlabel("Time (hrs)")
    plt.ylabel("Density (per cubic centimeter)")
    #plt.savefig("Antiinflammatory_ittr2")
    
    plt.figure(8)
    plt.plot(t[1:], output_arr[11])
    plt.title('A_l')
    plt.xlabel("Time (hrs)")
    plt.ylabel("A_l")
    
    plt.figure(9)
    plt.plot(t[1:], output_arr_nonnorm[13])
    plt.title('Anti-inflammation decay')
    plt.xlabel("Time (hrs)")
    plt.ylabel("d_a")
    '''
    #debugging code below
    print(len(t))
    print(len(output_arr[0]))
    print(len(output_arr_nonnorm[0]))
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
    plt.show()
    #plt.savefig("Graphs_of_ittr_2.jpg")'''
    
    print("HSPC[1]: " + str(output_arr_nonnorm[0][3]))
    print("Leukocyte[2]: " + str(output_arr_nonnorm[2][2]))
    print("R[1]: " + str(output_arr_nonnorm[7][1]))
    print("D[1]: " + str(output_arr_nonnorm[8][1]))
    print("A_l[1]: " + str(output_arr_nonnorm[10][1]))
    print("Leuko_decay[1]: " + str(output_arr_nonnorm[11][1]))
    #print(output_arr_nonnorm[6])

if __name__ == "__main__":
    main()
