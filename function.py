# -*- coding: utf-8 -*-
"""
This is the main function that calculates the values for each individual function from the derivatives using linear approximation
"""
import numpy as np

def simulate(init_values, rates, t_count):
    '''

    Parameters
    ----------
    init_values : Initial values for each function; f(T=0)
    
    rates : parameter values 
    
    t_count : Number of timeframes; as t_count -> 0, runtime and accuracy increase

    Returns
    -------
    output = ndarray of shape (#, t_count); output[0] = ndarray of timesteps, 
    output[i] denotes the py.list containing the values of function f_i for each t in output[0]

    '''
    #--------- 1. Initiliazing all values and parameters
    
    timesteps = np.linspace(0, t_count, t_count)
    
    H_0 = init_values[0]        # initial HSPC value
    N_0 = init_values[1]        # initial pathogen value
    T_0 = init_values[2]        # initial combined value of all leukocytes
    a_0 = init_values[3]        # initial alpha value
    b_0 = init_values[4]        # initial beta value
    e_0 = init_values[5]        # initial epsilon value
    P_0 = init_values[6]        # initial pro-inflammatory cytokines value
    A_0 = init_values[7]        # initial anti-inflammatory cytokines value
    
    H = H_0
    N = N_0
    T = T_0
    a = a_0
    b = b_0
    e = e_0
    P = P_0
    A = A_0
    
    Q = a * T
    S = b * T
    U = e * T
    
    H_output = [H]
    N_output = [N]
    T_output = [T]
    Q_output = [Q]
    S_output = [S]
    U_output = [U]
    P_output = [P]
    A_output = [A]
    T_output = [T]
    
    g_N = rates[0]          # pathogen growth rate per hour
    K_PS = rates[1]         #
    K_AQ = rates[2]
    K_AS = rates[3]
    K_PU = rates[4]
    k_nq = rates[5]
    k_ns = rates[6]
    k_tn = rates[7]
    w = rates[8]
    P_crit = rates[9]
    A_crit = rates[10]
    S_a = rates[11]
    S_n = rates[12]
    N_inf = rates[13]
    S_PQ = rates[14]
    S_PH = rates[15]
    S_PS = rates[16]
    S_AS = rates[17]
    S_AH = rates[18]
    S_AU = rates[19]
    theta_ps = rates[20]
    theta_ar = rates[21]
    theta_AS = rates[22]
    theta_UP = rates[23]
    Immune_start = rates[24]
    Active_start = rates[25]
    Immune_crit = rates[26]
    Active_crit = rates[27]
    y = rates[28]
    count = 1
    
    #----------- 2. Calculating derivative values for each t > t_0
    
    for x in timesteps[1:]:
        
            
        #---------- mechanism for repeating pathogenic insult ---------------------
        
        if count == 50 or count == 75:
            N += 1800
            count += 1
        else:
            count += 1
        
        #---------- 2a. Calculating individual terms ---------------
        "Function for E (exhaustion - dH/dt)"
        
        if (P + S_n * N - S_a * A) < P_crit:
            E = 1
        elif (P + S_n * N - S_a * a) >= P_crit:
            try:
                E = 2 - (2/(1+np.exp(-1 * y * (P + S_n * N - S_a * A - P_crit))))
            except FloatingPointError:
                E = 0
            
        if E < 0:
            E = 0
            
        "Function for R (self renewal - dH/dt)"
        
        if P + S_n * N - S_a * A <= 0:
            R = 0.1 * H
        elif P + S_n * N - S_a * A <= P_crit:
            R = (0.1 + 0.35 * ((P + S_n * N - S_a * A)/P_crit)) * H
        elif P + S_n * N - S_a * A > P_crit:
            R = 0.45 * H
            
        "Function for D (differentiation - dH/dt & dS/dt"
        
        if (P + S_n * N - S_a * A <= 0):
            D = 0.1*H
        elif ((P + S_n * N - S_a * A) <= P_crit):
            D = (0.1 + 0.35*((P + S_n * N - S_a * A)/P_crit))*H
        elif ((P + S_n * N - S_a * a) > P_crit):
            D = 0.45*H
            
        "BM Niche Renewal contribution (dH/dt)"
        
        if H_0 >= H:
            H_stable = 0.05 * H * (1-(H/H_0))
        
        else:
            H_stable = 0
            
        "Function for A_l (dS/dt - stable leukocyte kill rate by pathogens)"
        
        try:
            A_l = (k_tn * N) * ((-2/(1 + np.exp(w * S))) + 1)
        except FloatingPointError:
            A_l = k_tn * N

        "Function for mu_SP"
        
        if P + S_n*N - S_a*A <= Active_start:
            mu_SP = 0
            
        elif P + S_n*N - S_a*A <= Active_crit:
            mu_SP = 0.2 * ((P + N) / Active_crit) * S
            
        else:
            mu_SP = 0.2 * S
            
        "Function for mu_QA"

        mu_QA = (A*Q) / (theta_ar + A*Q)    # as of 10/4, not being used

            
        "Function for mu_SA"
        
        if P + N <= Immune_start:
            mu_SA = 0
            
        elif P + N <= Immune_crit:
            mu_SA = 0.2 * ((P + N) / Immune_crit) * S
            
        else:
            mu_SA = 0.2 * S
        
        "Function for mu_UP"

        mu_UP = (U*(P+N)) / (theta_UP + U*(P+N))    # as of 10/4, not being used
        
        "Function for D_P (dP/dt - pro-cytokine decay term"

        D_P = 0.25 * P
        
        "Function for D_A (dA/dt - anti-cytokine decay term"
        
        D_A = 0.25 * A
        
        "Function for D_S (S decay)"
        
        D_S = (1/70) * S        # chosen to balance out +D in steady state
        
        "Function for D_Q (Q decay)"
        
        D_Q = 0.25 * Q
        
        "Function for D_U (U decay)"
        
        D_U = 0.05 * U
        
        #---------- 2b. Calculating derivatives ---------------------
        
        dHdt = (E*R) + H_stable - D                                     # HSPC derivative
        
        dNdt = (g_N*N) * (1-(N/N_inf)) - (k_nq * Q) - (k_ns * S)        # Pathogen derivative
        
        dPdt = (S_PS * S) + (S_PQ * Q) + (S_PH * H) - D_P     # Pro-inflammatory derivative
        
        dAdt = (S_AU * U) + (S_AS * S) + (S_AH * H) - D_A     # Anti-inflammatory derivative
        
        dSdt = D - A_l - D_S - mu_SA - mu_SP                            # Stable Leukocytes derivative
        
        dQdt = mu_SP - D_Q                                  # Active Leukocytes derivative
        
        dUdt = mu_SA - D_U                                  # Immuno-suppressive derivative
        
        #--------- 2c. Diagnostics --------------------------------
        
        #fill with code for checking variable values as needed
        
        #--------- 3. Updating lists and functions with linear approximation -------------
        
        H = dHdt + H        # control flow may affect function values
        if H < 0:
            H = 0
        
        N = dNdt + N
        if N < 0:
            N = 0
        
        S = dSdt + S
        if S < 0:
            S = 0
        
        Q = dQdt + Q
        if Q < 0:
            Q = 0
        
        U = dUdt + U
        if U < 0:
            U = 0
                
        P = dPdt + P
        if P < 0:
            P = 0
        
        A = dAdt + A
        if A < 0:
            A = 0
        
        T = Q + U + S
        
        H_output.append(H)
        N_output.append(N)
        P_output.append(P)
        A_output.append(A)
        S_output.append(S)
        Q_output.append(Q)
        U_output.append(U)
        T_output.append(T)
        
    #------------- 4. Output --------------
    
    output = np.array([H_output, N_output, P_output, A_output, S_output, Q_output, U_output, T_output])
    
    return output
    
        
        
        
    
