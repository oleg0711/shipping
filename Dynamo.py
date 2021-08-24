# -*- coding: utf-8 -*-
"""
Created on Sat May  8 11:13:36 2021

"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize

# TODO: 1. Get real values from the paper Ge/Beullens 2021
# TODO: 2. Align formulas in the model paper
# TODO: 3. Run the algorithm



# O.U. parameters of the log rate model:
# the data is assuming T=1 it is one year.
r_bar = 8 #  e^8 ~ 3000  
r_lambda = 0.5
r_sigma = 0.1

alpha = 0.08 / 365 # decay factor of the system, internal rate of return of the shipping business or interest rates at which this business is financed
eta = 0.1 # risk aversion of the shipping agent

fut_curve_slope = -0.1 # negative means that there is a discount to the expectation built-in by the model
#r_state is the difference between the prevailing log shipping rate and its long-term average 
v_max = 17
v_min = 10

m = 2
S = [8000,80000] 
w = [0.3,0.3]

f_TCH = 20000

C_u_j = 0
C_h_j = 0
c_j_f = 1
# formula 5 Ge/Beullens

G0 = 5000000 # future profit potential

def gen_states(r,r_sigma):
    states,states_step = np.linspace(start=r-r_sigma*3,
                                 stop=r_bar+r_sigma*3,
                                 num=10)
    return states,states_step


def fuel_daily(S_j,T,w_j):
    """Cost of fuel per day depending on the speed of the ship and weight carried """
    v = S_j/(24*T)
    k = 3.91e-06
    p = 381
    g = 3.1
    h = 2/3
    A = 49000
    DWT = 50000
    return k * (p + v**g) * ((w_j * DWT + A)**h)

def cost_loading(S,w,T):
    """Determines the cost of loading the ship"""
    return C_j_u + c_j_f * fuel_daily(S,T,w) * T

def T_conv(T_days):
    T_years = T_days / 365
    return T_years

def E_exp(r0, T_days):
    """Calculates the expected value of the exonent of a random variable, In the O.U. setting the expected value of the function (exp) of random variable with normal distribution given (O.U.) obrained using Ito Lemma, can be checked again log-normal distribution"""
    T = T_conv(T_days)
    s2_T = 0.5 * r_sigma * r_sigma * (1 - np.exp(-2 * r_lambda * T))
    mu_T = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    R_T = np.exp(mu_T + 0.5*s2_T)
    return R_T


def V_exp(r0, T_days):
    """Calculates the expected value of the exonent of a random variable, In the O.U. setting the expected value of the function (exp) of random variable with normal distribution given (O.U.) obrained using Ito Lemma, can be checked again log-normal distribution"""
    T = T_conv(T_days)
    s2_T = 0.5 * r_sigma * r_sigma * (1 - np.exp(-2 * r_lambda * T))
    mu_T = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    V_T = (np.exp(s2_T)-1) * np.exp(2*mu_T + s2_T)
    return V_T



def shipping_futures(E_R,T_days):
    """Discounts the point on the curve at time T, given that the expecation of the rate at T is R"""
    T = T_conv(T_days)
    return E_R * np.exp(fut_curve_slope * T)

def H(j,d,r_state,h_a,T):
    """Calculates reward/risk of h_a,T decision based on current state r_state and d at step j"""
    r0 = r_state + r_bar
    R = h_a * shipping_futures(E_exp(r0,T),T) + (1-h_a) * E_exp(r0,T)
    C_l_j = cost_loading(S[j],w[j],T)
    int_f_TCH = f_TCH * (1-np.exp(-alpha*T))/alpha
    V_R =  ((1-h_a)**2) * V_exp(r0,T)
    return (R - C_u_j) * np.exp(-alpha * T) - C_l_j - int_f_TCH - eta * V_R

def H_wrapper(x,params):
    """Wrapper function to H function for optimization"""
    h_a = x[0]
    T = x[1]
    j = params[0]
    d = params[1]
    r_state = params[2]
    return -H(j,d,r_state,h_a,T)


def H_wrapper_G0(x,params):
    """Wrapper function to H function for optimization of the last step with the FFP of G0"""
    h_a = x[0]
    T = x[1]
    j = params[0]
    d = params[1]
    r_state = params[2]
    return -(H(j,d,r_state,h_a,T) + np.exp(-alpha * T) * G0)

def H_star(j,d,r_state):
    """Funcation which returns the optimized H value and decisions best_h_a and best_T which produce this optimal outcome on the last step of the optimization"""
    T_min = S[j]/(24*v_max)
    T_max = S[j]/(24*v_min)
    T0 = 0.5 * (T_min + T_max)
    x0 = [0.5,T0]
    
    res_2 = minimize(H_wrapper_G0, 
                     x0, 
                     args=[j,d,r_state], 
                     method='trust-constr', 
                     bounds=[(0,1),(T_min,T_max)], 
                     options={'disp': False, 'maxiter': 500})
    mH = res_2.fun
    best_h_a = res_2.x[0]
    best_T = res_2.x[1]
    return mH, best_h_a, best_T
    
def p(r_s,r_state,T_days,states_step):
    """Calculates conditional transition probability of state r_s given state r_state at time T"""
    T = T_conv(T_days)
    r0 = r_state + r_bar
    r_sT = r_s + r_bar
    E = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    s1 = (1 - np.exp(-2 * r_lambda * T))
    S = s1 * np.sqrt(r_sigma * r_sigma / (2 * r_lambda))
    prob = norm.pdf(r_sT,E,S) * states_step        
    return prob

def G(j,d,r_state,h_a,T):
    """Intermediary value function for the optimization process, where G_star is already an optimal choice by on the future step"""
    E = 0
    states,states_step = gen_states(r_state,r_sigma*np.sqrt(T))
    for r_s in states:
        E = E + p(r_s,r_state,T,states_step) * G_star(j+1,d,r_s)
    return H(j,d,r_state,h_a,T) + np.exp(-alpha * T) * E 

def G_wrapper(x,params):
    """Wrapper function for G"""
    h_a = x[0]
    T = x[1]
    j = params[0]
    d = params[1]
    r_state = params[2]
    return -G(j,d,r_state,h_a,T)


def G_star(j,d,r_state):
    """The main optimization function"""
    if j==m:
        # if this is the last step, just
        return H_star(j,d,r_state)
    else:
        T_min = S[j]/(24*v_max)
        T_max = S[j]/(24*v_min)
        T0 = 0.5 * (T_min + T_max)
        x0 = [0.5,T0]
        
        res_2 = minimize(G_wrapper, x0, args=[j,d,r_state], method='trust-constr', bounds=[(0,1),(T_min,T_max)], options={'disp': False, 'maxiter': 500})
        mH = res_2.fun
        best_h_a = res_2.x[0]
        best_T = res_2.x[1]
        return mH, best_h_a, best_T
