# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 21:16:54 2021

@author: oleg_
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 21:19:36 2021

"""
import numpy as np
from scipy.stats import norm

"""
# O.U. parameters of the log rate model:
# the data is assuming T=1 it is one year.
r_bar = 6.2 #  e^6.2 ~ 500 (Corresponds to BIDY Index, Baltic Dirty Tanker Index)  
r_lambda = 3 # 0.5
r_sigma = 0.3 #0.4
eta = 0.000001 # risk aversion of the shipping agent
grid_size = np.max([r_lambda,60])
fut_curve_slope = -0.1 # negative means that there is a discount to the expectation built-in by the model
half_life = np.log(2)/r_lambda
"""

global r_bar #= 6.2 #  e^6.2 ~ 500 (Corresponds to BIDY Index, Baltic Dirty Tanker Index)  
global r_lambda #= 3 # 0.5
global r_sigma #= 0.3 #0.4

global eta #= 0.000001 # risk aversion of the shipping agent
global grid_size #= np.max([r_lambda,60])
global fut_curve_slope #= -0.1 # negative means that there is a discount to the expectation built-in by the model
global half_life #= np.log(2)/r_lambda


z_score = 5
grid_size = np.max([r_lambda,60])

def T_conv(T_days):
    T_years = T_days / 365
    return T_years


def gen_states(T_days):
    T = T_conv(T_days)
    
    states = np.linspace(start = r_bar-r_sigma*z_score*np.sqrt(T),
                         stop = r_bar+r_sigma*z_score*np.sqrt(T),
                         num = grid_size)
    states_step = np.sqrt(T)*r_sigma*z_score*2 / (grid_size-1)
    return states,states_step


def gen_cond_states(r,T_days,G_star_x_max, G_star_x_min):
    T = T_conv(T_days)
    UB = np.min([G_star_x_max,r+r_sigma*z_score*np.sqrt(T)])
    LB = np.max([G_star_x_min,r-r_sigma*z_score*np.sqrt(T)])
    states = np.linspace(start = LB,
                         stop = UB,
                         num = grid_size)
    states_step = np.sqrt(T)*r_sigma*z_score*2 / (grid_size-1)
    return states,states_step



def E_exp(r0, T_days):
    """Calculates the expected value of the exponent of a random variable, 
    In the O.U. setting the expected value of the function (exp) of random 
    variable with normal distribution given (O.U.) obrained using Ito Lemma, 
    can be checked again log-normal distribution"""
    T = T_conv(T_days)
    s2_T = 0.5 * r_sigma * r_sigma * (1 - np.exp(-2 * r_lambda * T))
    mu_T = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    R_T = np.exp(mu_T + 0.5*s2_T)
    return R_T


def V_exp(r0, T_days):
    """Calculates the expected value of the exonent of a random variable, In the O.U. setting the expected value of the function (exp) of random variable with normal distribution given (O.U.) obrained using Ito Lemma, can be checked again log-normal distribution"""
    T = T_conv(T_days)
    s2_T = 0.5 * r_sigma * r_sigma * (1 - np.exp(-2 * r_lambda * T)) / r_lambda
    mu_T = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    V_T = (np.exp(s2_T)-1) * np.exp(2*mu_T + s2_T)
    return V_T



def shipping_futures(E_R,T_days):
    """Discounts the point on the curve at time T, given that the expecation of the rate at T is R"""
    T = T_conv(T_days)
    return E_R * np.exp(fut_curve_slope * T)


    
def p_vector(states,step,r_state,T_days):
    """Calculates conditional transition probability of state r_s given state r_state at time T"""
    T = T_conv(T_days)
    r0 = r_state
    r_sT = states
    E = r0 * np.exp(-r_lambda * T) + r_bar * (1-np.exp(-r_lambda * T))
    s1 = (1 - np.exp(-2 * r_lambda * T))
    S = r_sigma * np.sqrt(s1 / (2 * r_lambda))
    prob = norm.pdf(r_sT,E,S) * step        
    return prob

def generate_path(r,T_days):
    dT = T_conv(1)
    path = [r]
    for t in range(1,T_days+1):
        
        E = r * np.exp(-r_lambda * dT) + r_bar * (1-np.exp(-r_lambda * dT))
        s1 = (1 - np.exp(-2 * r_lambda * dT))
        S = r_sigma * np.sqrt(s1 / (2 * r_lambda))
        r_new = np.random.normal(E,S)
        r = r_new
        path.append(r)
    return np.round(np.array(path),5)
        

"""
import matplotlib.pyplot as plt
states, step = gen_states()
print(p_vector(states, step,r_bar,15))
plt.plot(p_vector(states, step,r_bar,15))
np.sum(p_vector(states, step,r_bar,15))

states,step = gen_cond_states(r_bar,365,10e+10, -10e+10)
plt.plot(states,p_vector(states,step,5.2,365))
print(sum(p_vector(states,step,5.2,365)))
"""
#plt.plot(generate_path(6,1000))

