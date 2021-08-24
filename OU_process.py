# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 13:42:56 2021

@author: oleg_
"""
import numpy as np
from scipy.stats import norm
from conventions import T_conv

class OU_process:
    
    def __init__(self, r_bar, r_lambda, r_sigma, r_start=None, r_end=None):
        
        self.r_bar = r_bar
        self.r_lambda = r_lambda
        self.r_sigma = r_sigma
        
        self.r_start=r_start
        self.r_end=r_end
                
        self.z_score = 6
        self.grid_size = np.max([r_lambda,60])
    
    def get_params(self):
        return {'r_bar': self.r_bar, 'r_lambda': self.r_lambda, 'r_sigma':self.r_sigma}
    
        
    def gen_states(self,T_days):
        T = T_conv(T_days)
        if self.r_start==None:
            start = self.r_bar-self.r_sigma*self.z_score*np.sqrt(T)
        else:
            start = self.r_start
            
        if self.r_end==None:
            end = self.r_bar+self.r_sigma*self.z_score*np.sqrt(T)
        else:
            end = self.r_end
        print("Starting state:"+str(start))
        print("Ending state:"+str(end))
            
        states = np.linspace(start = start,
                             stop = end,
                             num = int(self.grid_size))
        states_step = np.sqrt(T)*self.r_sigma*self.z_score*2 / (self.grid_size-1)
        return states,states_step

    
    def gen_cond_states(self,r,T_days,G_star_x_max, G_star_x_min):
        T = T_conv(T_days)
        #if self.r_start==None:
        LB = np.max([G_star_x_min,r-self.r_sigma*self.z_score*np.sqrt(T)])
        #else:
        #    LB = self.r_start

        #if self.r_end==None:
        UB = np.min([G_star_x_max,r+self.r_sigma*self.z_score*np.sqrt(T)])
        #else:
        #    UB = self.r_end

        states = np.linspace(start = LB,
                             stop = UB,
                             num = int(self.grid_size))
        states_step = (UB - LB) / (int(self.grid_size)-1)
        return states,states_step



    def E_exp(self, r0, T_days):
        """Calculates the expected value of the exponent of a random variable, 
        In the O.U. setting the expected value of the function (exp) of random 
        variable with normal distribution given (O.U.) obrained using Ito Lemma, 
        can be checked again log-normal distribution"""
        T = T_conv(T_days)
        if self.r_lambda==0:
            s2_T = self.r_sigma * self.r_sigma * T
        else:
            s2_T = 0.5 * self.r_sigma * self.r_sigma * (1 - np.exp(-2 * self.r_lambda * T)) / self.r_lambda
                
        mu_T = r0 * np.exp(-self.r_lambda * T) + self.r_bar * (1-np.exp(-self.r_lambda * T))
        R_T = np.exp(mu_T - 0*0.5*s2_T)
        return R_T

    
    def V_exp(self, r0, T_days):
        """Calculates the expected value of the exonent of a random variable, In the O.U. setting the expected value of the function (exp) of random variable with normal distribution given (O.U.) obrained using Ito Lemma, can be checked again log-normal distribution"""
        T = T_conv(T_days)
        if self.r_lambda==0:
            s2_T = self.r_sigma * self.r_sigma * T
        else:
            s2_T = 0.5 * self.r_sigma * self.r_sigma * (1 - np.exp(-2 * self.r_lambda * T)) / self.r_lambda

        mu_T = r0 * np.exp(-self.r_lambda * T) + self.r_bar * (1-np.exp(-self.r_lambda * T))
        V_T = (np.exp(s2_T)-1) * np.exp(2*mu_T + s2_T)
        return V_T


    """
    def shipping_futures(self, E_R, T_days):
        #Discounts the point on the curve at time T, given that the expecation of the rate at T is R
        T = self.T_conv(T_days)
        return E_R * np.exp(self.fut_curve_slope * T)
    """
    
        
    def p_vector(self,states,step,r_state,T_days):
        """Calculates conditional transition probability of state r_s given state r_state at time T"""
        T = T_conv(T_days)
        r0 = r_state
        r_sT = states
        E = r0 * np.exp(-self.r_lambda * T) + self.r_bar * (1-np.exp(-self.r_lambda * T))
        if self.r_lambda==0:
            S = self.r_sigma * np.sqrt(T)
        else:
            s1 = (1 - np.exp(-2 * self.r_lambda * T))
            S = self.r_sigma * np.sqrt(s1 / (2 * self.r_lambda))

        prob = norm.pdf(r_sT,E - 0.5*S*S,S) * step           
        return prob
    
    def generate_path(self,r,T_days):
        dT = T_conv(1)
        path = [r]
        for t in range(1,T_days+1):
            
            E = r * np.exp(-self.r_lambda * dT) + self.r_bar * (1-np.exp(-self.r_lambda * dT))
            if self.r_lambda==0:
                S = self.r_sigma * np.sqrt(dT)
            else:
                s1 = (1 - np.exp(-2 * self.r_lambda * dT))
                S = self.r_sigma * np.sqrt(s1 / (2 * self.r_lambda))
            

            r_new = np.random.normal(E,S)
            r = r_new
            path.append(r)
        return np.round(np.array(path),5)
            