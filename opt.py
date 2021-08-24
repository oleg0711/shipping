# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 22:12:40 2021

@author: oleg_
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

n = 30
p1 = np.zeros(n)
p2 = np.zeros(n)

X = np.linspace(90,110,num=n)
Y = np.linspace(81,121,num=n)

for i in range(0,n):
    p1[i] = norm.pdf(X[i],100,3)
    p2[i] = norm.pdf(Y[i],100,3*np.sqrt(2))

p1 = p1 / sum(p1)
p2 = p2 / sum(p2)

P = np.zeros(shape=(n,n))
P0 = P.reshape(1,n**2)
eps = 1


def func(x,params):
    
    p1 = params[0]
    p2 = params[1]
    eps = params[2]
    P = x.reshape(n,n)
    v = np.dot(P,p1) - p2 
    r = 10000*np.dot(v,v.T) +  eps* sum(sum(np.log(P+0.00000001) * P))
    return r  

def sum_one(P):
    return np.dot(P.reshape(n,n).T, np.ones(n)) - np.ones(n) 

def pos_prob(P):
    return P - 0.00000000001

cons = ({'type': 'eq', 'fun': sum_one},
        {'type': 'ineq', 'fun': pos_prob})

params = {}
params[0] = p1
params[1] = p2
params[2] = eps


res = minimize(func, P0, params, method='SLSQP',constraints=cons,
               options={'maxiter': 10000, 
                        'ftol': 1e-07, 
                        'iprint': 1, 
                        'disp': True, 
                        'eps': 1e-7})
    
P = res.x.reshape(n,n)
plt.plot(P[:,-15])




