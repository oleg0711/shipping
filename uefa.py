# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 11:14:51 2021

@author: oleg_
"""

def check_group(s):
    t1 = 0 in s
    t2 = 1 in s
    t3 = 2 in s
    tt = t1+t2+t3 
    if tt==0:
        return 0
    elif tt==3:
        return 1
    else:
        return 2
    
import numpy as np
t = 0
p = 0
for n in range(0,5000000):
    a = np.random.permutation(24)
    g = a.reshape(6,4)
    t = t + 1
    for i in range(0,6):
        c = check_group(g[i,:])
        if c==1:
            p = p + 1
            break
        elif c==2:
            break
        else:
            pass
        
print(p/t)