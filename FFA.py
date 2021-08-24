# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 14:06:39 2021

@author: oleg_
"""
from conventions import T_conv
import numpy as np

class FFA:
    def __init__(self,fut_curve_slope):
        self.fut_curve_slope = fut_curve_slope
    
    def shipping_futures(self, E_R, T_days):
        """Discounts the point on the curve at time T, given that the expecation of the rate at T is R"""
        T = T_conv(T_days)
        return E_R * np.exp(self.fut_curve_slope * T)
