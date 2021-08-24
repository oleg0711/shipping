# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 18:18:55 2021

@author: oleg_
"""

# -*- coding: utf-8 -*-
"""

Created on Sat May 29 08:13:59 2021
@author: Patrick (VBA Code), Oleg (Python, reprogramming the excel spreadsheet in Python)
@author: Oleg (adaptation of the deterministic model to the stochastic setting)

"""

import pandas as pd
from Logistics import CShip
from Logistics import CRoundtrip
from Logistics import CJourney
from Logistics import CPort
import deterministic_opt_func as det_opt

print("NEW RUN ******")


oPTList = det_opt.run_deterministic_opt()