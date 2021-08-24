# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 10:03:45 2021

@author: oleg_
"""

from Logistics_stoch import CRoundtrip
from scipy.interpolate import interp1d

class Optimization_Problem:
    
    def __init__(self, oJourney, oShip, oPortList, process, ffa, eta):
        self.oJourney = oJourney
        self.oShip = oShip
        self.oPortList = oPortList
        self.process = process
        self.ffa = ffa
        self.eta = eta
        self.oRTList = []
    
    def run_optimization(self):
                       
        
        self.oPortList[0].CalculateCargoImplicationsData(oJ=self.oJourney, oS=self.oShip, oPList=self.oPortList)
        
        
        TotalJourneyDistance_nm = 0
        for i in range(1, 1 + self.oJourney.NrOfRoundtrips):
            rt = CRoundtrip(self.oJourney, self.oShip, self.oPortList, TripNr=i)
            TotalJourneyDistance_nm = TotalJourneyDistance_nm + rt.TotalDistance_nm
            self.oRTList.append(rt)
        
        self.oJourney.TotalJourneyDistance_nm = TotalJourneyDistance_nm
        
        rt0 = CRoundtrip(self.oJourney, self.oShip, self.oPortList, TripNr=0)
        rt0.RoundtripGoodwill = self.oJourney.FutureProfitPotential_USDperDay * 365 / self.oJourney.OpprtCostCapitalRate
        self.oRTList.insert(0,rt0)
        
        self.oRTList[1].FutureGoodwill = interp1d([-1e+10,1e+10],[rt0.RoundtripGoodwill,rt0.RoundtripGoodwill])
        
        #For each roundtrip, for each leg ... find the best goodwill
        for i in range(1,1+self.oJourney.NrOfRoundtrips):
            self.oRTList[i].FindBestGoodwill(self.oRTList, self.process, self.ffa, self.eta)
            if i < self.oJourney.NrOfRoundtrips:
                self.oRTList[i + 1].FutureGoodwill = self.oRTList[i].RoundtripGoodwill
        
        
        #Calculate solution metrics
        for i in range(self.oJourney.NrOfRoundtrips,0,-1):
            self.oRTList[i].CalculateRoundtripTime(self.oJourney)
        print(self.oRTList[0].RoundtripTime_Days)
        
        self.oJourney.CalculateJourneyTime(self.oRTList)
        self.oJourney.CalculateProfits(self.oRTList)
        return self.oRTList
    