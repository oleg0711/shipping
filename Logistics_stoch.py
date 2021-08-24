# -*- coding: utf-8 -*-
"""
Created on Sat May 29 09:23:13 2021

@author: oleg_
"""
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import OU_process
import FFA

#CLASS MODULE: CShip
class CShip:

    Vmin:float
    Vmax:float
    
    DWTscantling:float
    DWTdesign:float #DWT carried by ship
    Lightweight:float
    k:float
    p:float
    g:float
    a:float #Fuel consumption function parameters
    ShipDischargeRate:float #m^3 per hour
    BallastCapacity:float #m^3 on ship
    MinFillRateShip:float
    AuxFuelConsumption_TonnePerDay:float


    def __init__(self, Vmin, Vmax, DWTscantling, DWTdesign, Lightweight, k, p, g, a, ShipDischargeRate, BallastCapacity,MinFillRateShip,AuxFuelConsumption_TonnePerDay):

        self.Vmin = Vmin #Sheet1.Cells(10, 3)
        self.Vmax = Vmax #Sheet1.Cells(10, 5)
        self.DWTscantling = DWTscantling #Sheet1.Cells(11, 3)
        self.DWTdesign = DWTdesign #Sheet1.Cells(11, 5)
        self.Lightweight = Lightweight #Sheet1.Cells(12, 3)
        self.k = k #Sheet1.Cells(13, 3)
        self.p = p #Sheet1.Cells(13, 5)
        self.g = g #Sheet1.Cells(13, 7)
        self.a = a #Sheet1.Cells(13, 9)
        self.ShipDischargeRate = ShipDischargeRate #Sheet1.Cells(14, 3)
        self.BallastCapacity = BallastCapacity #Sheet1.Cells(15, 3)
        self.MinFillRateShip = MinFillRateShip #Sheet1.Cells(16, 3)
        self.AuxFuelConsumption_TonnePerDay = AuxFuelConsumption_TonnePerDay #Sheet1.Cells(17, 3)
    
        #print("Ship data-------------------")
        #print("Vmin:" + str(self.Vmin))
        #print("Vmax:" + str(self.Vmax))
        #print("DWTscantling:" + str(self.DWTscantling))
        #print("DWTdesign:" + str(self.DWTdesign))
        #print("Lightweight:" + str(self.Lightweight))
        #print("k:" + str(self.k))
        #print("p:" + str(self.p))
        #print("g:" + str(self.g))
        #print("a:" + str(self.a))
        #print("ShipDischargeRate:" + str(self.ShipDischargeRate))
        #print("BallastCapacity:" + str(self.BallastCapacity))
        #print("MinFillRateShip:" + str(self.MinFillRateShip))
    
    def FuelConsumption_TonnePerDay(self, Speed:float, DWT:float):
        FuelConsumption_TonnePerDay = self.k * (self.p + (Speed) ** self.g) * (DWT + self.Lightweight) ** self.a
        return FuelConsumption_TonnePerDay



class CJourney:
    
    NrOfRoundtrips:int
    LegsPerRoundtrip:int
    OpprtCostCapitalRate:float
    DailyHire_USDperDay:float
    FutureProfitPotential_USDperDay:float
    TotalJourneyDistance_nm:float
    
    TotalJourneyTime_Days:float
    NPV_Profits_IncludingFPP_USD:float
    TCE_Journey_USDperDay:float
    JourneyProfits_USDperDay:float
    oRTList_det:list
    
    def __init__(self,
                 NrOfRoundtrips,
                 LegsPerRoundtrip,
                 OpprtCostCapitalRate,
                 DailyHire_USDperDay,
                 FutureProfitPotential_USDperDay):
        
        self.NrOfRoundtrips = NrOfRoundtrips #Sheet1.Cells(3, 3)
        self.LegsPerRoundtrip = LegsPerRoundtrip #Sheet1.Cells(4, 3)
        self.OpprtCostCapitalRate = OpprtCostCapitalRate #Sheet1.Cells(5, 3)
        self.DailyHire_USDperDay = DailyHire_USDperDay #Sheet1.Cells(6, 3)
        self.FutureProfitPotential_USDperDay = FutureProfitPotential_USDperDay #Sheet1.Cells(7, 3)
        
        self.TotalJourneyTime_Days = 0
        
        #print ("NrOfRoundtrips:" + str(NrOfRoundtrips))
        #print ("LegsPerRoundtrip:" + str(LegsPerRoundtrip))
        #print ("OpprtCostCapitalRate:" + str(OpprtCostCapitalRate))
        #print ("DailyHire_USDperDay:" + str(DailyHire_USDperDay))
        #print ("FutureProfitPotential_USDperDay:" + str(FutureProfitPotential_USDperDay))
        
    
    def CalculateJourneyTime(self, oRTList:list):
        #print('Calculating total journey time...')
        TotalJourneyTime_Days = 0
        for i in range(1,1+self.NrOfRoundtrips):
            TotalJourneyTime_Days = TotalJourneyTime_Days + oRTList[i].RoundtripTime_Days
        self.TotalJourneyTime_Days = TotalJourneyTime_Days
        return TotalJourneyTime_Days

    
    def CalculateProfits(self, oRTList:list):

        #NPV profits:
        # profits needs to be calculated per state:
        #states = oRTList[self.NrOfRoundtrips].RoundtripGoodwill.x
        
        self.NPV_Profits_IncludingFPP_USD = oRTList[self.NrOfRoundtrips].RoundtripGoodwill
        #TCE EXCLUDING Daily Hire AND EXCLUDING Future Profit Potential
        #This calculates the AS value (as a USD/year rate) received over the TotalJourneyTime_Days:
        
        #self.TCE_Journey_USDperDay = interp1d(states,(self.NPV_Profits_IncludingFPP_USD.y - (self.FutureProfitPotential_USDperDay * 365 / self.OpprtCostCapitalRate) * np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365)) * self.OpprtCostCapitalRate / (1 - np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365)))
        #self.TCE_Journey_USDperDay = interp1d(states,self.TCE_Journey_USDperDay.y / 365)
        
        #This calculates the profit over the journey as a per day measure:
        #self.JourneyProfits_USDperDay = self.TCE_Journey_USDperDay
        #This adjusts by adding back the Daily Hire (TCH), as not incurred in TCE calculation
        #self.TCE_Journey_USDperDay = interp1d(states,self.TCE_Journey_USDperDay.y + self.DailyHire_USDperDay)
        
   
class CLeg:

    PortNr_StartOfLeg:int
    PortNr_EndOfLeg:int

    RoundtripNr:int
    LegNr:int

    Distance_nm:float
    Speed_kn:float

    TimeLoading_Days:float
    TimeWaiting_Days:float
    TimeUnloading_Days:float

    TimeElapsedOnRoundtripSoFar_Days:float

    TimeAtSea_Days:float
    TimeAtSeaMax_Days:float
    TimeAtSeaMin_Days:float
    TimeStepSize_Days:float
    LegTime_Days:float


    CargoCarried_Tonne:float
    MainBunkerConsumption_TonnePerDay:float
    MainBunkerCarried_Tonne:float
    DWTCarried_Tonne:float
    AuxFuelCarried_Tonne:float

    LoadingCost_OtherOperations_USD:float
    LoadingCost_MainBunker_USD:float
    LoadingCost_AuxFuel_USD:float

    UnloadingCosts_USD:float
    Revenues_USD:float

    Profitability_ThisLeg_USD:float
    FutureGoodwill:float
    Goodwill_USD:float
    Best_Hedge_Ratio:float

    def __init__(self,
                 LegNr,
                 RoundtripNr,
                 Distance_nm,
                 PortNr_StartOfLeg,
                 PortNr_EndOfLeg,
                 TimeLoading_Days,
                 TimeWaiting_Days,
                 TimeUnloading_Days,
                 CargoCarried_Tonne,
                 LoadingCost_OtherOperations_USD,
                 UnloadingCosts_USD,
                 Revenues_Barrels):
    
        self.LegNr = LegNr
        self.RoundtripNr = RoundtripNr
        self.Distance_nm = Distance_nm
        self.PortNr_StartOfLeg = PortNr_StartOfLeg
        self.PortNr_EndOfLeg = PortNr_EndOfLeg
        self.TimeLoading_Days = TimeLoading_Days
        self.TimeWaiting_Days = TimeWaiting_Days
        self.TimeUnloading_Days = TimeUnloading_Days
        self.CargoCarried_Tonne = CargoCarried_Tonne
        self.LoadingCost_OtherOperations_USD = LoadingCost_OtherOperations_USD
        self.UnloadingCosts_USD = UnloadingCosts_USD
        self.Revenues_Barrels = Revenues_Barrels
        self.Best_Hedge_Ratio = 0
        self.Profitability = 0

    
    def CalcGoodwillWrapper(self, x, params):
        
        Goodwill_USD, Profitability = self.CalcGoodwill(x,params)
        return -Goodwill_USD
    
    def CalcGoodwill(self, x, params):
        
        hedge_ratio = x[0]
        TimeAtSea_Days = x[1]
        
        r = params[0]
        oS = params[1]
        oPList = params[2]
        oJ = params[3]
        ou_process = params[4]
        ffa = params[5]
        eta = params[6]
        
        LegTime_Days = self.TimeLoading_Days + TimeAtSea_Days + self.TimeWaiting_Days + self.TimeUnloading_Days
       
        AuxFuelCarried_Tonne = (self.TimeLoading_Days + self.TimeWaiting_Days + self.TimeUnloading_Days) * oS.AuxFuelConsumption_TonnePerDay
        LoadingCost_AuxFuel_USD = AuxFuelCarried_Tonne * oPList[self.PortNr_StartOfLeg].AuxFuelRate_USDperTonne
       
       
        self.Speed_kn = self.Distance_nm / (24 * TimeAtSea_Days)
       
        #Minimum of Cargo carried or ballast water carried:
        self.DWTCarried_Tonne = self.CargoCarried_Tonne
        if self.DWTCarried_Tonne < self.MinimumLoad:
            self.DWTCarried_Tonne = self.MinimumLoad
       
        #Get Fuelconsumption and fuel needed using estimate for DWT carried
        self.MainBunkerConsumption_TonnePerDay = oS.FuelConsumption_TonnePerDay(self.Speed_kn, self.DWTCarried_Tonne)
        self.MainBunkerCarried_Tonne = self.MainBunkerConsumption_TonnePerDay * TimeAtSea_Days
        self.LoadingCost_MainBunker_USD = self.MainBunkerCarried_Tonne * oPList[self.PortNr_StartOfLeg].MainBunkerRate_USDperTonne
        
        #Calculate goodwill and retain highest value solution
        
        df = np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))
        
        Cost_ThisLeg_USD = ( - self.UnloadingCosts_USD) * df - self.LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - self.LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 * (1 - df) / (oJ.OpprtCostCapitalRate)
        
        E_R = ou_process.E_exp(r, LegTime_Days)
        Future_Price = ffa.shipping_futures(E_R, LegTime_Days)
        expected_price = (1-hedge_ratio) * E_R + hedge_ratio * Future_Price
        variance_price = ((1-hedge_ratio)**2) * ou_process.V_exp(r,LegTime_Days)
        """
        print('State:' + str(r))
        print('E_R:' + str(E_R))
        print('Future_Price:' + str(Future_Price))
        print('expected_price:' + str(expected_price))
        print('vol_price:' + str(np.sqrt(variance_price)))
        print('Leg_time:' + str(LegTime_Days))
        
        if input()=='q':
            sys.exit()
        """
        
        
        Revenue_ThisLeg_USD = (expected_price/1000) * df * self.Revenues_Barrels
        Penalty_USD = -df * self.Revenues_Barrels * eta * (variance_price)
        
        
        Profitability_ThisLeg_USD = Revenue_ThisLeg_USD + Cost_ThisLeg_USD + Penalty_USD  
        G_star_x_max = 1e+100 #np.max(self.FutureGoodwill.x)
        G_star_x_min = -1e+100 #np.min(self.FutureGoodwill.x)
        
        states_next , step_next = ou_process.gen_cond_states(r,self.TimeAtSeaMax_Days,G_star_x_max, G_star_x_min)
        prob_space = ou_process.p_vector(states_next, step_next, r, LegTime_Days)
        
        #prob_space = prob_space * 0 + 1
        #prob_space = prob_space/sum(prob_space)
        #print("sum:"+str(np.sum(prob_space)))
        #print("max:"+str(np.max(prob_space)))
        #print("min:"+str(np.min(prob_space)))
        #prob_space = prob_space/np.sum(prob_space)
        
        Goodwill_USD = Profitability_ThisLeg_USD + df * np.dot(prob_space, 
                                                               self.FutureGoodwill(states_next))
        
        print("Prob space:"+str(np.sum(prob_space)))
        return Goodwill_USD,Profitability_ThisLeg_USD

    def FindBestGoodwill(self,
                         ou_process:OU_process,
                         ffa: FFA,
                         eta,
                         oJ:CJourney, 
                         oS:CShip, 
                         oRTList = [], 
                         oPList = []): #list of round trips, list of ports
   
   
        #If last leg, get future goodwill from next roundtrip
        if self.LegNr == oJ.LegsPerRoundtrip:
            self.FutureGoodwill = oRTList[self.RoundtripNr].FutureGoodwill
            #Else get it from the next leg of this roundtrip
        else:
            self.FutureGoodwill = oRTList[self.RoundtripNr].FindGoodwillOfNextLeg(self.LegNr)
            
        print("RoundtripNr:" + str(self.RoundtripNr))
        print("LegNr:" + str(self.LegNr))
        plt.plot(self.FutureGoodwill.x,self.FutureGoodwill.y)
        plt.show()
        print("FutureGoodwill:" + str(self.FutureGoodwill))
         
        #Figure out the range of Tmax and Tmin and set a stepsize
        self.TimeAtSeaMax_Days = self.Distance_nm / (24 * oS.Vmin)
        self.TimeAtSeaMin_Days = self.Distance_nm / (24 * oS.Vmax)
        #***************************************
        #CHOOSE A STEPSIZE SAY EVERY HOUR
        self.TimeStepSize_Days = 1 / (24 * 6)
        #***************************************
        
        Tapprx_days = 0
        for i in range(1,self.RoundtripNr+1):
            Tapprx_days = Tapprx_days + oJ.oRTList_det[i].RoundtripTime_Days
        
                
        states,step = ou_process.gen_states(Tapprx_days) 
        
        G_star = []
        BestTimeAtSea_Days_array = []    
        BestSpeed_kn_array = []
        Best_Hedge_Ratio_array = []
        Best_Profitability_array = []
            
        for state in states:
            print("Running optimization for state: " + str(state))
            r = state # log shipping spot rate

            self.BestGoodwill_USD:float = -999999999
            self.BestTimeAtSea_Days:float = 9999999
            self.BestSpeed_kn:float = 9999999
            self.MinimumLoad:float = oS.DWTdesign * oS.MinFillRateShip

            
            params = [r,oS, oPList, oJ, ou_process, ffa, eta]
            
            staring_point_for_hedge_ratio_search = 0.5
            staring_point_for_time_at_sea_search = 0.5 * (self.TimeAtSeaMin_Days + self.TimeAtSeaMax_Days)
            
            x0 = [staring_point_for_hedge_ratio_search,
                  staring_point_for_time_at_sea_search]
            max_hedge_ratio = 1
            min_hedge_ratio = 0
            bounds = ((min_hedge_ratio, max_hedge_ratio),
                      (self.TimeAtSeaMin_Days, self.TimeAtSeaMax_Days))
            
            res = minimize(self.CalcGoodwillWrapper,
                           x0, 
                           args=params, 
                           bounds=bounds,
                           method='L-BFGS-B',
                           options={'disp':True,'maxiter':1000})
            
            self.Best_Hedge_Ratio = res.x[0]
            self.BestGoodwill_USD = -res.fun
            self.BestTimeAtSea_Days = res.x[1]
            self.BestSpeed_kn = self.Distance_nm / (24 * res.x[1])
            g,p = self.CalcGoodwill(res.x,params)
            self.Profitability = p
            
            """
            for hedge_ratio in np.linspace(0,1,11):
                for TimeAtSea_Days in np.linspace(start=self.TimeAtSeaMin_Days,stop=self.TimeAtSeaMax_Days,num=num):
                   
                    Goodwill_USD = self.CalcGoodwill([hedge_ratio,TimeAtSea_Days], params)
                    
                    if Goodwill_USD > self.BestGoodwill_USD:
                         self.Best_Hedge_Ratio = hedge_ratio
                         self.BestGoodwill_USD = Goodwill_USD
                         self.BestTimeAtSea_Days = TimeAtSea_Days
                         self.BestSpeed_kn = self.Speed_kn
                         self.Profitability = 0 
            """
            
            #print("Best speed: " + str(self.BestSpeed_kn))
            #print("Best hedge ratio: " + str(self.Best_Hedge_Ratio))
            #print("Best time at sea: " + str(self.BestTimeAtSea_Days))
            #print("-------------------------------------------------")
            G_star.append(self.BestGoodwill_USD)
            BestTimeAtSea_Days_array.append(self.BestTimeAtSea_Days)
            BestSpeed_kn_array.append(self.BestSpeed_kn)
            Best_Hedge_Ratio_array.append(self.Best_Hedge_Ratio)
            Best_Profitability_array.append(self.Profitability)
            
        
        self.Goodwill_USD = interp1d(states,np.array(G_star),fill_value='extrapolate',kind='cubic')
        self.Profitability_ThisLeg_USD = interp1d(states,Best_Profitability_array,fill_value='extrapolate',kind='cubic')
        
        self.TimeAtSea_Days = interp1d(states,np.array(BestTimeAtSea_Days_array),kind='cubic')
        self.Speed_kn = interp1d(states,BestSpeed_kn_array,kind='cubic')
        self.Hedge_Ratio = interp1d(states,Best_Hedge_Ratio_array,kind='cubic')
        
        #print("LoadingCost_OtherOperations_USD:" + str(self.LoadingCost_OtherOperations_USD))
        #print("LoadingCost_MainBunker_USD:" + str(self.LoadingCost_MainBunker_USD))
        #print("UnloadingCosts_USD:" + str(self.UnloadingCosts_USD))

   
       #If this leg is the first of the roundtrip, store the Goodwill also into this roundtrips goodwill
        if self.LegNr == 1:
            oRTList[self.RoundtripNr].RoundtripGoodwill = self.Goodwill_USD
            #print("oRT(RoundtripNr).RoundtripGoodwill=:" + str(oRTList[self.RoundtripNr].RoundtripGoodwill))
       
   
class CPort:

    DistancePreviousPort_nm:float
    LoadingRate_QbmetresperHr:float
    WaitingTime_Hrs:float

    
    CargoIntake_Barrels:float 
    CargoIntake_Tonne:float
    CargoIntake_QBmetresperBarrel:float
    CargoIntake_QbmetresPerTonne:float
    
    UnloadingTime_Hrs:float
    LoadingTime_Hrs:float
    
    FixedPortAccessCosts_USD:float
    UnloadingCharge_USDperHr:float
    LoadingCharge_USDperHr:float
    
    UnloadingCosts_USD:float
    LoadingCosts_USD:float
    
    CargoRevenueRate_USDperBarrelper1000nm:float    
    LegRevenue_USD:float #Associated with unloading
    MainBunkerRate_USDperBarrel:float
    MainBunkerRate_USDperTonne:float
    MainBunker_QBmetresperBarrel:float
    MainBunker_QbmetresPerTonne:float
    
    AuxFuelRate_USDperTonne:float
    
    def __init__(self, 
                 PortNr, 
                 DistancePreviousPort_nm,
                 LoadingRate_QbmetresperHr,
                 WaitingTime_Hrs,
                 CargoIntake_Barrels,
                 CargoIntake_QBmetresperBarrel,
                 CargoIntake_QbmetresPerTonne,
                 CargoRevenueRate_USDperBarrelper1000nm,
                 FixedPortAccessCosts_USD,
                 UnloadingCharge_USDperHr,
                 LoadingCharge_USDperHr,
                 MainBunkerRate_USDperBarrel,
                 MainBunker_QBmetresperBarrel,
                 MainBunker_QbmetresPerTonne,
                 AuxFuelRate_USDperTonne,
                 
                 UnloadingCosts_USD,
                 LoadingCosts_USD,
                 CargoIntake_Tonne,
                 UnloadingTime_Hrs,
                 LoadingTime_Hrs,
                 
                 LegRevenue_Barrels,
                 LegRevenue_USD):

        if PortNr > 0: 
            self.DistancePreviousPort_nm = DistancePreviousPort_nm #Sheet1.Cells(21, 4 + PortNr) 
        else: 
            self.DistancePreviousPort_nm = 0
        
        self.LoadingRate_QbmetresperHr = LoadingRate_QbmetresperHr #Sheet1.Cells(22, 4 + PortNr)
        self.WaitingTime_Hrs = WaitingTime_Hrs #Sheet1.Cells(23, 4 + PortNr)
         
        self.CargoIntake_Barrels = CargoIntake_Barrels #Sheet1.Cells(26, 4 + PortNr)
        self.CargoIntake_QBmetresperBarrel = CargoIntake_QBmetresperBarrel #Sheet1.Cells(27, 4 + PortNr)
        self.CargoIntake_QbmetresPerTonne = CargoIntake_QbmetresPerTonne #Sheet1.Cells(28, 4 + PortNr)
         
        self.FixedPortAccessCosts_USD = FixedPortAccessCosts_USD #Sheet1.Cells(41, 4 + PortNr)
        if PortNr > 0:
            self.UnloadingCharge_USDperHr = UnloadingCharge_USDperHr #Sheet1.Cells(42, 4 + PortNr) 
        else: 
            self.UnloadingCharge_USDperHr = 0
        
        self.LoadingCharge_USDperHr = LoadingCharge_USDperHr #Sheet1.Cells(43, 4 + PortNr)   
        self.CargoRevenueRate_USDperBarrelper1000nm = CargoRevenueRate_USDperBarrelper1000nm #Sheet1.Cells(36, 4 + PortNr)
        self.MainBunkerRate_USDperBarrel = MainBunkerRate_USDperBarrel #Sheet1.Cells(49, 4 + PortNr)
        self.MainBunker_QBmetresperBarrel = MainBunker_QBmetresperBarrel #Sheet1.Cells(51, 4 + PortNr)
        self.MainBunker_QbmetresPerTonne = MainBunker_QbmetresPerTonne #Sheet1.Cells(52, 4 + PortNr)
        self.AuxFuelRate_USDperTonne = AuxFuelRate_USDperTonne #Sheet1.Cells(53, 4 + PortNr)
         
        self.CargoIntake_Tonne = CargoIntake_Tonne
        self.UnloadingTime_Hrs = UnloadingTime_Hrs
        self.LoadingTime_Hrs = LoadingTime_Hrs
        self.LoadingCosts_USD = LoadingCosts_USD
        self.UnloadingCosts_USD = UnloadingCosts_USD
        self.LegRevenue_Barrels = LegRevenue_Barrels    
        self.LegRevenue_USD = LegRevenue_USD    
        
        #print("Port data:" + str(PortNr))
        #print("DistancePreviousPort_nm:" + str(self.DistancePreviousPort_nm))
        
        #print("LoadingRate_QbmetresperHr:" + str(self.LoadingRate_QbmetresperHr))
        #print("WaitingTime_Hrs:" + str(self.LoadingRate_QbmetresperHr))
        
        #print("CargoIntake_Barrels:" + str(self.CargoIntake_Barrels))
        #print("CargoIntake_QBmetresperBarrel:" + str(self.CargoIntake_QBmetresperBarrel))
        #print("CargoIntake_QbmetresPerTonne:" + str(self.CargoIntake_QbmetresPerTonne))
          
        #print("FixedPortAccessCosts_USD:" + str(self.FixedPortAccessCosts_USD))
        #print("UnloadingCharge_USDperHr:" + str(self.UnloadingCharge_USDperHr))
        #print("LoadingCharge_USDperHr:" + str(self.LoadingCharge_USDperHr))
        
        #print("CargoRevenueRate_USDperBarrelper1000nm:" + str(self.CargoRevenueRate_USDperBarrelper1000nm))
        #print("MainBunkerRate_USDperBarrel:" + str(self.MainBunkerRate_USDperBarrel))
        #print("MainBunker_QBmetresperBarrel:" + str(self.MainBunker_QBmetresperBarrel))
        #print("MainBunker_QbmetresPerTonne:" + str(self.MainBunker_QbmetresPerTonne))
        
        #print("AuxFuelRate_USDperTonne:" + str(self.AuxFuelRate_USDperTonne))
    
    
    
    
    def CalculateCargoImplicationsData(self, oJ, oS, oPList): #list of ports
        print ("********* Running CalculateCargoImplicationsData")
        for i in range(0,1+oJ.LegsPerRoundtrip):
            #print("i="+str(i))
            #print("oPList[i].CargoIntake_Barrels="+str(oPList[i].CargoIntake_Barrels))
            
            if oPList[i].CargoIntake_Barrels > 0:
                oPList[i].CargoIntake_Tonne = oPList[i].CargoIntake_Barrels * oPList[i].CargoIntake_QBmetresperBarrel / oPList[i].CargoIntake_QbmetresPerTonne
            else:
                oPList[i].CargoIntake_Tonne = 0
            
            
            if i > 0:
                oPList[i].UnloadingTime_Hrs = oPList[i - 1].CargoIntake_Barrels * oPList[i - 1].CargoIntake_QBmetresperBarrel / oS.ShipDischargeRate
            else:
               oPList[i].UnloadingTime_Hrs = 0
            
            
            oPList[i].LoadingTime_Hrs = oPList[i].CargoIntake_Barrels * oPList[i].CargoIntake_QBmetresperBarrel / oPList[i].LoadingRate_QbmetresperHr
         
            oPList[i].UnloadingCosts_USD = oPList[i].UnloadingTime_Hrs * oPList[i].UnloadingCharge_USDperHr
         
         
            oPList[i].LoadingCosts_USD = oPList[i].LoadingTime_Hrs * oPList[i].LoadingCharge_USDperHr
         
            if i > 0:
                # oPList[i - 1].CargoRevenueRate_USDperBarrelper1000nm
                oPList[i].LegRevenue_Barrels = (oPList[i].DistancePreviousPort_nm / 1000) * oPList[i - 1].CargoIntake_Barrels
            else:
               oPList[i].LegRevenue_Barrels = 0
            
         
            if i < oJ.LegsPerRoundtrip:
               oPList[i].MainBunkerRate_USDperTonne = (oPList[i].MainBunkerRate_USDperBarrel / oPList[i].MainBunker_QBmetresperBarrel) * oPList[i].MainBunker_QbmetresPerTonne
            else:
               oPList[i].MainBunkerRate_USDperTonne = 0
            
         
            #print("Port: " + str(i))
            #print("CargoIntake_Tonne:" +str(oPList[i].CargoIntake_Tonne))
            #print("UnloadingTime_Hrs:" +str(oPList[i].UnloadingTime_Hrs))
            #print("LoadingTime_Hrs:" +str(oPList[i].LoadingTime_Hrs))
            #print("UnloadingCosts_USD:" +str(oPList[i].UnloadingCosts_USD))
            #print("LoadingCosts_USD:" +str(oPList[i].LoadingCosts_USD))
            #print("LegRevenue_Barrels:" +str(oPList[i].LegRevenue_Barrels))
            #print("MainBunkerRate_USDperTonne:" +str(oPList[i].MainBunkerRate_USDperTonne))
            #print("---------------------------------------------")
        
        return oPList





class CRoundtrip():

    TripNr:int
    RoundtripGoodwill:float
    oLegList:list
    FutureGoodwill:float
    RoundtripTime_Days:float
    oJ:CJourney
    oS:CShip
    oPList:list
    TotalDistance_nm:float    
    def __init__(self, oJ:CJourney, oS:CShip, oPList, TripNr):
        
        self.oJ = oJ
        self.oS = oS
        self.oPList = oPList
        self.oLegList = []
        self.RoundtripTime_Days = 0
        self.TripNr = TripNr
        self.TotalDistance_nm = 0
        
        SPLIT = 1
        
        for i in range(1,1+oJ.LegsPerRoundtrip):
            #print("Define leg" + str(i))            
            oLeg = CLeg(LegNr=i,
                        RoundtripNr=TripNr,
                        Distance_nm = oPList[i].DistancePreviousPort_nm,
                        PortNr_StartOfLeg = i - 1,
                        PortNr_EndOfLeg = i,
                        TimeLoading_Days = oPList[i - 1].LoadingTime_Hrs / 24,
                        TimeWaiting_Days = oPList[i].WaitingTime_Hrs / 24,
                        TimeUnloading_Days = oPList[i].UnloadingTime_Hrs / 24,
                        CargoCarried_Tonne = oPList[i-1].CargoIntake_Tonne,
                        LoadingCost_OtherOperations_USD = oPList[i - 1].LoadingCosts_USD + SPLIT * (oPList[i - 1].FixedPortAccessCosts_USD),
                        UnloadingCosts_USD = oPList[i].UnloadingCosts_USD + (1 - SPLIT) * (oPList[i].FixedPortAccessCosts_USD),
                        Revenues_Barrels = oPList[i].LegRevenue_Barrels)
                        
            self.oLegList.append(oLeg)
            self.TotalDistance_nm = self.TotalDistance_nm + oPList[i].DistancePreviousPort_nm      
            #print("Leg data: " + str(i))
            #print("LegNr:" + str(oLeg.LegNr))
            #print("RoundtripNr:" + str(oLeg.RoundtripNr))
            #print("Distance_nm:" + str(oLeg.Distance_nm))
            #print("PortNr_StartOfLeg:" + str(oLeg.PortNr_StartOfLeg))
            #print("PortNr_EndOfLeg:" + str(oLeg.PortNr_EndOfLeg))
            #print("TimeLoading_Days:" + str(oLeg.TimeLoading_Days))
            #print("TimeWaiting_Days:" + str(oLeg.TimeWaiting_Days))
            #print("TimeUnloading_Days:" + str(oLeg.TimeUnloading_Days))
            #print("CargoCarried_Tonne:" + str(oLeg.CargoCarried_Tonne))
            #print("LoadingCost_OtherOperations_USD:" + str(oLeg.LoadingCost_OtherOperations_USD))
            #print("UnloadingCosts_USD:" + str(oLeg.UnloadingCosts_USD))
            #print("Revenues_Barrels:" + str(oLeg.Revenues_Barrels))


    def FindBestGoodwill(self, oRTList, ou_process, FFA, eta):
        for j in range(self.oJ.LegsPerRoundtrip-1, -1, -1):
            self.oLegList[j].FindBestGoodwill(ou_process, FFA, eta, self.oJ, self.oS, oRTList, self.oPList)
    
    
    def FindGoodwillOfNextLeg(self,CurrentLegNr:int):
        #print("Size of leg list:" + str(len(self.oLegList)))
        #print("CurrentLegNr" + str(CurrentLegNr+1))
        return self.oLegList[CurrentLegNr].Goodwill_USD

    def PrintOptimalLegSpeeds(self, oJ):
        df = pd.DataFrame(columns=['Leg Nr','Days at Sea','Speed(knots)','Roundtrip Time Days'])
        for j in range(0,oJ.LegsPerRoundtrip):
            df = df.append({'Leg Nr':j,
                            'Days at Sea':self.oLegList[j].TimeAtSea_Days,
                            'Speed(knots)':self.oLegList[j].Speed_kn,
                            'Roundtrip Time Days':self.RoundtripTime_Days},
                ignore_index=True)
            #print("Leg Nr: " + str(j))
            #print("Days at Sea (Days):" + str(self.oLegList[j].TimeAtSea_Days))
            #print("Speed(knots):" + str(self.oLegList[j].Speed_kn))
            #print("RoundtripTime_Days:", self.RoundtripTime_Days)
        print(df)
    
    def CalculateRoundtripTime(self, oJ):
        """
        RoundtripTime_Days = 0
        for j in range(0,oJ.LegsPerRoundtrip):
            RoundtripTime_Days = RoundtripTime_Days + self.oLegList[j].LegTime_Days
        self.RoundtripTime_Days = RoundtripTime_Days
        print("RoundtripTime_Days:" + str(RoundtripTime_Days))
        """
        pass
