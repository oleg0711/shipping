# -*- coding: utf-8 -*-
"""
Created on Sat May 29 09:23:13 2021

@author: oleg_
"""
import numpy as np
import pandas as pd


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
        """
        print("Ship data-------------------")
        print("Vmin:" + str(self.Vmin))
        print("Vmax:" + str(self.Vmax))
        print("DWTscantling:" + str(self.DWTscantling))
        print("DWTdesign:" + str(self.DWTdesign))
        print("Lightweight:" + str(self.Lightweight))
        print("k:" + str(self.k))
        print("p:" + str(self.p))
        print("g:" + str(self.g))
        print("a:" + str(self.a))
        print("ShipDischargeRate:" + str(self.ShipDischargeRate))
        print("BallastCapacity:" + str(self.BallastCapacity))
        print("MinFillRateShip:" + str(self.MinFillRateShip))
        """
    def FuelConsumption_TonnePerDay(self, Speed:float, DWT:float):
        FuelConsumption_TonnePerDay = self.k * (self.p + (Speed) ** self.g) * (DWT + self.Lightweight) ** self.a
        return FuelConsumption_TonnePerDay





class CJourney:
    
    NrOfRoundtrips:int
    LegsPerRoundtrip:int
    OpprtCostCapitalRate:float
    DailyHire_USDperDay:float
    FutureProfitPotential_USDperDay:float
    
    TotalJourneyTime_Days:float
    NPV_Profits_IncludingFPP_USD:float
    TCE_Journey_USDperDay:float
    JourneyProfits_USDperDay:float

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
        """
        print ("NrOfRoundtrips:" + str(NrOfRoundtrips))
        print ("LegsPerRoundtrip:" + str(LegsPerRoundtrip))
        print ("OpprtCostCapitalRate:" + str(OpprtCostCapitalRate))
        print ("DailyHire_USDperDay:" + str(DailyHire_USDperDay))
        print ("FutureProfitPotential_USDperDay:" + str(FutureProfitPotential_USDperDay))
        """
    
    def CalculateJourneyTime(self, oRTList:list):
        
        TotalJourneyTime_Days = 0
        for i in range(1,1+self.NrOfRoundtrips):
            TotalJourneyTime_Days = TotalJourneyTime_Days + oRTList[i].RoundtripTime_Days
        self.TotalJourneyTime_Days = TotalJourneyTime_Days
        return TotalJourneyTime_Days

    
    def CalculateProfits(self, oRTList:list):

        #NPV profits:
        self.NPV_Profits_IncludingFPP_USD = oRTList[self.NrOfRoundtrips].RoundtripGoodwill
        #TCE EXCLUDING Daily Hire AND EXCLUDING Future Profit Potential
        #This calculates the AS value (as a USD/year rate) received over the TotalJourneyTime_Days:
        self.TCE_Journey_USDperDay = (self.NPV_Profits_IncludingFPP_USD - (self.FutureProfitPotential_USDperDay * 365 / self.OpprtCostCapitalRate) * np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365)) * self.OpprtCostCapitalRate / (1 - np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365))
        self.TCE_Journey_USDperDay = self.TCE_Journey_USDperDay / 365
        
        #This calculates the profit over the journey as a per day measure:
        self.JourneyProfits_USDperDay = self.TCE_Journey_USDperDay
        #This adjusts by adding back the Daily Hire (TCH), as not incurred in TCE calculation
        self.TCE_Journey_USDperDay = self.TCE_Journey_USDperDay + self.DailyHire_USDperDay
        
   
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
                 Revenues_USD):
    
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
        self.Revenues_USD = Revenues_USD
    
    def CalcGoodwill(self,
                         TimeAtSea_Days,
                         oJ:CJourney, 
                         oS:CShip, 
                         oRTList = [], 
                         oPList = []):
        
        #If last leg, get future goodwill from next roundtrip
        if self.LegNr == oJ.LegsPerRoundtrip:
            self.FutureGoodwill = oRTList[self.RoundtripNr].FutureGoodwill
            #Else get it from the next leg of this roundtrip
        else:
            self.FutureGoodwill = oRTList[self.RoundtripNr].FindGoodwillOfNextLeg(self.LegNr)
            
        print("RoundtripNr:" + str(self.RoundtripNr))
        print("LegNr:" + str(self.LegNr))
        print("FutureGoodwill:" + str(self.FutureGoodwill))
         
        #Figure out the range of Tmax and Tmin and set a stepsize
        self.TimeAtSeaMax_Days = self.Distance_nm / (24 * oS.Vmin)
        self.TimeAtSeaMin_Days = self.Distance_nm / (24 * oS.Vmax)
        #***************************************
        #CHOOSE A STEPSIZE SAY EVERY HOUR
        self.TimeStepSize_Days = np.min([1 / (24 * 6), np.max([1/10000,oJ.OpprtCostCapitalRate])])
        #***************************************
          
          
        self.BestGoodwill_USD:float = -999999999
        self.BestTimeAtSea_Days:float = 9999999
        self.BestSpeed_kn:float = 9999999
        self.MinimumLoad:float = oS.DWTdesign * oS.MinFillRateShip
        
          
        #Optimise over T range
        num:int = int(np.round((self.TimeAtSeaMax_Days - self.TimeAtSeaMin_Days) / self.TimeStepSize_Days,0))
        
        
        print(str(self.LegNr)+'-!-!-!-!----Leg-------------------------!!!!!!')
        print(str(self.RoundtripNr)+'-!-!-!-!----RoundTrip-------------------------!!!!!!')
        
        
          
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
        Profitability_ThisLeg_USD = (self.Revenues_USD - self.UnloadingCosts_USD) * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365)) - self.LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - self.LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 * (1 - np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))) / (oJ.OpprtCostCapitalRate)

        Goodwill_USD = Profitability_ThisLeg_USD + self.FutureGoodwill * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))        
        return Goodwill_USD,Profitability_ThisLeg_USD
            
    def FindBestGoodwill(self,
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
        print("FutureGoodwill:" + str(self.FutureGoodwill))
         
        #Figure out the range of Tmax and Tmin and set a stepsize
        self.TimeAtSeaMax_Days = self.Distance_nm / (24 * oS.Vmin)
        self.TimeAtSeaMin_Days = self.Distance_nm / (24 * oS.Vmax)
        #***************************************
        #CHOOSE A STEPSIZE SAY EVERY HOUR
        self.TimeStepSize_Days = np.min([1 / (24 * 6), np.max([1/10000,oJ.OpprtCostCapitalRate])])
        #***************************************
          
          
        self.BestGoodwill_USD:float = -999999999
        self.BestTimeAtSea_Days:float = 9999999
        self.BestSpeed_kn:float = 9999999
        self.MinimumLoad:float = oS.DWTdesign * oS.MinFillRateShip
        
          
        #Optimise over T range
        num:int = 5*int(np.round((self.TimeAtSeaMax_Days - self.TimeAtSeaMin_Days) / self.TimeStepSize_Days,0))
        
        
        print(str(self.LegNr)+'-!-!-!-!----Leg-------------------------!!!!!!')
        print(str(self.RoundtripNr)+'-!-!-!-!----RoundTrip-------------------------!!!!!!')
        
        
        iteration_nr = 0
        
        for TimeAtSea_Days in np.linspace(start=self.TimeAtSeaMin_Days,stop=self.TimeAtSeaMax_Days,num=num):
           
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
            Profitability_ThisLeg_USD = (self.Revenues_USD - self.UnloadingCosts_USD) * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365)) - self.LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - self.LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 * (1 - np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))) / (oJ.OpprtCostCapitalRate)
            
            Goodwill_USD = Profitability_ThisLeg_USD + self.FutureGoodwill * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))

            
                
                
            """
            if iteration_nr==0:
                print("---------------------------------")
                print("MinimumLoad:" +str(self.MinimumLoad))
                print("CargoCarried_Tonne:" +str(self.CargoCarried_Tonne))
                
                print("self.TimeLoading_Days:" + str(self.TimeLoading_Days))
                print("TimeAtSea_Days:" + str(TimeAtSea_Days))
                print("self.TimeWaiting_Days:" +str(self.TimeWaiting_Days))
                print("self.TimeUnloading_Days:" + str(self.TimeUnloading_Days))
               
                print("LegTime_Days:" +str(LegTime_Days))
                print("AuxFuelCarried_Tonne:" +str(AuxFuelCarried_Tonne))
                print("LoadingCost_AuxFuel_USD:" +str(LoadingCost_AuxFuel_USD))
                print("Speed_kn:" +str(self.Speed_kn))
               
                print("DWTCarried_Tonne:" +str(self.DWTCarried_Tonne))
                print("MainBunkerConsumption_TonnePerDay:" +str(self.MainBunkerConsumption_TonnePerDay))
                print("MainBunkerCarried_Tonne:" +str(self.MainBunkerCarried_Tonne))
               
                print("LoadingCost_MainBunker_USD:" +str(self.LoadingCost_MainBunker_USD))
                
                print("Revenues_USD" + str(self.Revenues_USD))
                print("UnloadingCosts_USD" + str(self.UnloadingCosts_USD))
                print("oJ.OpprtCostCapitalRate" + str(oJ.OpprtCostCapitalRate))
                
                
                print("Profitability_ThisLeg_USD:" +str(Profitability_ThisLeg_USD))
                
                print("Goodwill_USD:" +str(Goodwill_USD))
                print("---------------------------------")
                
                
               
                iteration_nr=1
           
           """
           
            if Goodwill_USD > self.BestGoodwill_USD:
                 self.BestGoodwill_USD = Goodwill_USD
                 self.BestTimeAtSea_Days = TimeAtSea_Days
                 self.BestSpeed_kn = self.Speed_kn
                 self.BestProfitability_ThisLeg_USD = Profitability_ThisLeg_USD
           
           
        self.Goodwill_USD = self.BestGoodwill_USD
        self.TimeAtSea_Days = self.BestTimeAtSea_Days
        self.Speed_kn = self.BestSpeed_kn
        self.MainBunkerConsumption_TonnePerDay = oS.FuelConsumption_TonnePerDay(self.Speed_kn, self.DWTCarried_Tonne)
        self.MainBunkerCarried_Tonne = self.MainBunkerConsumption_TonnePerDay * TimeAtSea_Days
        self.LoadingCost_MainBunker_USD = self.MainBunkerCarried_Tonne * oPList[self.PortNr_StartOfLeg].MainBunkerRate_USDperTonne
        self.LegTime_Days = self.TimeLoading_Days + self.TimeAtSea_Days + self.TimeWaiting_Days + self.TimeUnloading_Days
        self.Profitability_ThisLeg_USD = self.BestProfitability_ThisLeg_USD
#        (self.Revenues_USD - self.UnloadingCosts_USD) * np.exp(-1 * oJ.OpprtCostCapitalRate * (self.LegTime_Days / 365)) - #self.LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - self.LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 #* (1 - np.exp(-1 * oJ.OpprtCostCapitalRate * (self.LegTime_Days / 365))) / (oJ.OpprtCostCapitalRate)

   
        print("LoadingCost_OtherOperations_USD:" + str(self.LoadingCost_OtherOperations_USD))
        print("LoadingCost_MainBunker_USD:" + str(self.LoadingCost_MainBunker_USD))
        print("UnloadingCosts_USD:" + str(self.UnloadingCosts_USD))

   
       #If this leg is the first of the roundtrip, store the Goodwill also into this roundtrips goodwill
        if self.LegNr == 1:
            oRTList[self.RoundtripNr].RoundtripGoodwill = self.Goodwill_USD
            print("oRT(RoundtripNr).RoundtripGoodwill=:" + str(oRTList[self.RoundtripNr].RoundtripGoodwill))
       
   
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
        self.LegRevenue_USD = LegRevenue_USD    
        """
        print("Port data:" + str(PortNr))
        print("DistancePreviousPort_nm:" + str(self.DistancePreviousPort_nm))
        
        print("LoadingRate_QbmetresperHr:" + str(self.LoadingRate_QbmetresperHr))
        print("WaitingTime_Hrs:" + str(self.LoadingRate_QbmetresperHr))
        
        print("CargoIntake_Barrels:" + str(self.CargoIntake_Barrels))
        print("CargoIntake_QBmetresperBarrel:" + str(self.CargoIntake_QBmetresperBarrel))
        print("CargoIntake_QbmetresPerTonne:" + str(self.CargoIntake_QbmetresPerTonne))
          
        print("FixedPortAccessCosts_USD:" + str(self.FixedPortAccessCosts_USD))
        print("UnloadingCharge_USDperHr:" + str(self.UnloadingCharge_USDperHr))
        print("LoadingCharge_USDperHr:" + str(self.LoadingCharge_USDperHr))
        
        print("CargoRevenueRate_USDperBarrelper1000nm:" + str(self.CargoRevenueRate_USDperBarrelper1000nm))
        print("MainBunkerRate_USDperBarrel:" + str(self.MainBunkerRate_USDperBarrel))
        print("MainBunker_QBmetresperBarrel:" + str(self.MainBunker_QBmetresperBarrel))
        print("MainBunker_QbmetresPerTonne:" + str(self.MainBunker_QbmetresPerTonne))
        
        print("AuxFuelRate_USDperTonne:" + str(self.AuxFuelRate_USDperTonne))
        """
    
    
    
    def CalculateCargoImplicationsData(self, oJ, oS, oPList): #list of ports
        print ("********* Running CalculateCargoImplicationsData")
        for i in range(0,1+oJ.LegsPerRoundtrip):
            print("i="+str(i))
            print("oPList[i].CargoIntake_Barrels="+str(oPList[i].CargoIntake_Barrels))
            
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
                oPList[i].LegRevenue_USD = oPList[i - 1].CargoRevenueRate_USDperBarrelper1000nm * (oPList[i].DistancePreviousPort_nm / 1000) * oPList[i - 1].CargoIntake_Barrels
            else:
               oPList[i].LegRevenue_USD = 0
            
         
            if i < oJ.LegsPerRoundtrip:
               oPList[i].MainBunkerRate_USDperTonne = (oPList[i].MainBunkerRate_USDperBarrel / oPList[i].MainBunker_QBmetresperBarrel) * oPList[i].MainBunker_QbmetresPerTonne
            else:
               oPList[i].MainBunkerRate_USDperTonne = 0
            
            """
            print("Port: " + str(i))
            print("CargoIntake_Tonne:" +str(oPList[i].CargoIntake_Tonne))
            print("UnloadingTime_Hrs:" +str(oPList[i].UnloadingTime_Hrs))
            print("LoadingTime_Hrs:" +str(oPList[i].LoadingTime_Hrs))
            print("UnloadingCosts_USD:" +str(oPList[i].UnloadingCosts_USD))
            print("LoadingCosts_USD:" +str(oPList[i].LoadingCosts_USD))
            print("LegRevenue_USD:" +str(oPList[i].LegRevenue_USD))
            print("MainBunkerRate_USDperTonne:" +str(oPList[i].MainBunkerRate_USDperTonne))
            print("---------------------------------------------")
            """
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
        
    def __init__(self, oJ:CJourney, oS:CShip, oPList, TripNr):
        
        self.oJ = oJ
        self.oS = oS
        self.oPList = oPList
        self.oLegList = []
        self.RoundtripTime_Days = 0
        self.TripNr = TripNr
        
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
                        Revenues_USD = oPList[i].LegRevenue_USD)
                        
            self.oLegList.append(oLeg)
            """
            print("Leg data: " + str(i))
            print("LegNr:" + str(oLeg.LegNr))
            print("RoundtripNr:" + str(oLeg.RoundtripNr))
            print("Distance_nm:" + str(oLeg.Distance_nm))
            print("PortNr_StartOfLeg:" + str(oLeg.PortNr_StartOfLeg))
            print("PortNr_EndOfLeg:" + str(oLeg.PortNr_EndOfLeg))
            print("TimeLoading_Days:" + str(oLeg.TimeLoading_Days))
            print("TimeWaiting_Days:" + str(oLeg.TimeWaiting_Days))
            print("TimeUnloading_Days:" + str(oLeg.TimeUnloading_Days))
            print("CargoCarried_Tonne:" + str(oLeg.CargoCarried_Tonne))
            print("LoadingCost_OtherOperations_USD:" + str(oLeg.LoadingCost_OtherOperations_USD))
            print("UnloadingCosts_USD:" + str(oLeg.UnloadingCosts_USD))
            print("Revenues_USD:" + str(oLeg.Revenues_USD))
            """

    def FindBestGoodwill(self, oRTList):
        for j in range(self.oJ.LegsPerRoundtrip-1, -1, -1):
            self.oLegList[j].FindBestGoodwill(self.oJ, self.oS, oRTList, self.oPList)
    
    
    def FindGoodwillOfNextLeg(self,CurrentLegNr:int):
        #print("Size of leg list:" + str(len(self.oLegList)))
        #print("CurrentLegNr" + str(CurrentLegNr+1))
        return self.oLegList[CurrentLegNr].Goodwill_USD

    def PrintOptimalLegSpeeds(self, oJ):
        df = pd.DataFrame(columns=['Leg Nr','Days at Sea','Speed(knots)','Roundtrip Time Days'],dtype=object)
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
        
        RoundtripTime_Days = 0
        for j in range(0,oJ.LegsPerRoundtrip):
            print("Delta T = "+str(self.oLegList[j].LegTime_Days))
            RoundtripTime_Days = RoundtripTime_Days + self.oLegList[j].LegTime_Days
        
        self.RoundtripTime_Days = RoundtripTime_Days
        print("RoundtripTime_Days:" + str(RoundtripTime_Days))
    
