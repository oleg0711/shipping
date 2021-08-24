# -*- coding: utf-8 -*-
"""
Created on Sat May 29 09:23:13 2021

@author: oleg_
"""
import numpy as np



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
    
        print("Ship data:")
        print("Vmin" + str(self.Vmin))
        print("Vmax" + str(self.Vmax))
        print("DWTscantling" + str(self.DWTscantling))
        print("DWTdesign" + str(self.DWTdesign))
        print("Lightweight" + str(self.Lightweight))
        print("k" + str(self.k))
        print("p" + str(self.p))
        print("g" + str(self.g))
        print("a" + str(self.a))
        print("ShipDischargeRate" + str(self.ShipDischargeRate))
        print("BallastCapacity" + str(self.BallastCapacity))
        print("MinFillRateShip" + str(self.MinFillRateShip))
    
    def FuelConsumption_TonnePerDay(self, Speed:float, DWT:float):
        self.FuelConsumption_TonnePerDay = self.k * (self.p + (self.Speed) ^ self.g) * (self.DWT + self.Lightweight) ^ self.a
        return self.FuelConsumption_TonnePerDay





class Cleg:

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

'Methods (functions/subs)
    def __init__(self,
                 LegNr,
                 RoundTripNr,
                 Distance_nm,
                 PortNr_StartOfLeg,
                 PortNr_EndOfLeg,
                 TimeLoading_Days,
                 TimeWaiting_Days,
                 TimeUnloading_Days,
                 CargoCarried_Tonne,
                 LoadingCost_OtherOperations_USD,
                 UnloadingCosts_USD,
                 Revenues_USD)
        
    def FindBestGoodwill(oJ:CJourney, 
                         oS:CShip, 
                         oRTList = [], 
                         oPList = []) 'list of round trips, list of ports
   
   
        'If last leg, get future goodwill from next roundtrip
        if LegNr == oJ.LegsPerRoundtrip:
            FutureGoodwill = oRTList[RoundtripNr].FutureGoodwill
            'Else get it from the next leg of this roundtrip
        else:
            FutureGoodwill = oRTList[RoundtripNr].FindGoodwillOfNextLeg(LegNr)
            
        print("RoundtripNr:" + str(RoundtripNr))
        print("LegNr:" + str(LegNr))
        print("FutureGoodwill:" + str(FutureGoodwill))
         
        'Figure out the range of Tmax and Tmin and set a stepsize
        TimeAtSeaMax_Days = Distance_nm / (24 * oS.Vmin)
        TimeAtSeaMin_Days = Distance_nm / (24 * oS.Vmax)
        '***************************************
        'CHOOSE A STEPSIZE SAY EVERY HOUR
        TimeStepSize_Days = 1 / (24 * 6)
        '***************************************
          
          
        BestGoodwill_USD:float = -999999999
        BestTimeAtSea_Days:float = 9999999
        BestSpeed_kn:float = 9999999
        MinimumLoad:float = oS.DWTdesign * oS.MinFillRateShip
        
          
        'Optimise over T range
        num:int = int(np.round(TimeAtSeaMax_Days - TimeAtSeaMin_Days) / TimeStepSize_Days,0))
        for TimeAtSea_Days in np.linspace(start=TimeAtSeaMin_Days,stop=TimeAtSeaMax_Days,num=num):
           LegTime_Days = TimeLoading_Days + TimeAtSea_Days + TimeWaiting_Days + TimeUnloading_Days
           
           AuxFuelCarried_Tonne = (TimeLoading_Days + TimeWaiting_Days + TimeUnloading_Days) * oS.AuxFuelConsumption_TonnePerDay
           LoadingCost_AuxFuel_USD = AuxFuelCarried_Tonne * oP[PortNr_StartOfLeg].AuxFuelRate_USDperTonne
           
           
           Speed_kn = Distance_nm / (24 * TimeAtSea_Days)
           
           'Minimum of Cargo carried or ballast water carried:
           DWTCarried_Tonne = CargoCarried_Tonne
           if DWTCarried_Tonne < MinimumLoad:
               DWTCarried_Tonne = MinimumLoad
           
           'Get Fuelconsumption and fuel needed using estimate for DWT carried
           MainBunkerConsumption_TonnePerDay = oS.FuelConsumption_TonnePerDay(Speed_kn, DWTCarried_Tonne)
           MainBunkerCarried_Tonne = MainBunkerConsumption_TonnePerDay * TimeAtSea_Days
           LoadingCost_MainBunker_USD = MainBunkerCarried_Tonne * oP[PortNr_StartOfLeg].MainBunkerRate_USDperTonne
           
           'Calculate goodwill and retain highest value solution
           Profitability_ThisLeg_USD = (Revenues_USD - UnloadingCosts_USD) * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365)) - LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 * (1 - np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))) / (oJ.OpprtCostCapitalRate)
           Goodwill_USD = Profitability_ThisLeg_USD + FutureGoodwill * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))
            
           if Goodwill_USD > BestGoodwill_USD:
                BestGoodwill_USD = Goodwill_USD
                BestTimeAtSea_Days = TimeAtSea_Days
                BestSpeed_kn = Speed_kn
           
           
        
        Goodwill_USD = BestGoodwill_USD
        TimeAtSea_Days = BestTimeAtSea_Days
        Speed_kn = BestSpeed_kn
        MainBunkerConsumption_TonnePerDay = oS.FuelConsumption_TonnePerDay(Speed_kn, DWTCarried_Tonne)
        MainBunkerCarried_Tonne = MainBunkerConsumption_TonnePerDay * TimeAtSea_Days
        LoadingCost_MainBunker_USD = MainBunkerCarried_Tonne * oP[PortNr_StartOfLeg].MainBunkerRate_USDperTonne
        LegTime_Days = TimeLoading_Days + TimeAtSea_Days + TimeWaiting_Days + TimeUnloading_Days
        Profitability_ThisLeg_USD = (Revenues_USD - UnloadingCosts_USD) * np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365)) - LoadingCost_MainBunker_USD - LoadingCost_AuxFuel_USD - LoadingCost_OtherOperations_USD - oJ.DailyHire_USDperDay * 365 * (1 - np.exp(-1 * oJ.OpprtCostCapitalRate * (LegTime_Days / 365))) / (oJ.OpprtCostCapitalRate)

   
        print("LoadingCost_OtherOperations_USD:" + str(LoadingCost_OtherOperations_USD))
        print("LoadingCost_MainBunker_USD:" + str(LoadingCost_MainBunker_USD))
        print("UnloadingCosts_USD" + str(UnloadingCosts_USD))

   
       'If this leg is the first of the roundtrip, store the Goodwill also into this roundtrips goodwill
       if LegNr == 1:
         oRTList[RoundtripNr].RoundtripGoodwill = Goodwill_USD
         print("oRT(RoundtripNr).RoundtripGoodwill=" str(oRTList[RoundtripNr].RoundtripGoodwill))
       
   
class CPort:
'CLASS MODULE: CPort
'Member variables
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
    LegRevenue_USD:float 'Associated with unloading
    MainBunkerRate_USDperBarrel:float
    MainBunkerRate_USDperTonne:float
    MainBunker_QBmetresperBarrel:float
    MainBunker_QbmetresPerTonne:float
    
    AuxFuelRate_USDperTonne:float
    
'Methods (functions/subs)
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
                 AuxFuelRate_USDperTonne
                 )

        if PortNr > 0: 
            self.DistancePreviousPort_nm = DistancePreviousPort_nm 'Sheet1.Cells(21, 4 + PortNr) 
        else: 
            self.DistancePreviousPort_nm = 0
        
        self.LoadingRate_QbmetresperHr = LoadingRate_QbmetresperHr 'Sheet1.Cells(22, 4 + PortNr)
        self.WaitingTime_Hrs = WaitingTime_Hrs 'Sheet1.Cells(23, 4 + PortNr)
         
        self.CargoIntake_Barrels = CargoIntake_Barrels ' Sheet1.Cells(26, 4 + PortNr)
        self.CargoIntake_QBmetresperBarrel = CargoIntake_QBmetresperBarrel 'Sheet1.Cells(27, 4 + PortNr)
        self.CargoIntake_QbmetresPerTonne = CargoIntake_QbmetresPerTonne 'Sheet1.Cells(28, 4 + PortNr)
         
        self.FixedPortAccessCosts_USD = FixedPortAccessCosts_USD 'Sheet1.Cells(41, 4 + PortNr)
        if PortNr > 0:
            self.UnloadingCharge_USDperHr = UnloadingCharge_USDperHr 'Sheet1.Cells(42, 4 + PortNr) 
        else: 
            self.UnloadingCharge_USDperHr = 0
        
        self.LoadingCharge_USDperHr = LoadingCharge_USDperHr 'Sheet1.Cells(43, 4 + PortNr)   
        self.CargoRevenueRate_USDperBarrelper1000nm = CargoRevenueRate_USDperBarrelper1000nm 'Sheet1.Cells(36, 4 + PortNr)
        self.MainBunkerRate_USDperBarrel = MainBunkerRate_USDperBarrel 'Sheet1.Cells(49, 4 + PortNr)
        self.MainBunker_QBmetresperBarrel = MainBunker_QBmetresperBarrel 'Sheet1.Cells(51, 4 + PortNr)
        self.MainBunker_QbmetresPerTonne = MainBunker_QbmetresPerTonne 'Sheet1.Cells(52, 4 + PortNr)
        self.AuxFuelRate_USDperTonne = AuxFuelRate_USDperTonne 'Sheet1.Cells(53, 4 + PortNr)
         
        
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
    
    
    
'CALCULATE BLUE CELLS
'Public Sub CalculateCargoImplicationsData(ByRef oJ As CJourney, ByRef oS As CShip)
    def CalculateCargoImplicationsData(oJ:CJourney, oS:CShip, oPList =[] ) 'list of ports

        for i in range(0,oJ.LegsPerRoundtrip):
            
            if oPList[i].CargoIntake_Barrels > 0:
                oPList[i].CargoIntake_Tonne = oPList[i].CargoIntake_Barrels * oPList[i].CargoIntake_QBmetresperBarrel / oP(i).CargoIntake_QbmetresPerTonne
            else:
                oPList[i].CargoIntake_Tonne = 0
            
            
            if i > 0:
                oPList[i].UnloadingTime_Hrs = oPList[i - 1].CargoIntake_Barrels * oP[i - ].CargoIntake_QBmetresperBarrel / oS.ShipDischargeRate
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
            
         
            print("Port: " + str(i))
            print("CargoIntake_Tonne" +str(oPList[i].CargoIntake_Tonne))
            print("UnloadingTime_Hrs" +str(oPList[i].UnloadingTime_Hrs))
            print("LoadingTime_Hrs" +str(oPList[i].LoadingTime_Hrs))
            print("UnloadingCosts_USD" +str(oPList[i].UnloadingCosts_USD))
            print("LoadingCosts_USD" +str(oPList[i].LoadingCosts_USD))
            print("LegRevenue_USD" +str(oPList[i].LegRevenue_USD))
            print("MainBunkerRate_USDperTonne" +str(oPList[i].MainBunkerRate_USDperTonne))
            print("---------------------------------------------")
        
        return oPList





class CRoundtrip():
'CLASS MODULE: CRoundtrip
'Member variables
    TripNr:int
    RoundtripGoodwill:float
    oLegList:list
    FutureGoodwill:float
    RoundtripTime_Days:float
    oJ:CJorney
    oS:CShip
    oPList:list
        
    def __init__(self, oJ:CJourney, oS:CShip, oPList=[], TripNr)
        
        self.oJ = oJ
        self.oS = oS
        self.oPList = oPLst
        
        SPLIT = 1
        
        For i in range(1,1+oJ.LegsPerRoundtrip):
            
            oLeg = CLeg(LegNr=LegNr,
                        RoundTripNr=TripNr,
                        Distance_nm = oPList[i].DistancePreviousPort_nm
                        PortNr_StartOfLeg = i - 1,
                        PortNr_EndOfLeg = i,
                        TimeLoading_Days = oPList[i - 1].LoadingTime_Hrs / 24,
                        TimeWaiting_Days = oPList[i].WaitingTime_Hrs / 24,
                        TimeUnloading_Days = oPList[i].UnloadingTime_Hrs / 24,
                        CargoCarried_Tonne = oPList[i-1].CargoIntake_Tonne,
                        LoadingCost_OtherOperations_USD = oP(i - 1).LoadingCosts_USD + SPLIT * (oP(i - 1).FixedPortAccessCosts_USD),
                        UnloadingCosts_USD = oP(i).UnloadingCosts_USD + (1 - SPLIT) * (oP(i).FixedPortAccessCosts_USD),
                        Revenues_USD = oP(i).LegRevenue_USD)
                        
            self.oLegList.append(oLeg)
            
            print("Leg data: " + str(i))
            print("LegNr" + str(oLeg.LegNr))
            print("RoundtripNr" + str(oLeg.RoundtripNr))
            print("Distance_nm" + str(oLeg.Distance_nm))
            print("PortNr_StartOfLeg" + str(oLeg.PortNr_StartOfLeg))
            print("PortNr_EndOfLeg" + str(oLeg.PortNr_EndOfLeg))
            print("TimeLoading_Days" + str(oLeg.TimeLoading_Days))
            print("TimeWaiting_Days" + str(oLeg.TimeWaiting_Days))
            print("TimeUnloading_Days" + str(oLeg.TimeUnloading_Days))
            print("CargoCarried_Tonne" + str(oLeg.CargoCarried_Tonne
            print("LoadingCost_OtherOperations_USD" + str(oLeg.LoadingCost_OtherOperations_USD))
            print("UnloadingCosts_USD" + str(oLeg.UnloadingCosts_USD))
            print("Revenues_USD" + str(oLeg.Revenues_USD))


    def FindBestGoodwill(self, oRTList)
        for j in range(self.oJ.LegsPerRoundtrip, 0, -1):
            oLegList[j].FindBestGoodwill(self.oJ, self.oS, oRTList, self.oP)
    
    
    def FindGoodwillOfNextLeg(self,CurrentLegNr:int):
        self.FindGoodwillOfNextLeg = self.oLegList[CurrentLegNr + 1].Goodwill_USD


    def PrintOptimalLegSpeeds(self):
        for j in range(1,1+self.oJ.LegsPerRoundtrip):
            print("Leg Nr: " + str(j))
            print("Days at Sea (Days):" + str(self.oLegList[j].TimeAtSea_Days))
            print("Speed(knots):" + str(self.oLegList[j].Speed_kn))
            print("RoundtripTime_Days:", self.RoundtripTime_Days)
   
    
    def CalculateRoundtripTime(self):
        
        RoundtripTime_Days = 0
        for j in range(1,1+self.oJ.LegsPerRoundtrip):
            RoundtripTime_Days = RoundtripTime_Days + oLegList[j].LegTime_Days
        
        self.RoundtripTime_Days = RoundtripTime_Days
        print("RoundtripTime_Days" + str(RoundtripTime_Days))
    


class CJourney:
    
    'CLASS MODULE: CJourney
'Member variables

    NrOfRoundtrips:int
    LegsPerRoundtrip:int
    OpprtCostCapitalRate:float
    DailyHire_USDperDay:float
    FutureProfitPotential_USDperDay:float
    
    TotalJourneyTime_Days:float
    NPV_Profits_IncludingFPP_USD:float
    TCE_Journey_USDperDay:float
    JourneyProfits_USDperDay:float


'Methods (functions/subs)
    def __init__(self,
                 NrOfRoundtrips,
                 LegsPerRoundtrip,
                 OpprtCostCapitalRate,
                 DailyHire_USDperDay,
                 FutureProfitPotential_USDperDay):
        
        self.NrOfRoundtrips = NrOfRoundtrips 'Sheet1.Cells(3, 3)
        self.LegsPerRoundtrip = LegsPerRoundtrip'Sheet1.Cells(4, 3)
        self.OpprtCostCapitalRate = OpprtCostCapitalRate'Sheet1.Cells(5, 3)
        self.DailyHire_USDperDay = DailyHire_USDperDay'Sheet1.Cells(6, 3)
        self.FutureProfitPotential_USDperDay = FutureProfitPotential_USDperDay'Sheet1.Cells(7, 3)
    
        print ("NrOfRoundtrips": + str(NrOfRoundtrips))
        print ("LegsPerRoundtrip": + str(LegsPerRoundtrip))
        print ("OpprtCostCapitalRate": + str(OpprtCostCapitalRate))
        print ("DailyHire_USDperDay": + str(DailyHire_USDperDay))
        print ("FutureProfitPotential_USDperDay": + str(FutureProfitPotential_USDperDay))
        
    
    def CalculateJourneyTime(oRTList:list):
        
        TotalJourneyTime_Days = 0
        for i in range(1,1+NrOfRoundtrips):
            TotalJourneyTime_Days = TotalJourneyTime_Days + oRTList[i].RoundtripTime_Days
        self.TotalJourneyTime_Days = TotalJourneyTime_Days

    
    def CalculateProfits(oRTList:list):

        'NPV profits:
        self.NPV_Profits_IncludingFPP_USD = oRTList[NrOfRoundtrips].RoundtripGoodwill
        'TCE EXCLUDING Daily Hire AND EXCLUDING Future Profit Potential
        'This calculates the AS value (as a USD/year rate) received over the TotalJourneyTime_Days:
        self.TCE_Journey_USDperDay = (self.NPV_Profits_IncludingFPP_USD - (self.FutureProfitPotential_USDperDay * 365 / self.OpprtCostCapitalRate) * np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365)) * self.OpprtCostCapitalRate / (1 - np.exp(-1 * self.OpprtCostCapitalRate * self.TotalJourneyTime_Days / 365))
        self.TCE_Journey_USDperDay = self.TCE_Journey_USDperDay / 365
        
        'This calculates the profit over the journey as a per day measure:
        self.JourneyProfits_USDperDay = self.TCE_Journey_USDperDay
        'This adjusts by adding back the Daily Hire (TCH), as not incurred in TCE calculation
        self.TCE_Journey_USDperDay = self.TCE_Journey_USDperDay + self.DailyHire_USDperDay
        

    

