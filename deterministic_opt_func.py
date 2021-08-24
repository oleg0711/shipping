# -*- coding: utf-8 -*-
"""

Created on Sat May 29 08:13:59 2021
@author: Patrick (VBA Code), Oleg (Python, reprogramming the excel spreadsheet in Python)
@author: Oleg (adaptation of the deterministic model to the stochastic setting)

"""


from Logistics import CShip
from Logistics import CRoundtrip
from Logistics import CJourney
from Logistics import CPort




def run_deterministic_opt(oJ,oS,oPL):


    oJourney = CJourney(NrOfRoundtrips = oJ.NrOfRoundtrips,
                        LegsPerRoundtrip = oJ.LegsPerRoundtrip,
                        OpprtCostCapitalRate = oJ.OpprtCostCapitalRate,
                        DailyHire_USDperDay = oJ.DailyHire_USDperDay,
                        FutureProfitPotential_USDperDay = oJ.FutureProfitPotential_USDperDay)
    
    
    oShip = CShip(Vmin=oS.Vmin, 
                  Vmax=oS.Vmax, 
                  DWTscantling=oS.DWTscantling, 
                  DWTdesign=oS.DWTdesign, 
                  Lightweight=oS.Lightweight, 
                  k=oS.k, 
                  p=oS.p, 
                  g=oS.g, 
                  a=oS.a, 
                  ShipDischargeRate=oS.ShipDischargeRate, 
                  BallastCapacity=oS.BallastCapacity,
                  MinFillRateShip=oS.MinFillRateShip,
                  AuxFuelConsumption_TonnePerDay=oS.AuxFuelConsumption_TonnePerDay)
    
    oPortList = []
    
    for i in range(0,oJourney.LegsPerRoundtrip+1):
        
        oPort = CPort(PortNr=i, 
                      DistancePreviousPort_nm=oPL[i].DistancePreviousPort_nm,
                      LoadingRate_QbmetresperHr=oPL[i].LoadingRate_QbmetresperHr,
                      WaitingTime_Hrs=oPL[i].WaitingTime_Hrs,              
                      CargoIntake_Barrels=oPL[i].CargoIntake_Barrels,
                      CargoIntake_QBmetresperBarrel=oPL[i].CargoIntake_QBmetresperBarrel,
                      CargoIntake_QbmetresPerTonne=oPL[i].CargoIntake_QbmetresPerTonne,
                      FixedPortAccessCosts_USD = oPL[i].FixedPortAccessCosts_USD,
                      UnloadingCharge_USDperHr = oPL[i].UnloadingCharge_USDperHr,
                      LoadingCharge_USDperHr = oPL[i].LoadingCharge_USDperHr,
                      CargoRevenueRate_USDperBarrelper1000nm = oPL[i].CargoRevenueRate_USDperBarrelper1000nm,
                      MainBunkerRate_USDperBarrel = oPL[i].MainBunkerRate_USDperBarrel,
                      MainBunker_QBmetresperBarrel = oPL[i].MainBunker_QBmetresperBarrel,
                      MainBunker_QbmetresPerTonne = oPL[i].MainBunker_QbmetresPerTonne,
                      AuxFuelRate_USDperTonne = oPL[i].AuxFuelRate_USDperTonne,
                      UnloadingCosts_USD = oPL[i].UnloadingCosts_USD,
                      LoadingCosts_USD = oPL[i].LoadingCosts_USD,
                      CargoIntake_Tonne = oPL[i].CargoIntake_Tonne,
                      UnloadingTime_Hrs = oPL[i].UnloadingTime_Hrs,
                      LoadingTime_Hrs = oPL[i].LoadingTime_Hrs,
                      LegRevenue_USD = oPL[i].LegRevenue_USD)
        oPortList.append(oPort)
    
        
    oPortList[0].CalculateCargoImplicationsData(oJ=oJourney, oS=oShip, oPList=oPortList)
    
    oRTList = []
    
    TotalJourneyDistance_nm = 0
    for i in range(1, 1 + oJourney.NrOfRoundtrips):
        rt = CRoundtrip(oJourney, oShip, oPortList, TripNr=i)
        oRTList.append(rt)
    
    oJourney.TotalJourneyDistance_nm = TotalJourneyDistance_nm
    
    rt0 = CRoundtrip(oJourney, oShip, oPortList, TripNr=0)
    rt0.RoundtripGoodwill = oJourney.FutureProfitPotential_USDperDay * 365 / oJourney.OpprtCostCapitalRate
    oRTList.insert(0,rt0)
    oRTList[1].FutureGoodwill = rt0.RoundtripGoodwill
    
    #For each roundtrip, for each leg ... find the best goodwill
    for i in range(1,1+oJourney.NrOfRoundtrips):
       oRTList[i].FindBestGoodwill(oRTList)
       if i < oJourney.NrOfRoundtrips:
           oRTList[i + 1].FutureGoodwill = oRTList[i].RoundtripGoodwill
    
    
    #Calculate solution metrics
    for i in range(oJourney.NrOfRoundtrips,0,-1):
       oRTList[i].CalculateRoundtripTime(oJourney)
    #print(oRTList[0].RoundtripTime_Days)
    
    oJourney.CalculateJourneyTime(oRTList)
    oJourney.CalculateProfits(oRTList)
    
    
    #Print solution
    #print("NPV Future Profit Potential (USD):" +str(oRTList[0].RoundtripGoodwill))
    #print("NPV all future operations of ship (USD):" +str( oRTList[oJourney.NrOfRoundtrips].RoundtripGoodwill))
    #print("AS all future operations of ship (USD):" +str(oRTList[oJourney.NrOfRoundtrips].RoundtripGoodwill * oJourney.OpprtCostCapitalRate))
    
    #print("oJourney.TotalJourneyTime_Days:" +str(oJourney.TotalJourneyTime_Days))
    #print("oJourney.NPV_Profits_IncludingFPP_USD:"+str(oJourney.NPV_Profits_IncludingFPP_USD))
    #print("oJourney.TCE_Journey_USDperDay:"+str(oJourney.TCE_Journey_USDperDay))
    #print("oJourney.DailyHire_USDperDay:"+str(oJourney.DailyHire_USDperDay))
    #print("oJourney.FutureProfitPotential_USDperDay:"+str(oJourney.FutureProfitPotential_USDperDay))
    #print("oJourney.JourneyProfits_USDperDay:"+str(oJourney.JourneyProfits_USDperDay))
    
    for i in range(oJourney.NrOfRoundtrips,0,-1):
    #  print("Roundtrip Journey Nr: " + str(oRTList[i].TripNr))
    #  print("RoundtripTime_Days:" +  str(oRTList[i].RoundtripTime_Days))
       oRTList[i].PrintOptimalLegSpeeds(oJourney)
    
    return oRTList  


