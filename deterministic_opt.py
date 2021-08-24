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


print("NEW RUN ******")

oJourney = CJourney(NrOfRoundtrips = 1,
                    LegsPerRoundtrip = 4,
                    OpprtCostCapitalRate = 0.08,
                    DailyHire_USDperDay = 30000,
                    FutureProfitPotential_USDperDay = 12968)


oShip = CShip(Vmin=10, 
              Vmax=17, 
              DWTscantling=157880, 
              DWTdesign=145900, 
              Lightweight=49000, 
              k=0.00000391, 
              p=381, 
              g=3.1, 
              a=0.666667, 
              ShipDischargeRate=3000, 
              BallastCapacity=54500,
              MinFillRateShip=0.3,
              AuxFuelConsumption_TonnePerDay=5)

oPortList = []

df_data = pd.read_excel('journey_data.xlsx',sheet_name='Ports')
for i in range(0,oJourney.LegsPerRoundtrip+1):
    for index,row in df_data.iterrows():
        v_name = row['Variable']
        v_value = row[i]
        print("Loading from excel, leg:"+str(i)+", Executing statement("+str(v_name)+"="+str(v_value)+")")
        
        exec(v_name + '=' + str(v_value))
    
    
    
    oPort = CPort(PortNr=i, 
                  DistancePreviousPort_nm=DistancePreviousPort_nm,
                  LoadingRate_QbmetresperHr=LoadingRate_QbmetresperHr,
                  WaitingTime_Hrs=WaitingTime_Hrs,              
                  CargoIntake_Barrels=CargoIntake_Barrels,
                  CargoIntake_QBmetresperBarrel=CargoIntake_QBmetresperBarrel,
                  CargoIntake_QbmetresPerTonne=CargoIntake_QbmetresPerTonne,
                  FixedPortAccessCosts_USD = FixedPortAccessCosts_USD,
                  UnloadingCharge_USDperHr = UnloadingCharge_USDperHr,
                  LoadingCharge_USDperHr = LoadingCharge_USDperHr,
                  CargoRevenueRate_USDperBarrelper1000nm = CargoRevenueRate_USDperBarrelper1000nm,
                  MainBunkerRate_USDperBarrel = MainBunkerRate_USDperBarrel,
                  MainBunker_QBmetresperBarrel = MainBunker_QBmetresperBarrel,
                  MainBunker_QbmetresPerTonne = MainBunker_QbmetresPerTonne,
                  AuxFuelRate_USDperTonne = AuxFuelRate_USDperTonne,
                  UnloadingCosts_USD = UnloadingCosts_USD,
                  LoadingCosts_USD = LoadingCosts_USD,
                  CargoIntake_Tonne = CargoIntake_Tonne,
                  UnloadingTime_Hrs = UnloadingTime_Hrs,
                  LoadingTime_Hrs = LoadingTime_Hrs,
                  LegRevenue_USD = LegRevenue_USD)
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
print(oRTList[0].RoundtripTime_Days)

oJourney.CalculateJourneyTime(oRTList)
oJourney.CalculateProfits(oRTList)


#Print solution
print("NPV Future Profit Potential (USD):" +str(oRTList[0].RoundtripGoodwill))
print("NPV all future operations of ship (USD):" +str( oRTList[oJourney.NrOfRoundtrips].RoundtripGoodwill))
print("AS all future operations of ship (USD):" +str(oRTList[oJourney.NrOfRoundtrips].RoundtripGoodwill * oJourney.OpprtCostCapitalRate))

print("oJourney.TotalJourneyTime_Days:" +str(oJourney.TotalJourneyTime_Days))
print("oJourney.NPV_Profits_IncludingFPP_USD:"+str(oJourney.NPV_Profits_IncludingFPP_USD))
print("oJourney.TCE_Journey_USDperDay:"+str(oJourney.TCE_Journey_USDperDay))
print("oJourney.DailyHire_USDperDay:"+str(oJourney.DailyHire_USDperDay))
print("oJourney.FutureProfitPotential_USDperDay:"+str(oJourney.FutureProfitPotential_USDperDay))
print("oJourney.JourneyProfits_USDperDay:"+str(oJourney.JourneyProfits_USDperDay))

for i in range(oJourney.NrOfRoundtrips,0,-1):
  print("Roundtrip Journey Nr: " + str(oRTList[i].TripNr))
  print("RoundtripTime_Days:" +  str(oRTList[i].RoundtripTime_Days))
  oRTList[i].PrintOptimalLegSpeeds(oJourney)



