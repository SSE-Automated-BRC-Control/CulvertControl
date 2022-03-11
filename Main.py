"""
Created on Fri Mar 11 09:51:14 2022

@author: nettiewallace
"""
import numpy as np

#Biomass Variables
pm = 25 #Maximum net daily photosynthesis rate at 20°C (gCO2/m2day)
kco=0.65 #conversion of carbon dioxide to AGB (g/g)
o=1.09 #Temperature constant 
tbase=5 #Base air temperature (C)
Knp=1  #Nutrient availability , maximum of 1 when no nutrient limitation 
IPAR=5.8 # Daily total incident photosynthetically active radiation (moles/m2day)
Ya=0.0025 #Leaf and shoot mortality constant 

#Data file containing [day;percip;ETo;Temp;Radiation Data]
WeatherData = np.loadtxt(fname = 'Weather_data.csv', delimiter=',')

day= WeatherData[119:275,0]
precip= WeatherData[119:275,1] #mm
Temp= WeatherData[119:275,3] 
RAD=WeatherData[119:275,4]

def biomass (pm,kco,o,tbase,Knp,IPAR,Ya):
  
# This funciton predicts the quantity of biomass produced 
# The cumulative biomass produced is the total biomass  
# The daily biomass production is the AGB_Daily 
  
#Growing Degree Days  
  GDD = [] 
  for num in range(len(day)):
    if Temp[num]-tbase <0:
        GDD=np.append(GDD,0)
    else:
        GDD=np.append(GDD,Temp[num]-tbase)
  GDD_cumulative = np.cumsum(GDD)
    
#LAI;
  LAI=[]
  for num in range(len(day)):
    x=day[num]
    if GDD_cumulative[num]>0:
      LAI=np.append(LAI,((-0.0002000*(x**2))+(0.08000*x)-5.000))
    if LAI[num]<=1:
        LAI=np.append(LAI,0)
        
#Changes Radiation Data to PAR data 
  PAR=[0]*len(day)
  for num in range(len(RAD)):
    PAR[num]=RAD[num]*4.57*0.4*3600/1000000
    PAR[num]=RAD[num]/(IPAR+RAD[num])

#Biomass calcuation 
  AGB_Daily = [] #This will hold the daily production of above ground biomass (g/m2)
  for num in range(len(day)):
    if GDD[num]>0:
      AGB_Daily=np.append(AGB_Daily,(Knp*pm*kco*PAR[num]*(o**(Temp[num]-20))*LAI[num]))
    else: 
      AGB_Daily=np.append(AGB_Daily,0)
  
#Mortality Rate of Leaves and Shoots
  mortality=[]
  for num in range(len(day)):
      mortality=np.append(mortality,(Ya*AGB_Daily[num]*(o**(Temp[num]-20))))

  AGB_cumulative = []
  mortality_cumulative = np.cumsum(mortality); #cumulative mortality 
  AGB_cumulative = np.cumsum(AGB_Daily); #Biomass production without mortality (g/m2)
  totalbiomass = AGB_cumulative-mortality_cumulative; #Total cumulative biomass (g/m2)
    
  return(AGB_Daily,LAI,PAR,RAD,AGB_cumulative)

AGB_Daily, lai, par, rad, AGB_cumulative= biomass (pm,kco,o,tbase,Knp,IPAR,Ya)
