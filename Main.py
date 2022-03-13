#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

def biomass (pm,kco,o,tbase,Knp,IPAR,Ya,day):
  
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


#Calibration Data ______________________________

def ETo_calc (day):
    calibration_data= np.loadtxt(fname = 'gleam_data.csv', delimiter=',')
    Ep_wpg=calibration_data

#HARGREAVES DATA ________________________________________

#Radiaiton Variables
    lat = 0.85172 #radians 
#Extraterrestial Radiation Variables
    data_file = np.loadtxt(fname ='Sample data file2.csv', delimiter=',')
    KET=0.00236
    RA=[0]*len(day)
    ETo=[0]*len(day)
    mean_temp=[0]*len(day)
    
    high_temp = data_file[119:275,4]
    low_temp= data_file[119:275,5]
  
    for num in range(len(day)):
      mean_temp[num]=(high_temp[num]+low_temp[num])/2

#Extraterrestial Radiation Calculation
    solar_declination=[0]*len(day)
    realitive_distance=[0]*len(day)
    sunset_hour=[0]*len(day)
  
    for num in range(len(day)): 
      solar_declination[num]=(0.409*(np.sin((0.0172*num)-1.39)));
      sunset_hour[num] = np.arccos(-np.tan(lat)*np.tan(solar_declination[num]));
      realitive_distance[num] = 1+0.033*np.cos(0.0172*num);
      RA[num]=0.408*37.6*realitive_distance[num]*((sunset_hour[num]*np.sin(lat)*np.sin(solar_declination[num]))+np.cos(lat)*np.cos(solar_declination[num])*np.sin(sunset_hour[num]))
      #RA in mm/day

#ETo Hargraves Method
    for num in range(len(day)):
        ETo[num]=(high_temp[num]-low_temp[num])        
        
        ETo[num]=KET*RA[num]*(mean_temp[num]+17.8)*np.sqrt((high_temp[num]-low_temp[num]))

#Calibrating Hargraves Method Using Room Mean Square Deviation (RMSD)
#KET variable was manually adjusted based on below calcuations 
   # NUM =(Ep_wpg-ETo)**2;
   # RMSD = (sum[NUM]/(len[ETo]))**0.5
    
    return (ETo) 

ETO_Daily = ETo_calc(day)

AGB_Daily, lai, par, rad, AGB_cumulative= biomass (pm,kco,o,tbase,Knp,IPAR,Ya, day)
