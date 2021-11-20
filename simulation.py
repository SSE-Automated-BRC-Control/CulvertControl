#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:51:12 2021

@author: nettiewallace
"""

# this script will contain two simulations
# the first simulation will be to set the environmental
# conditions. This includes water level, quality, and soil
# moisture. The second simulation will be the system's
# automated reaction. This includes pulling 'sensor' data
# and triggering reactions to data inputs.

# the entire system will be looped to simulate a reaction
# to environmental inputs:

import matplotlib.pyplot as plt
import numpy as np
my_data = np.loadtxt(fname = 'Output_data.csv', delimiter=',')

#For simplicity Days starts on April 30th until October 1st, Growing season can be adjusted
days=[]
for num in range(156):
    days=np.append(days,num) 
days=days.astype(int)

#Data From CSV For Later
bio_cumlative=my_data[119:275,1]
AGB_Daily=my_data[119:275,2]

#Reading Weather data CSV 
weather_data = np.loadtxt(fname = 'Weather_data.csv', delimiter=',')

#Assigning Columns of CSV to Arrays 
percip=weather_data[119:275,1] #mm
percip=percip/1000 #m
ETo=weather_data[119:275,2] #mm
ETo=ETo/1000 #m
Temp=weather_data[119:275,3] 

#Variables for Calcuaitons Based on Oakburn Estimations  
#Simplified to rectangles, should be updated to trapazodial prisms 
cell1_SA=9000 #m^2
cell2_SA=9000
cell3_SA=6000
Lagoon_SA=13000

#Lagoon Water Level Range 
lagoon_max_level=1.5 #m
lagoon_min_level=0.3 #m

#Lagoon Level Assumed to be max in Spring, Adjusted for ETo and Percip
LagoonLevel=[lagoon_max_level]*len(days)
for day in range(1,len(days)):
    LagoonLevel[day]=LagoonLevel[day-1]-ETo[day]+percip[day]

#Optium Depth and Variance of The Cells
opt_depth=0.5 #m
variance = 0.02 #m

#Assume Water Levels are to Optium on fill date (April 31)
WaterLevel1=[opt_depth]*len(days)
WaterLevel2=[opt_depth]*len(days)
WaterLevel3=[opt_depth]*len(days)

WaterLevel1[0]=WaterLevel1[0]-ETo[0]+percip[0]
WaterLevel2[0]=WaterLevel2[0]-ETo[0]+percip[0]
WaterLevel3[0]=WaterLevel3[0]-ETo[0]+percip[0]
LagoonLevel[0]=LagoonLevel[0]-ETo[0]+percip[0]

#How much each cell needs to be filled up by
fill_vol1=[0]*len(days)
fill_vol2=[0]*len(days)
fill_vol3=[0]*len(days)

for day in range(1,len(days)):
    #Changing Orignial Water Levels Based on ETo and Percip 
    WaterLevel1[day]=WaterLevel1[day-1]-ETo[day]+percip[day]
    WaterLevel2[day]=WaterLevel2[day-1]-ETo[day]+percip[day]    
    WaterLevel3[day]=WaterLevel3[day-1]-ETo[day]+percip[day]
    
    #Changing water levels based on if they need water 
    if WaterLevel3[day]>=(opt_depth-variance):
        fill_vol3[day]=0
    elif WaterLevel3[day]<(opt_depth-variance):  
        fill_vol3[day]=(opt_depth-WaterLevel3[day])*cell3_SA
        WaterLevel2[day]=(((WaterLevel2[day]*cell2_SA)-(fill_vol3[day]))/cell2_SA)
    WaterLevel3[day]= WaterLevel3[day]+(fill_vol3[day]/cell3_SA) 
        
    if WaterLevel2[day]>=(opt_depth-variance):
        fill_vol2[day]=0
    elif WaterLevel2[day]<(opt_depth-variance):    
        fill_vol2[day]=(opt_depth-WaterLevel2[day])*cell2_SA
        WaterLevel1[day]=((WaterLevel1[day]*cell1_SA)-((fill_vol2[day])))/cell1_SA
    WaterLevel2[day]=WaterLevel2[day]+(fill_vol2[day]/cell2_SA)  
      
    if WaterLevel1[day]>=(opt_depth-variance):
        fill_vol1[day]=0
    elif WaterLevel1[day]<(opt_depth-variance):    
        fill_vol1[day]=(opt_depth-WaterLevel1[day])*cell1_SA
        LagoonLevel[day]=((LagoonLevel[day]*Lagoon_SA)-fill_vol1[day])/Lagoon_SA
    WaterLevel1[day]=WaterLevel1[day]+(fill_vol1[day]/cell1_SA)
 
Gates_Operated=[]
for num in range(len(days)):
    if fill_vol1[num]  or fill_vol2[num]  or fill_vol3[num] > 0:
        Gates_Operated=np.append(Gates_Operated,days[num])    

#Plot of Water Levels        
plt.plot(days, WaterLevel3)
#plt.plot(days, WaterLevel2,linestyle="--")
#plt.plot(days, WaterLevel1,linestyle="--")
      
#_____________________________________________________________
#Alternative Appraoch if we wanted to put the logic in a funciton     

def ReadWaterLevel(opt_depth,variance, WaterLevel1,WaterLevel2,WaterLevel3):
    #This function reads each water level and determines if it needs water 
    #This is mainly a place holder for when sesnors need to be read 
    #In a real situation sensor data would be called from the Pi
    #by this fucntion 
    
    if WaterLevel1 < opt_depth-variance :
        water_level_1 = 1
    else:
        water_level_1=0
    if WaterLevel2 < opt_depth-variance :
        water_level_2 = 1
    else:
         water_level_2=0
    if WaterLevel3 < opt_depth-variance :
        water_level_3 = 1
    else:
         water_level_3=0  
    return(water_level_1,water_level_2,water_level_1)       


       

    








