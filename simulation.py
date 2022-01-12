#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:51:12 2021

@author: nettiewallace


this script will contain two simulations
the first simulation will be to set the environmental
conditions. This includes water level, quality, and soil
moisture. The second simulation will be the system's
automated reaction. This includes pulling 'sensor' data
and triggering reactions to data inputs.
"""


import matplotlib.pyplot as plt
import numpy as np




#_________________________________________
# Defining classes

class cell1():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.NiL = NiLevel
        self.PhL = PhLevel
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea

class cell2():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.NiL = NiLevel
        self.PhL = PhLevel
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea
    
class cell3():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.NiL = NiLevel
        self.PhL = PhLevel
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea

class Lagoon():
    def __init__ (self, LagoonLevel, SurfaceArea, MaxLevel, MinLevel):
        self.LL = LagoonLevel
        self.SA = SurfaceArea
        self.Max = MaxLevel
        self.Min = MinLevel

        
        
#_______________________________________________
#Environmental Simulation Setup

def environmentalSimulation(MainData, WeatherData):
    #For simplicity Days starts on April 30th until October 1st, Growing season can be adjusted
    days=[]
    for num in range(156):
        days=np.append(days,num) 
    days=days.astype(int)
    
    #Data From CSV For Later
#    bio_cumlative= MainData[119:275,1]
#    AGB_Daily= MainData[119:275,2]
    
    #Assigning Columns of CSV to Arrays 
    precip= WeatherData[119:275,1] #mm
    precip= precip/1000 #m
    ETo= WeatherData[119:275,2] #mm
    ETo=ETo/1000 #m
#    Temp= WeatherData[119:275,3] 
    
    #Variables for Calcuaitons Based on Oakburn Estimations  
    #Simplified to rectangles, should be updated to trapazodial prisms 
    cell1.SA=9000 #m^2
    cell2.SA=9000
    cell3.SA=6000
    Lagoon.SA=13000
    
    #Lagoon Water Level Range 
    Lagoon.Max =1.5 #m
    Lagoon.Min =0.3 #m
    
    #Lagoon Level Assumed to be max in Spring, Adjusted for ETo and Percip
    Lagoon.LL = Lagoon.Max
#    for day in range(1,len(days)):
#        Lagoon.LL[day]= Lagoon.LL[day-1]-ETo[day]+percip[day]
    
    #Optium Depth of The Cells
    cell1.OD =0.5 #m
    cell2.OD =0.5 #m
    cell3.OD =0.5 #m

    #Assume Water Levels are optimal on fill date (April 31)
    cell1.WL = cell1.OD-ETo[0]+precip[0]
    cell2.WL = cell2.OD-ETo[0]+precip[0]
    cell3.WL = cell3.OD-ETo[0]+precip[0]
    
    Lagoon.LL= Lagoon.LL-ETo[0]+precip[0] 
    
    #How much each cell needs to be filled up by initially
    cell1.FV = 0
    cell2.FV = 0
    cell3.FV = 0

    return (ETo, precip, days)

    
    
#_________________________________________________
#Control System Simulation

def controlSystem(): 
    # simulates the repsponce of the control system for obe single moment in time 
    # (can be one day or one reading)
    #Changes water levels based on if they need water 
    
    variance = 0.02 #m for all measurements
    
    if cell3.WL >=(cell3.OD-variance):
        cell3.FV =0
    else:                                                                     
        cell3.FV =(cell3.OD-cell3.WL)*cell3.SA
        cell2.WL =(((cell2.WL*cell2.SA)-(cell3.FV))/cell2.SA)
    cell3.WL = cell3.WL+(cell3.FV/cell3.SA) 
        
    if cell2.WL >=(cell2.OD-variance):
        cell2.FV =0
    else:                                                                    
        cell2.FV =(cell2.OD-cell2.WL)*cell2.SA
        cell1.WL =((cell1.WL*cell1.SA)-((cell2.FV)))/cell1.SA
    cell2.WL = cell2.WL +(cell2.FV/cell2.SA)  
      
    if cell1.WL >=(cell1.OD-variance):
        cell1.FV =0
    else:                                                            
        cell1.FV =(cell1.OD-cell1.WL)*cell1.SA
        Lagoon.LL =((Lagoon.LL*Lagoon.SA)-cell1.FV)/Lagoon.SA
    cell1.WL = cell1.WL +(cell1.FV/cell1.SA)

    
    

#__________________________________________
# Putting Everything Together:

def Main(MainData, WeatherData):
    
    ETo, precip, days = environmentalSimulation(MainData, WeatherData)
    
    WaterLevel1 = np.array(cell1.WL)
    WaterLevel2 = np.array(cell2.WL)
    WaterLevel3 = np.array(cell3.WL)
    
    FillVolume1 = np.array(cell1.FV)
    FillVolume2 = np.array(cell2.FV)
    FillVolume3 = np.array(cell3.FV)
    
    for day in range(1,len(days)):
        #Changing Orignial Water Levels Based on ETo and Percip 
        cell1.WL= cell1.WL-ETo[day]+precip[day]
        cell2.WL= cell2.WL-ETo[day]+precip[day]    
        cell3.WL= cell3.WL-ETo[day]+precip[day]
        
        #setting up daily value lists to graph afterwards
        WaterLevel1 = np.append(WaterLevel1,cell1.WL)
        WaterLevel2 = np.append(WaterLevel2,cell2.WL)
        WaterLevel3 = np.append(WaterLevel3,cell3.WL)
        FillVolume1 = np.append(FillVolume1,cell1.FV)
        FillVolume2 = np.append(FillVolume2,cell2.FV)
        FillVolume3 = np.append(FillVolume3,cell3.FV)
        
        controlSystem()
     
    Gates_Operated=[]
    for num in range(len(days)):
        if FillVolume1[num]  or FillVolume2[num]  or FillVolume3[num] > 0:
            Gates_Operated=np.append(Gates_Operated,days[num])    
    
    #Plot of Water Levels     
    maxWL3 = [cell3.OD]*len(days)
    minWL3 = [cell3.OD - 0.02]*len(days)
    plt.plot(days, WaterLevel3, linestyle='solid', color='black')
    plt.plot(days, minWL3, linestyle='dashed', color='red')
    plt.plot(days, maxWL3, linestyle='dashed', color='red')
    plt.ylabel('Water Level')
    plt.xlabel('Day')
    plt.show()
#    plt.plot(days, WaterLevel2,linestyle="--")
#    plt.plot(days, WaterLevel1,linestyle="--")
    

    
#Read main data CSV
MainData = np.loadtxt(fname = 'Output_data.csv', delimiter=',')
#Read Weather data CSV 
WeatherData = np.loadtxt(fname = 'Weather_data.csv', delimiter=',')

    
Main(MainData, WeatherData) 


    








