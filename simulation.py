#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:51:12 2021

@author: nettiewallace & Arlo Reichert


this script will contain two simulations
the first simulation will be to set the environmental
conditions. This includes water level, quality, and soil
moisture. The second simulation will be the system's
automated reaction. This includes pulling 'sensor' data
and triggering reactions to data inputs.
"""


import matplotlib.pyplot as plt
import numpy as np
import yagmail
# start a connection
yag = yagmail.SMTP("strategicsystemspi@gmail.com")




#_________________________________________
# Defining classes

class cell1():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea, Volume):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.TN = NiLevel
        self.TP = PhLevel
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea
        self.CV = Volume

class cell2():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea, Volume):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.TN = NiLevel #TN = total nitrogen 
        self.TP = PhLevel #TP = total phosphorus 
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea
        self.CV = Volume
    
class cell3():
    def __intit__ (self, WaterLevel, SoilMoisture, NiLevel, PhLevel, FillVolume, OptimalDepth, SurfaceArea, Volume):
        self.WL = WaterLevel
        self.SM = SoilMoisture
        self.TN = NiLevel
        self.TP = PhLevel
        self.FV = FillVolume
        self.OD = OptimalDepth
        self.SA = SurfaceArea
        self.CV = Volume

class Lagoon():
    def __init__ (self, LagoonLevel, SurfaceArea, NiLevel, PhLevel, MaxLevel, MinLevel):
        self.LL = LagoonLevel
        self.SA = SurfaceArea
        self.TN = NiLevel
        self.TP = PhLevel
        self.Max = MaxLevel
        self.Min = MinLevel

        
        
#_______________________________________________
#Environmental Simulation Setup

def environmentalSimulation(MainData, WeatherData):
    #For simplicity Days starts on April 30th until October 1st, Growing season can be adjusted
    # using CSV data to assign precipitation and evapotranspiration values for each day in the simulation period
    
    days=[]
    for num in range(156):
        days=np.append(days,num) 
    days=days.astype(int)
    
    #Assigning Columns of MainData CSV to Arrays
    bio_cumlative= MainData[119:275,1]
    AGB_Daily= MainData[119:275,2]
    
    #Assigning Columns of WeatherData CSV to Arrays 
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
    Lagoon.Min =0.3 #m  # WE could potentially use this minimum value as a warning indicator during a period of insufficient water
    
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
    
    cell1.CV = cell1.SA*cell1.WL
    cell2.CV = cell2.SA*cell2.WL
    cell3.CV = cell3.SA*cell3.WL

    return (ETo, precip, days, bio_cumlative, AGB_Daily)

    
    
#_________________________________________________
#Control System Simulation

def controlSystem(): 
    # simulates the repsponce of the control system for obe single moment in time 
    # (can be one day or one reading)
    #Changes water levels based on if they need water 
    # In the case of this simulation, the culvert control is assumed to be functioning properly, so the culvert control 
    # is not explicitly defined in this code. 
    
    to = "recipient_email_address" # where to send email alerts
    subject = "ALERT from your automated retention cells!"
    
    variance = 0.02 #m for all measurements
    # First, check to see if water of overfilled anywhere:
    if cell1.WL > (cell1.OD+variance):
        cell1.FV = (cell1.OD-cell1.WL)*cell1.SA
        cell2.WL =(((cell2.WL*cell2.SA)-(cell1.FV))/cell2.SA)
        cell1.WL = cell1.WL+(cell1.FV/cell1.SA) 
        body = "Retention cell 1 has exceeded its allowable water depth. Drainage of {} to cell 2 has begun.".format((-cell1.FV))
        yag.send(to, subject, body)  # sending the email

    if cell2.WL > (cell2.OD+variance):
        cell2.FV = (cell2.OD-cell2.WL)*cell2.SA
        cell3.WL =(((cell3.WL*cell3.SA)-(cell2.FV))/cell3.SA)
        cell2.WL = cell2.WL+(cell2.FV/cell2.SA)         
        body = "Retention cell 2 has exceeded its allowable water depth. Drainage of {} to cell 3 has begun.".format((-cell2.FV))
        yag.send(to, subject, body)  # sending the email

    if cell3.WL > (cell3.OD+variance):
        cell3.FV = (cell3.OD-cell3.WL)*cell3.SA # the excess water from cell3 gets ejected from the entire system
        cell3.WL = cell3.WL+(cell3.FV/cell3.SA) 
        body = "Retention cell 3 has exceeded its allowable water depth. Drainage of {} out of cell 3 has begun.".format((-cell3.FV))
        yag.send(to, subject, body)  # sending the email

    # Second, check to see if water is underfilled anywhere:
    if cell3.WL >=(cell3.OD-variance):
        cell3.FV =0
    else:                                                                     
        cell3.FV =(cell3.OD-cell3.WL)*cell3.SA
        cell2.WL =(((cell2.WL*cell2.SA)-(cell3.FV))/cell2.SA)
        cell3.WL = cell3.WL+(cell3.FV/cell3.SA) 
        body = "Retention cell 3 has dropped below its allowable water depth. Drainage of {} to cell 3 from cell 2 has begun.".format((-cell3.FV))
        yag.send(to, subject, body)  # sending the email

    if cell2.WL >=(cell2.OD-variance):
        cell2.FV =0
    else:                                                                    
        cell2.FV =(cell2.OD-cell2.WL)*cell2.SA
        cell1.WL =((cell1.WL*cell1.SA)-((cell2.FV)))/cell1.SA
        cell2.WL = cell2.WL +(cell2.FV/cell2.SA)  
        body = "Retention cell 2 has dropped below its allowable water depth. Drainage of {} to cell 2 from cell 1 has begun.".format((-cell2.FV))
        yag.send(to, subject, body)  # sending the email
      
    if cell1.WL >=(cell1.OD-variance):
        cell1.FV =0
    else:                                                            
        cell1.FV =(cell1.OD-cell1.WL)*cell1.SA
        Lagoon.LL =((Lagoon.LL*Lagoon.SA)-cell1.FV)/Lagoon.SA
        cell1.WL = cell1.WL +(cell1.FV/cell1.SA)
        body = "Retention cell 1 has dropped below its allowable water depth. Drainage of {} to cell 1 from the lagoon source has begun.".format((-cell1.FV))
        yag.send(to, subject, body)  # sending the email


#__________________________________________
# Putting Everything Together:

def Main(MainData, WeatherData):
    
    ETo, precip, days, bio_cumlative, AGB_Daily = environmentalSimulation(MainData, WeatherData) # take the ETO and precip values for each day of the simulation
    
    TN_extract = 0.0124 #Total Nitrogen accumulation rate (percent of above ground biomass)
    TP_extract = 0.0026 #Total Phosphorus accumulation rate (percent of above ground biomass)

    Lagoon.TP = 0.55 #Should be changed to reflect actual values
    Lagoon.TN = 7.07 #Lagoon Total Nitrogen concentration (g/m3) 

    WaterLevel1 = np.array(cell1.WL) # initiating the first value for each cell (day 1). FOR USE IN GRAPHS
    WaterLevel2 = np.array(cell2.WL)
    WaterLevel3 = np.array(cell3.WL)
    
    CellVolume1 = np.array(cell1.CV) # initiating the first value for each cell (day 1). FOR USE IN GRAPHS
    CellVolume2 = np.array(cell2.CV)
    CellVolume3 = np.array(cell3.CV)
    
    FillVolume1 = np.array(cell1.FV) # initiating the first value for each cell (day 1). FOR USE IN GRAPHS
    FillVolume2 = np.array(cell2.FV)
    FillVolume3 = np.array(cell3.FV)
    
    cell1_TN = np.array(Lagoon.TN)   # initiating the first value for each cell (day 1). FOR USE IN GRAPHS
    cell2_TN = np.array(Lagoon.TN)
    cell3_TN = np.array(Lagoon.TN)

    cell1_TP = np.array(Lagoon.TP)   # initiating the first value for each cell (day 1). FOR USE IN GRAPHS
    cell2_TP = np.array(Lagoon.TP)
    cell3_TP = np.array(Lagoon.TP)

    cell1.TN = Lagoon.TN
    cell2.TN = Lagoon.TN
    cell3.TN = Lagoon.TN

    cell1.TP = Lagoon.TP
    cell2.TP = Lagoon.TP
    cell3.TP = Lagoon.TP
    
    #TN removed is amount accumulated in biomass
    #AGB_Daily is g/m2 multiplied by SA equals grams 
    cell1_TN_removed = [0]*len(days) #[grams] || initiating lists for the amount of Ni removed each day
    cell2_TN_removed = [0]*len(days) 
    cell3_TN_removed = [0]*len(days) 

    cell1_TP_removed = [0]*len(days) #[grams] || initiating lists for the amount of Ph removed each day
    cell2_TP_removed = [0]*len(days) 
    cell3_TP_removed = [0]*len(days) 

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
        CellVolume1 = np.append(CellVolume1,cell1.CV) 
        CellVolume2 = np.append(CellVolume2,cell2.CV)
        CellVolume3 = np.append(CellVolume3,cell3.CV)
        
        controlSystem() # calling the control system allows the simulation to update key variables (Primarily Water Levels) based on the daily input states
        
        cell1_TN_removed[day]=(AGB_Daily[day]*cell1.SA*TN_extract) #g
        cell2_TN_removed[day]=(AGB_Daily[day]*cell2.SA*TN_extract)
        cell3_TN_removed[day]=(AGB_Daily[day]*cell3.SA*TN_extract) 

        cell1_TP_removed[day]=(AGB_Daily[day]*cell1.SA*TP_extract) #g
        cell2_TP_removed[day]=(AGB_Daily[day]*cell2.SA*TP_extract)
        cell3_TP_removed[day]=(AGB_Daily[day]*cell3.SA*TP_extract) 
        
        cell1.TN=((cell1.TN*cell1.CV)+(cell1.FV*Lagoon.TN)-(cell1_TN_removed[day]))/(cell1.CV)
        cell2.TN=((cell2.TN*cell2.CV)+(cell2.FV*cell1.TN)-(cell2_TN_removed[day]))/(cell2.CV)
        cell3.TN=((cell3.TN*cell3.CV)+(cell3.FV*cell2.TN)-(cell3_TN_removed[day]))/(cell3.CV)
        
        cell1.TP=((cell1.TP*cell1.CV)+(cell1.FV*Lagoon.TP)-(cell1_TP_removed[day]))/(cell1.CV)
        cell2.TP=((cell2.TP*cell2.CV)+(cell2.FV*cell1.TP)-(cell2_TP_removed[day]))/(cell2.CV)
        cell3.TP=((cell3.TP*cell3.CV)+(cell3.FV*cell2.TP)-(cell3_TP_removed[day]))/(cell3.CV)

        cell1_TN = np.append(cell1_TN,cell1.TN)
        cell2_TN = np.append(cell2_TN,cell2.TN)
        cell3_TN = np.append(cell3_TN,cell3.TN)

        cell1_TP = np.append(cell1_TP,cell1.TP)
        cell2_TP = np.append(cell2_TP,cell2.TP)
        cell3_TP = np.append(cell3_TP,cell3.TP)   
        
        
    Gates_Operated=[]
    for num in range(len(days)):
        if FillVolume1[num]  or FillVolume2[num]  or FillVolume3[num] > 0:
            Gates_Operated=np.append(Gates_Operated,days[num])    
    
            
    def plotWaterLevel(cell):
        #Plot of Water Levels     
        maxWL3 = [cell3.OD + 0.02]*len(days)
        minWL3 = [cell3.OD - 0.02]*len(days)
        
        if cell == '3':
            plt.plot(days, WaterLevel3, linestyle='solid', color='black')
            plt.plot(days, minWL3, linestyle='dashed', color='red')
            plt.plot(days, maxWL3, linestyle='dashed', color='red')
            plt.ylabel('Cell 3: Water Level')
            plt.xlabel('Day')
            plt.show()
        if cell == '2': 
            plt.plot(days, WaterLevel2,linestyle='solid', color='blue')
            plt.plot(days, minWL3, linestyle='dashed', color='red')
            plt.plot(days, maxWL3, linestyle='dashed', color='red')
            plt.ylabel('Cell 2: Water Level')
            plt.xlabel('Day')
            plt.show()
        if cell == '1':
            plt.plot(days, WaterLevel1,linestyle='solid', color='green')
            plt.plot(days, minWL3, linestyle='dashed', color='red')
            plt.plot(days, maxWL3, linestyle='dashed', color='red')
            plt.ylabel('Cell 1: Water Level')
            plt.xlabel('Day')
            plt.show()

    def plotNi(cell):
        #Plot of Nitrogen Levels
        zero= [0]*len(days)
        if cell == '1':
            plt.plot(days, cell1_TN, linestyle='solid', color='green')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 1: Total Nitrogen Per Day')
            plt.xlabel('Day')
            plt.show()
        if cell == '2':
            plt.plot(days, cell2_TN, linestyle='solid', color='blue')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 2: Total Nitrogen Per Day')
            plt.xlabel('Day')
            plt.show()    
        if cell == '3':
            plt.plot(days, cell3_TN, linestyle='solid', color='black')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 3: Total Nitrogen Per Day')
            plt.xlabel('Day')
            plt.show()
            
    def plotPh(cell):
        #Plot of Phosphorus Levels
        zero= [0]*len(days)
        if cell == '1':
            plt.plot(days, cell1_TP, linestyle='solid', color='green')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 1: Total Phosphorus Per Day')
            plt.xlabel('Day')
            plt.show()            
        if cell == '2':
            plt.plot(days, cell3_TN, linestyle='solid', color='blue')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 2: Total Phosphorus Per Day')
            plt.xlabel('Day')
            plt.show()
        if cell == '3':
            plt.plot(days, cell3_TN, linestyle='solid', color='black')
            plt.plot(days, zero, linestyle='solid', color='red')
            plt.ylabel('Cell 3: Total Phosphorus Per Day')
            plt.xlabel('Day')
            plt.show()
            
    def plotFillVolume(cell):
        #Plot of Fill Volumes for each day
        if cell == '1':
            plt.scatter(days, FillVolume1, color='blue')
            plt.ylabel('Cell 1: Daily Fill Volume')
            plt.xlabel('Day')
            plt.show()
        if cell == '2':
            plt.scatter(days, FillVolume2, color='blue')
            plt.ylabel('Cell 2: Daily Fill Volume')
            plt.xlabel('Day')
            plt.show()
        if cell == '3':
            plt.scatter(days, FillVolume3, color='blue')
            plt.ylabel('Cell 3: Daily Fill Volume')
            plt.xlabel('Day')
            plt.show()        
            
    def plotCellVolume(cell):
        #Plot of Cell Volumes for each day
        #We might want to adjust the scale or the type of graph for better presentation
        if cell == '1':
            plt.plot(days, CellVolume1,linestyle='solid', color='blue')
            plt.ylabel('Cell 1: Daily Cell Volume')
            plt.xlabel('Day')
            plt.show()  
        if cell == '2':
            plt.plot(days, CellVolume2,linestyle='solid', color='blue')
            plt.ylabel('Cell 2: Daily Cell Volume')
            plt.xlabel('Day')
            plt.show()  
        if cell == '3':
            plt.plot(days, CellVolume3,linestyle='solid', color='blue')
            plt.ylabel('Cell 3: Daily Cell Volume')
            plt.xlabel('Day')
            plt.show()       
    
    cell = '1' # change this variable depending on which cell graphs you want to see
    plotWaterLevel(cell)
#    plotNi(cell)
#    plotPh(cell)
#    plotFillVolume(cell)
#    plotCellVolume(cell)
    
    
    
#Read main data CSV
MainData = np.loadtxt(fname = 'Output_data.csv', delimiter=',')
#Read Weather data CSV 
WeatherData = np.loadtxt(fname = 'Weather_data.csv', delimiter=',')

    
Main(MainData, WeatherData) 

''' END OF PROGRAM '''
