#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:51:12 2021

@author: nettiewallace & Arlo Reichert & Owen Thoresby


this script will contain two simulations
the first simulation will be to set the environmental
conditions. This includes water level, quality, and soil
moisture. The second simulation will be the system's
automated reaction. This includes pulling 'sensor' data
and triggering reactions to data inputs.
"""


import matplotlib.pyplot as plt
import numpy as np
#import smtplib, ssl
import math
import pandas as pd

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

def environmentalSimulation(data_file):
    #For simplicity Days starts on April 30th until October 1st, Growing season can be adjusted
    # using CSV data to assign precipitation and evapotranspiration values for each day in the simulation period

    
    #Biomass Variables
    pm = 25 #Maximum net daily photosynthesis rate at 20°C (gCO2/m2day)
    kco=0.65 #conversion of carbon dioxide to AGB (g/g)
    o=1.09 #Temperature constant 
    tbase=5 #Base air temperature (C)
    Knp=1  #Nutrient availability , maximum of 1 when no nutrient limitation 
    IPAR=5.8 # Daily total incident photosynthetically active radiation (moles/m2day)
    Ya=0.0025 #Leaf and shoot mortality constant 

    #Data file containing [day;percip;ETo;Temp;Radiation Data]
    #WeatherData = np.loadtxt(fname = 'weather_data.csv', delimiter=',')
    
    #Day, Temp, PAR	, Percip (mm), Average High Temp , Average Low Temp)

    day= data_file[119:275,0]
    precip= data_file[119:275,3] #mm
    precip=precip/1000
    Temp= data_file[119:275,1] 
    RAD=data_file[119:275,2]
    
    # This section predicts the quantity of biomass produced 
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
        
    #LAI
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

#    AGB_cumulative = []
#    totalbiomass = []
    mortality_cumulative = np.cumsum(mortality); #cumulative mortality 
    bio_cumlative = np.cumsum(AGB_Daily); #Biomass production without mortality (g/m2)
    bio_cumlative = bio_cumlative-mortality_cumulative; #Total cumulative biomass (g/m2)

    #Calibration Data ______________________________

#    calibration_data= np.loadtxt(fname = 'gleam_data.csv', delimiter=',')
#    Ep_wpg=calibration_data

    #HARGREAVES DATA ________________________________________

    #Radiaiton Variables
    lat = 0.85172 #radians 
    #Extraterrestial Radiation Variables
    
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
      solar_declination[num]=(0.409*(np.sin((0.0172*(num+119)-1.39))));
      sunset_hour[num] = np.arccos(-np.tan(lat)*np.tan(solar_declination[num]));
      realitive_distance[num] = 1+0.033*np.cos(0.0172*(num+119));
      RA[num]=0.408*37.6*realitive_distance[num]*((sunset_hour[num]*np.sin(lat)*np.sin(solar_declination[num]))+np.cos(lat)*np.cos(solar_declination[num])*np.sin(sunset_hour[num]))
      #RA in mm/day

    #ETo Hargraves Method
    for num in range(len(day)):
        ETo[num]=(high_temp[num]-low_temp[num])        
        ETo[num]=(KET*RA[num]*(mean_temp[num]+17.8)*np.sqrt((high_temp[num]-low_temp[num])))
        ETo[num]=ETo[num]/1000
    

    #Calibrating Hargraves Method Using Room Mean Square Deviation (RMSD)
    #KET variable was manually adjusted based on below calcuations 
    # NUM =(Ep_wpg-ETo)**2;
    # RMSD = (sum[NUM]/(len[ETo]))**0.5 
    
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

    return (ETo, precip, day, bio_cumlative, AGB_Daily)


#_____________________________________________________
# Time to fill function

def timeToFill():
    
    variance = 0.02 #m
    
    # Culvert Dimensions, Slop and Friction Coefficient
    Diam = 0.6096
    R = Diam / 2
    SLOPE = 0.01
    n = 0.022
    
    # Initial Depths at start of transfer 
    waterLevel1 = cell1.WL
    waterLevel2 = cell2.WL
    waterLevel3 = cell3.WL
    lagoonLevel = Lagoon.LL
    
    # Empty lists that store the depths for each cell
    DepthL = []
    Depth1 = []
    Depth2 = []
    Depth3 = []
    FlowL = []
    Flow1 = []
    Flow2 = []
    Flow3 = []
    
    # Transfer between cells 2 and 3
    # Continues loop until Cell 3 Depth reaches 0.5 m
    # Alternates between 2 open flow equations based on Cell 2 Depth
    
    if waterLevel3 <= (cell3.OD - variance):          # Case where cells are low on water
    
        while waterLevel3 < cell3.OD:
            while waterLevel2 > 0.1:
                if waterLevel2 > R:
                    H1 = (2 * R - waterLevel2)
                    Theta = 2 * math.acos((R - H1) / R)
                    A = math.pi * (R ** 2) - (((R ** 2) * (Theta - math.sin(Theta))) / 2)
                    P = (2 * math.pi * R) - R * Theta
                    Rh = A / P
                    Q = (A * (Rh ** (2 / 3)) * (SLOPE ** 0.5)) / n
                elif waterLevel2 < R:
                    H1 = waterLevel2
                    Theta = 2 * math.acos((R - H1) / R)
                    A = (((R ** 2) * (Theta - math.sin(Theta))) / 2)
                    P = R * Theta
                    Rh = A / P
                    Q = (A * (Rh ** (2 / 3)) * (SLOPE ** 0.5)) / n
                Depth2.append(waterLevel2)
                Flow2.append(Q)
                waterLevel2 += (Q * (-60)) / cell2.SA
                Depth3.append(waterLevel3)
                waterLevel3 += (Q * 60) / cell3.SA
                if waterLevel3 >= cell3.OD:      #potentially not necessary to specify within WHILE loop
                    break
        
        # Transfer between cells 1 and 2
        # Continues loop until Cell 2 Depth reaches 0.5 m
        # Alternates between 2 open flow equations based on Cell 2 Depth
        while waterLevel2 < cell2.OD:
            while waterLevel1 > 0.1:
                if waterLevel1 > R:
                    H1 = (2 * R - waterLevel1)
                    Theta = 2 * math.acos((R - H1) / R)
                    A = math.pi * (R ** 2) - (((R ** 2) * (Theta - math.sin(Theta))) / 2)
                    P = (2 * math.pi * R) - R * Theta
                    Rh = A / P
                    Q = (A * (Rh ** (2 / 3)) * (SLOPE ** 0.5)) / n
                elif waterLevel1 < R:
                    H1 = waterLevel1
                    Theta = 2 * math.acos((R - H1) / R)
                    A = (((R ** 2) * (Theta - math.sin(Theta))) / 2)
                    P = R * Theta
                    Rh = A / P
                    Q = (A * (Rh ** (2 / 3)) * (SLOPE ** 0.5)) / n
                Depth1.append(waterLevel1)
                Flow1.append(Q)
                waterLevel1 += (Q * (-60)) / cell1.SA
                Depth2.append(waterLevel2)
                waterLevel2 += (Q * 60) / cell1.SA
                if waterLevel2 >= cell2.OD:     #potentially not necessary to specify within WHILE loop
                    break
        
        
        # Transfer between cells 1 and lagoon
        # Continues loop until Cell 1 Depth reaches 0.5 m
        # Assumed an average velocity for when the water is higher than the culvert diameter
        # Changes to open channel flow if Lagoon height drops below Diam
        
        while waterLevel1 < cell1.OD:
            if lagoonLevel > Diam:
                QL = math.sqrt(lagoonLevel * 2 * 9.81) * (math.pi * (R**2))
            elif lagoonLevel < Diam:
                H1 = (2 * R - lagoonLevel)
                Theta = 2 * math.acos((R - H1) / R)
                A = math.pi * (R ** 2) - (((R ** 2) * (Theta - math.sin(Theta))) / 2)
                P = (2 * math.pi * R) - R * Theta
                Rh = A / P
                QL = (A * (Rh ** (2 / 3)) * (SLOPE ** 0.5)) / n
            DepthL.append(lagoonLevel)
            FlowL.append(QL)
            lagoonLevel += (Q * (-60)) / Lagoon.SA
            Depth1.append(waterLevel1)
            waterLevel1 += (Q * 60) / cell1.SA
            if waterLevel1 >= cell1.OD:     #potentially not necessary to specify within WHILE loop
                break
        
        TotalTime = (len(Depth3) + len(Depth2) + len(Depth1)) / 60
    
    elif waterLevel1 >= (cell1.OD + variance):            # Case where cells have too much water
        while waterLevel1 > cell1.OD:
            HL1 = (waterLevel1 + 0.2) - waterLevel2
            H1 = (2 * R - waterLevel1)
            Theta1 = 2 * math.acos((R - H1) / R)
            A1 = math.pi * (R ** 2) - (((R ** 2) * (Theta1 - math.sin(Theta1))) / 2)  # Cell1 Wetted Area
            P1 = (2 * math.pi * R) - R * Theta1  # Cell 1 Wetted Perimeter
            Rh1 = A1 / P1  # Cell 1 Hydraulic Radius
            Q1 = (A1 * (Rh1 ** (2 / 3)) * (SLOPE ** 0.5)) / n  # Cell1 Flow Rate
            Depth1.append(waterLevel1)
            waterLevel1 += (Q1 * (-60)) / cell1.SA  # Water Leaving Cell3
            waterLevel2 += (Q1 * 60) / cell2.SA
            
            while waterLevel2 > cell2.OD:
                H2 = (2 * R - waterLevel2)
                Theta2 = 2 * math.acos((R - H2) / R)
                A2 = math.pi * (R ** 2) - (((R ** 2) * (Theta2 - math.sin(Theta2))) / 2)  # Cell1 Wetted Area
                P2 = (2 * math.pi * R) - R * Theta2 # Cell 1 Wetted Perimeter
                Rh2 = A2 / P2  # Cell 1 Hydraulic Radius
                Q2 = (A2 * (Rh1 ** (2 / 3)) * (SLOPE ** 0.5)) / n  # Cell1 Flow Rate
                Depth2.append(waterLevel2)
                waterLevel2 += (Q2 * (-60)) / cell2.SA  # Water Leaving Cell3
                waterLevel3 += (Q2 * 60) / cell3.SA
                
                while waterLevel3 > cell3.OD:
                    if waterLevel3 > R:  # Cell3.WL > R (Inlet Control)
                        Cd = 1.0  # Coefficient of Discharge (Round)
                        A3 = math.pi * (R ** 2)  # Cell3 Wetted Area
                        h = waterLevel3 - R
                        Q3 = Cd * A3 * math.sqrt(9.81 * h)  # Cell3 Flow Rate
                        if waterLevel3 < R:  # Cell3.WL < R (Open Channel Flow)
                            H3 = (2 * R - waterLevel3)
                            Theta3 = 2 * math.acos((R - H3) / R)
                            A3 = math.pi * (R ** 2) - (((R ** 2) * (Theta3 - math.sin(Theta3))) / 2)  # Cell3 Wetted Area
                            P3 = (2 * math.pi * R) - R * Theta3  # Cell3 Wetted Perimeter
                            Rh3 = A3 / P3  # Cell3 Hydraulic Radius
                            Q3 = (A3 * (Rh3 ** (2 / 3)) * (SLOPE ** 0.5)) / n  # Cell3 Flow Rate
                        Depth3.append(waterLevel3)
                        waterLevel3 += (Q3 * (-60)) / cell3.SA
            
        TotalTime = (len(Depth2) + len(Depth1) + len(Depth3)) / 60
        
    else:
        return

    return TotalTime
 
'''
#_____________________________________________________
# Email Alert function 

def emailAlert(Date, Cell, Problem, Solution, Volume):
    
    port = 465  # For SSL
    smtp_server = "smtp.gmail.com"
    sender_email = "strategicsystemspi@gmail.com"  # Enter your address
    receiver_email = "strategicsystemspi@gmail.com"  # Enter receiver address
    password = input("Type your password and press enter: ")
    
    for num in range(1):
        message = """\
        Subject: Update From Automated Retention Cell.
        

        This message is sent from Python.
    
        Update from {date}: Retention cell {cellnum} is {problem} its allowable water depth. {solution} of {volume:.2f} m^3 has begun."""
        
    
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
            server.login(sender_email, password)
            server.sendmail(sender_email, receiver_email, message.format(date = Date, cellnum =Cell, problem = Problem, solution = Solution , volume = Volume)) 
'''
 
#_________________________________________________
#Control System Simulation

def controlSystem(date): 
    # simulates the repsponce of the control system for obe single moment in time 
    # (can be one day or one reading)
    #Changes water levels based on if they need water 
    # In the case of this simulation, the culvert control is assumed to be functioning properly, so the culvert control 
    # is not explicitly defined in this code. 

    date_range = pd.date_range(start= '2018,04,30', end='2018,10,01')
    date_range = date_range.date


    variance = 0.02 #m for all measurements
    # First, check to see if water of overfilled anywhere:
    if cell1.WL > (cell1.OD+variance):
        cell1.FV = (cell1.OD-cell1.WL)*cell1.SA
        cell2.WL =(((cell2.WL*cell2.SA)-(cell1.FV))/cell2.SA)
        cell1.WL = cell1.WL +(cell1.FV/cell1.SA)
        '''Day = str(date_range(date))
        Cell = 1
        Problem = 'Above'
        Solution = 'Discharge'
        Volume = cell1.FV
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        
    if cell2.WL > (cell2.OD+variance):
        cell2.FV = (cell2.OD-cell2.WL)*cell2.SA
        cell3.WL =(((cell3.WL*cell3.SA)-(cell2.FV))/cell3.SA)
        cell2.WL = cell2.WL +(cell2.FV/cell2.SA) 
        '''Day = str(date_range(date))
        Cell = 2
        Problem = 'Above'
        Solution = 'Discharge'
        Volume = cell2.FV
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        
    if cell3.WL > (cell3.OD+variance):
        cell3.FV = (cell3.OD-cell3.WL)*cell3.SA # the excess water from cell3 gets ejected from the entire system
        cell3.WL = cell3.WL+(cell3.FV/cell3.SA) 
        '''Day = str(date_range(date))
        Cell = 3
        Problem = 'Above'
        Solution = 'Discharge'
        Volume = cell3.FV
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        
    # Second, check to see if water is underfilled anywhere:
    if cell3.WL >=(cell3.OD-variance):
        cell3.FV =0
    else:                                                                     
        cell3.FV =(cell3.OD-cell3.WL)*cell3.SA
        cell2.WL =(((cell2.WL*cell2.SA)-(cell3.FV))/cell2.SA)
        '''Day = str(date_range(date))
        Cell = 3
        Problem = 'below'
        Solution = 'Addition'
        Volume = abs(cell3.FV)
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        cell3.WL = cell3.WL+(cell3.FV/cell3.SA) 
        
    if cell2.WL >=(cell2.OD-variance):
        cell2.FV =0
    else:                                                                    
        cell2.FV =(cell2.OD-cell2.WL)*cell2.SA
        cell1.WL =((cell1.WL*cell1.SA)-((cell2.FV)))/cell1.SA
        '''Day = str(date_range(date))
        Cell = 2
        Problem = 'below'
        Solution = 'Addition'
        Volume = abs(cell2.FV)
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        cell2.WL = cell2.WL +(cell2.FV/cell2.SA)  
      
    if cell1.WL >=(cell1.OD-variance):
        cell1.FV =0
    else:                                                            
        cell1.FV =(cell1.OD-cell1.WL)*cell1.SA
        Lagoon.LL =((Lagoon.LL*Lagoon.SA)-cell1.FV)/Lagoon.SA
        '''Day = str(date_range(date))
        Cell = 1
        Problem = 'below'
        Solution = 'Addition'
        Volume = cell1.FV
        emailAlert(Day, Cell, Problem, Solution, Volume)'''
        cell1.WL = cell1.WL +(cell1.FV/cell1.SA)


#__________________________________________
# Putting Everything Together:

def Main(data_file):
    
    ETo, precip, days, bio_cumlative, AGB_Daily = environmentalSimulation(data_file) # take the ETO and precip values for each day of the simulation
    
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
    
    fillTimes = np.array([])

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
        
        fillTimes = np.append(fillTimes,timeToFill())
        controlSystem(day)              # calling the control system allows the simulation to update key variables (Primarily Water Levels) based on the daily input states
        
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
    
            
            
    def plotFillTimes(fillTimes):
        #Plot of fill times for each day when intervention is necessary
        fakedays = np.array([])
        datalength = len(fillTimes)
        for i in range(0,datalength):
            fakedays = np.append(fakedays, i)
        plt.scatter(fakedays, fillTimes, linestyle='solid', color='green')
        plt.ylabel('Duration of Adjustment (hours)')
        plt.xlabel('Day')
        plt.show()
            
        
        
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
    
#    cell = '1' # change this variable depending on which cell graphs you want to see
#    plotWaterLevel(cell)
#    plotNi(cell)
#    plotPh(cell)
#    plotFillVolume(cell)
#    plotCellVolume(cell)
#    plotFillTimes(fillTimes)
    
    
#Read Weather data CSV 
data_file = np.loadtxt(fname ='Sample data file2.csv', delimiter=',')

ETo, precip, days, bio_cumlative, AGB_Daily = environmentalSimulation(data_file)   


Main( data_file) 

''' END OF PROGRAM '''






    










    








