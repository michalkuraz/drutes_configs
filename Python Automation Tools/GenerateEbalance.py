## ---------------------------------------------------------------------------------------------------------------
## Script name: Ebalance
## Purpose of script: Generate Ebalance.in file
## Author: Arij Chmeis
## Date Created: 2021-03-29
## Copyright (c) Arij Chmeis, 2021
## Email: achmeis1949@gmail.com
## ----------------------------------------------------------------------------------------------------------------
##


import csv
import math
from datetime import datetime, timedelta
from calendar import monthrange


direct = '~ownCloud/Thesis/Case study'
out_direct = '~ownCloud/Thesis/Case study/drutes.conf/evaporation/'


def is_leap(year):
    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
                result = True
            else:
                result = False
        else:
            result = True
    else:
        result = False

    return result



def solar_rad(day, mon, year, lat, t):
    a = 0.25
    b = 0.5

    r_lat = - lat * math.pi / 180

    G_sc = 1360

    j = int((275 * mon/9) - 30 + day) - 2
    if mon < 3:
        j = j + 2
    elif mon > 2 and is_leap(year) == True:
        j = j + 1

    delta = 0.409 * math.sin( (2 * math.pi * j/ 365) - 1.39)
    
    ws = math.acos(-math.tan(r_lat) * math.tan(delta))

    dr = 1 + 0.033 * math.cos(2*math.pi*j / 365)
    
    R_a = G_sc * dr * (ws * math.sin(r_lat) * math.sin(delta) + math.cos(r_lat) * math.cos(delta) * math.sin(ws))/ math.pi
        
    if R_a < 0:
        R_a = 0

    sin_e =  math.sin(r_lat) * math.sin(delta) + math.cos(r_lat) * math.cos(delta) * math.cos( 2 * math.pi * (t - 12) / 24)

    R_s = max(sin_e * R_a, 0)

    S_t = (a + b) * R_s
    
    return S_t




with open(direct + 'meteo_data.out', 'r') as inputfile:
    reader = csv.DictReader(inputfile)

    ebalance_file = open(out_direct + 'ebalance.in' , 'w')
    writer = csv.writer(ebalance_file)

    col_names = ['#t [s]', 'S_t [W.m-2]', 'R_a [W.m-2]', 'T_a [ºC]',
                 'wind speed [m/s]', 'fraction of cloud cover [-]', 'relative humidity [-]']
    count = 0
    time_step = 0
    
   
    lat = 31.5375
    long = 68.5386
    

    for line in reader:
        rowlines = []
        if count == 0:
            writer.writerow(col_names)
        count+=1

        date = line['Date']
        d = int(date[0:2])
        m = int(date[2:4])
        y = int(date[4:8])
        t = int(line['Hour'])
        
        S_t = solar_rad(d, m, y, lat, t)

        rowlines.append(time_step)
        time_step = time_step + 3600
        

        
        rowlines.append(S_t)
        
        rowlines.append('1')
        rowlines.append(line['Temp [ºC]'])

        ws = line['FF [km/hr]']
        
        if ws == 'SAN':
            u = 0
        else:
            u = int(ws)
            u = (u * 5) / 18

        rowlines.append(u)
        rowlines.append('0')
        H = float(line['Humidity [%]'])/100
        rowlines.append(H)
        
        writer.writerow(rowlines)



ebalance_file.close()











    

   
