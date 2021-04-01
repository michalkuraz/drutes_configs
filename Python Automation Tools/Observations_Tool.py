## ---------------------------------------------------------------------------------------------------------------
## Script name: Observation tools
## Purpose of script: Construct observations record
## Author: Arij Chmeis
## Date Created: 2021-03-29
## Copyright (c) Arij Chmeis, 2021
## Email: achmeis1949@gmail.com
## ----------------------------------------------------------------------------------------------------------------
##


import json
import urllib.request
import csv
from datetime import datetime, timedelta


#-------------- Daily access ----------------#


today = datetime.today().strftime('%Y-%m-%d')
yesterday = (datetime.today() - timedelta(days=1)).strftime('%Y-%m-%d')

url = "http://foodie.lesprojekt.cz:8080/senslogOT2/rest/observation?unit_id=4151559&from_time=" + yesterday + "&to_time=" + today



#-------------- User's choice access ----------------#
##str_date = input('Enter start date (YYYY-MM-DD):\n')
##end = input('Enter end date (YYYY-MM-DD):\n')
##
##url = "http://foodie.lesprojekt.cz:8080/senslogOT2/rest/observation?unit_id=4151559&from_time=" + str_date + "&to_time=" + end



direct =  = '~ownCloud/Evaporation Model/San Juan/'

response = urllib.request.urlopen(url)
data = json.load(response)  


with open(direct + 'observations.json', 'w') as f1:
    json.dump(data, f1)
    
with open(direct + 'observations.json') as input_file:
    data1 = json.load(input_file)
    

output_file = open(direct + 'observations.out', 'w')

csv_writer = csv.writer(output_file)


#--------------------- Parsing through JSON file ------------------#

count = 0
col_names = ['Step','Psi-80cm','h-80cm','Temp-80cm','Psi-35cm','h-35cm','Temp-35cm', 'Time']
counter = 0

for line in data1:
    if count == 0:
        csv_writer.writerow(col_names)
    count+=1

    sensors = line['sensors']  
    rowlines = []
    datetimeFormat = '%Y-%m-%d %H:%M:%S+%f'


   
    for i in range(len(sensors)):
        for key in sensors[i].items():
           
            if sensors[i]['sensor_id'] == 710010002:
                value72 = (sensors[i]['value'])
                h2 = (value72 * 1000)/(999.8395 * 9.80665)
                
            elif sensors[i]['sensor_id'] == 340260002:
                value32 = (sensors[i]['value'])

            elif sensors[i]['sensor_id'] == 710010001:
                value71 = (sensors[i]['value'])
                h1 = (value71 * 1000)/(999.8395 * 9.80665)
               
            elif sensors[i]['sensor_id'] == 340260001:
                value31 = (sensors[i]['value'])
               
            if all(var in globals() for var in ('value72', 'h2', 'value32', 'value71','h1','value31')):
                rowlines.append(value72)
                rowlines.append(h2)
                rowlines.append(value32)
                rowlines.append(value71)
                rowlines.append(h1)
                rowlines.append(value31)
        
            if(len(rowlines) == 6):
                
                if counter == 0:
                    start = line['time_stamp']
                    
                counter+=1
                
    
                diff = datetime.strptime(line['time_stamp'], datetimeFormat) - datetime.strptime(start, datetimeFormat)
                if diff.days == 0:
                    time_step = diff.seconds
                else:
                    time_step = diff.seconds + (86400 * diff.days)
                
                rowlines.insert(0, time_step)
                rowlines.insert(7, line['time_stamp'])
                
                csv_writer.writerow(rowlines)
                

output_file.close()
f1.close()
input_file.close()

