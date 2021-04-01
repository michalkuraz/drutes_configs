## ---------------------------------------------------------------------------------------------------------------
## Script name: Input tools
## Purpose of script: Automate meteo data input
## Author: Arij Chmeis
## Date Created: 2021-03-29
## Copyright (c) Arij Chmeis, 2021
## Email: achmeis1949@gmail.com
## ----------------------------------------------------------------------------------------------------------------
##


import requests,zipfile, io
import csv
from datetime import datetime, timedelta
from pytz import timezone    


#--- Finding specific data::

def find_string(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append(line.rstrip())
                
    return list_of_results


#--- Cleaning matched data for export::

def clean_up(in_string):
 
    split = in_string.split(" ")
    cleanlist = []
    for elem in split:
        if elem != "":
            cleanlist.append(elem)
    return cleanlist


#--- Creating csv file::

def csv_generate(input_list, col_names, file_name):

    with open(temp_direct + file_name + '.txt', 'w')as out_file:
        writer = csv.writer(out_file)
        count = 0
        col_n = len(col_names)
        row_lines = [None] * col_n 

        for line in input_list:
        
            if count == 0:
                writer.writerow(col_names)
            count+=1

            for i in range(col_n):
                row_lines[i] = line[i]  
        
            writer.writerow(row_lines)

    out_file.close()

#--- Combining meteo data from several files::
def combine(f, col_names, n_days):
    with open(out_direct + 'meteo_data.out', 'w') as output_file:
        writer = csv.writer(output_file)
        
        
        for i in range(n_days):
            if i == 0:
                writer.writerow(col_names)
                file = open(temp_direct + f + '.txt', 'r')
                start_day = datetime.strptime(f, '%Y%m%d')
                start_day = (start_day + timedelta(days=1)).strftime('%Y%m%d')

            else:
                file = open(temp_direct + start_day + '.txt', 'r')
                start_day = datetime.strptime(start_day, '%Y%m%d')
                start_day = (start_day + timedelta(days=1)).strftime('%Y%m%d')
            


            reader = csv.DictReader(file)
                
            for line in reader:
                rowlines = []
    
                rowlines.append(line['Date'])
                rowlines.append(line['Hour'])
                rowlines.append(line['Temp [ºC]'])
                rowlines.append(line['Humidity [%]'])
                rowlines.append(line['PNM [hPa]'])
                rowlines.append(line['DD [gr]'])
                rowlines.append(line['FF [km/hr]'])
                writer.writerow(rowlines)
        
            
    
    output_file.close()
    file.close()





#Downloading daily meteorological data (one file)

url = 'https://ssl.smn.gob.ar/dpd/zipopendata.php?dato=datohorario'
r = requests.get(url)
z = zipfile.ZipFile(io.BytesIO(r.content))
result_file = z.extractall()



#Extracting San Juan info and combining downloaded meteo files

direct = '~ownCloud/Evaporation Model/San Juan/'


day = input("Enter start date (YYYYMMDD):\n")
n_days = int(input("Enter number of days:\n"))
start_date = day

col_names = ['Date','Hour','Temp [ºC]','Humidity [%]', 'PNM [hPa]', 'DD [gr]', 'FF [km/hr]']



for j in range(n_days):
    f_name = 'datohorario' + day
    matched_lines = find_string(in_direct + f_name + '.txt', 'SAN JUAN AERO')
    finallist = []
    for i in range(len(matched_lines)):
        temp = clean_up(matched_lines[i])
        finallist.append(temp)
    csv_generate(finallist, col_names, day)
    day = datetime.strptime(day, '%Y%m%d')  
    day = (day + timedelta(days=1)).strftime('%Y%m%d')
    
combine(start_date, col_names, n_days)



