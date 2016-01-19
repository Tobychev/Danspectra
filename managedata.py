#encoding: utf-8
import pandas as pd

data = pd.read_csv("list-of-data.csv",header=0)
wavebands = data["ångström"].unique(); wavebands.sort()

def find_series():
    ser = {}
    for wave in wavebands:
        ser[wave] = data["serie"][data["ångström"]==wave].unique()

    return ser

def count_series_files(wave,series):
    return len(data[(data["ångström"]==wave) & (data["serie"]==series) & (data["typ"] == "cor")])

def series_date(wave,series):
    return data["datum"][(data["ångström"]==wave) & (data["serie"]==series)].unique()

def normalize_series(wave,series,date):
    yr,mon,day = date.split("-")
    new_series = series+mon[0]+day[-1]
    data["serie"][(data["ångström"]==wave)&(data["serie"]==series)&(data["datum"] == date)] = new_series

def make_normalized_series(wave,series,dates):
    for date in dates:
        normalize_series(wave,series,date) 

def make_metadata_rows():
    series_collection = find_series()
    for wave in series_collection.keys():
        for series in series_collection[wave]:
            dates = series_date(wave,series)
            make_normalized_series(wave,series,dates)

    
    series_collection = find_series()
    for wave in series_collection.keys():
        for series in series_collection[wave]:
            count = count_series_files(wave,series)
            yr,mon,day = series_date(wave,series)[0].split("-")
            print "| {} {} |   {}   |   {}   |  {:>3} |".format(day,mon,wave,series,count)
 
