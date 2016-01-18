#encoding: utf-8
import pandas as pd

data = pd.read_csv("list-of-data.csv",header=0)
wavebands = data["ångström"].unique(); wavebands.sort()

def find_series():
    ser = {}
    for wave in wavebands:
        ser[wave] = data["serie"][data["ångström"]==wave].unique()

    return ser

def count_series_files(wave,serie):
    return len(data[(data["ångström"]==wave) & (data["serie"]==serie)])

def make_metadata_rows(series_collection):
    for wave in series_collection.keys():
        for series in series_collection[wave]:
            count = count_series_files(wave,series)
            print "|   {}   |   {}   |  {:>3} |".format(wave,series,count)
 
