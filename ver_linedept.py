import matplotlib.pyplot as pl
import pickle as pic
import numpy as np

regnames = ["5053","5215","5654","6405","6449"]

for region in regnames:
    rg = __import__(region)
    for typ in ["qu1","qu2","spt"]:
        
        lines = zip( getattr(rg,typ+"lines"),pic.load(open("bin/{}_{}.lin".format(region,typ),"rb")) )
        print("Type "+typ)
        for load,lbin in lines:
            print("{:<18} {:.6f}. {:<18} {:.6f}".format(load.name,load.dept,lbin.name,lbin.dept))
    
