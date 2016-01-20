#encoding: utf-8
import pyfits as f
import os
import re

Dir = "data/"

def parse_dan_filename(name):
    reg = (r"(?P<dag>21|22)"
           r"(?P<manad>Sep|Jul)99_"
           r"(?P<angstrom>\d{4})"
           r"(?P<serie>[a-z])"
           r"(?:_|m16b_im(?P=dag)(?P=manad)99\.)"
           r"(?P<last>\w+)")
    pat = re.compile(reg)
    dag,manad,angstrom,serie,last = pat.match(name).groups()
    if last.isdigit():
        runnr = int(last)
        typ   = "cor"
    else:
        typ   = last
        runnr = ""
    return (manad,dag,typ,serie,angstrom,runnr)


def rename_data():
    files = os.listdir(Dir)
    files.sort()
    for name in files:
        (manad,dag,typ,serie,angstrom,runnr) = parse_dan_filename(name)
        new_series = serie+manad[0]+dag[-1]
        new_name   = "{}_{}_{}_{}.fits".format(angstrom,new_series,runnr,typ)
        print "Filnamn: {:<50} Nytt filnamn {:<50}".format(Dir+name,Dir+new_name)
        print "Renaming disabled for safety"
        #os.rename(Dir+name,Dir+new_name)

