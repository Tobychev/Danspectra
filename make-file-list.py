#encoding: utf-8
import pyfits as f
import os
import re

Dir = "data/"

def make_list_of_files():
    files = os.listdir(Dir)
    files.sort()
    with open("allfiles.list","w") as out:
        print "Found {} files".format( len(files) )
        line = "{:<8},{:<43}\n".format

        for fil in files:
            try:
                fits = f.getheader(Dir+"/"+fil)
                out.write( line(fits["SST_HEAD"],fil) )
            except:
                print "Error happened at:" + fil 

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

def parse_filenames():
    with open("allfiles.list","r") as infile:
        for line in infile:
            filename = line.split(",")[-1][:43].strip()
            (manad,dag,typ,serie,angstrom,runnr) = parse_dan_filename(filename)
            print """Datum 1999-{}-{}, typ: {:>10}, serie {}, Ångström = {}, nr: {}""".format(
                manad,dag,typ,serie,angstrom,runnr)

def parse_filenames_to_file(outfile):
    with open("allfiles.list","r") as infile, open(outfile,"w") as out:
        out.write("datum,typ,serie,ångström,nr,filnamn\n")
        for line in infile:
            filename = line.split(",")[-1][:43].strip()
            (manad,dag,typ,serie,angstrom,runnr) = parse_dan_filename(filename)
            out.write(
                "1999-{}-{},{},{},{},{},{}\n".format(
                 manad,dag,typ,serie,angstrom,runnr,filename))
