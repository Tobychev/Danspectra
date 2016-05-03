import numpy as np

regnames = ["5053","5215","5654","6405","6449"]

for regname in regnames:
    region = __import__(regname)

    spec1  = region.qu1
    lines1 = region.qu1lines
    res = {} 
    for i,line in enumerate(lines1):
        name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
        try:
            res[name] = line.measure(spec1)
        except Exception as err:
            print("Line {} failed to measure:".format(name))
            print(err)
    np.savez_compressed("bin/{}_qu1".format(regname),**res)

    spec2  = region.qu2
    lines2 = region.qu2lines
    res = {} 
    for i,line in enumerate(lines2):
        name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
        try:
            res[name] = line.measure(spec2)
        except Exception as err:
            print("Line {} failed to measure:".format(name))
            print(err)
    np.savez_compressed("bin/{}_qu2".format(regname),**res)

    spec3  = region.spt
    lines3 = region.sptlines
    res = {} 
    for i,line in enumerate(lines3):
        name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
        try:
            res[name] = line.measure(spec3)
        except Exception as err:
            print("Line {} failed to measure:".format(name))
            print(err)
    np.savez_compressed("bin/{}_spt".format(regname),**res)


