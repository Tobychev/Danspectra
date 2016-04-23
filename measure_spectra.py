import matplotlib.pyplot as pl
import visualize as vis
region = __import__("5654")

as2res = {} 
for line in region.as2lines:
    as2res[line.name.split(" ")[0]] = line.measure(region.as2)
