import matplotlib.pyplot as pl
import matplotlib.cm as cm
import visualize as vis
region = __import__("5654")

sf_spot  = region.sf56_bs2
spotcon  = region.bs2con
spot     = region.bs2
umbra    = region.umbra
wall     = region.wall
penumbra = region.penumbra
quiet    = region.quiet
spotlines = region.bs2lines

#as2res = {} 
#for line in region.as2lines:
#    as2res[line.name.split(" ")[0]] = line.measure(region.as2)

spotum = sf_spot.make_spectra_subset(spot,rowsubset=( spotcon < umbra),desc="Umbra subset")
spotwl = sf_spot.make_spectra_subset(spot,rowsubset=((spotcon >= wall[0]) & (spotcon <= wall[1])),
                                                      desc="Umbra/penumbra wall")
spotpn = sf_spot.make_spectra_subset(spot,rowsubset=((spotcon > penumbra[0]) & (spotcon < penumbra[1])),desc="Penumbra")
spotqu = sf_spot.make_spectra_subset(spot,rowsubset=( spotcon >= quiet),desc="Quiet sun")

spotlist = [spotum,spotwl,spotpn,spotqu]


umres = {} 
colours = cm.Greys(np.linspace(0.2,0.7,6)); colours[1,:] = cm.Greys(1.)
umfig,axs = pl.subplots(3,3)
label = []
for i,line in enumerate(spotlines):
    name = line.name.split(" ")[0]
    umres[name] = line.measure(spotum)
    label.append(name)
    umfig = vis.addreg(umfig,umres[name],line,colour=colours[i])

wlres = {} 
for line in spotlines:
    wlres[line.name.split(" ")[0]] = line.measure(spotwl)
pnres = {} 
for line in spotlines:
    pnres[line.name.split(" ")[0]] = line.measure(spotpn)
qures = {} 
for line in spotlines:
    qures[line.name.split(" ")[0]] = line.measure(spotqu)
