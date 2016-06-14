import matplotlib.pyplot as pl
import scipy.interpolate as si

regnames = ["5053","5215","5654","6405","6449"]

fg,ax = pl.subplots(1)

for regname in regnames:
    region = __import__(regname)
    myst  = region.qu1Myst
    spec  = region.qu1[:,:].mean(axis=0)
    myst.recenter(spec)
    lmbd  = region.qu1.lmbd[myst.idx]; lmbd = lmbd - myst.cent # Shifting
    ax.plot(lmbd,spec[myst.idx],label="{} line".format(int(regname)/10))

ax.legend(loc="best")
ax.set_xlabel("Distance from line center [nm]")
ax.set_ylabel("Relative inensity")
ax.set_ylim(0.7997473945700595, 1.0121068790234533)
ax.set_xlim(-0.053167457081707142, 0.04466243368945573)
fg.savefig("../thesis/figures/Mystery_comparison.png")
fg.show()
