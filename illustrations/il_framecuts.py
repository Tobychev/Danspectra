import numpy as np
import matplotlib.pyplot as pl

regnames = ["5053","5215","5654","6405","6449"]

def edgeplot(frame):
    pl.subplot(4,2,1)
    pl.title("First col")
    pl.plot(frame[:,0])
    pl.plot(frame[:,1])
    pl.subplot(4,2,2)
    pl.title("Last col")
    pl.plot(frame[:,-1])
    pl.plot(frame[:,-2])
    pl.subplot(4,2,3)
    pl.plot(frame[:,1])
    pl.plot(frame[:,2])
    pl.subplot(4,2,4)
    pl.plot(frame[:,-2])
    pl.plot(frame[:,-3])
    pl.subplot(4,2,5)
    pl.title("First row")
    pl.plot(frame[0,:])
    pl.plot(frame[1,:])
    pl.subplot(4,2,6)
    pl.title("Last row")
    pl.plot(frame[-1,:])
    pl.plot(frame[-2,:])
    pl.subplot(4,2,7)
    pl.plot(frame[1,:])
    pl.plot(frame[2,:])
    pl.subplot(4,2,8)
    pl.plot(frame[-2,:])
    pl.plot(frame[-3,:])
    pl.show()

def middleplot(frame,dwn=650,up=700):
    pl.plot(frame[:,:].mean(axis=1))
    pl.xlim(dwn,up)
    pl.show()

for regname in regnames:
    region = __import__(regname)
    qu1 = region.sf_qu1.rawstack()
    print(qu1.shape)
    edgeplot(qu1)
    middleplot(qu1)
    qu2 = region.sf_qu2.rawstack()
    print(qu2.shape)
    edgeplot(qu2)
    middleplot(qu2)
    spt = region.sf_spt.rawstack()
    print(spt.shape)
    edgeplot(spt)
    middleplot(spt)

