#encoding: utf8
from wipfile import * 
import numpy.polynomial.polynomial as pol
s6405_seg = dan.frameseries("data/6405_aS1","segments")

#####
#
# Conclusions
# ===========
# Just do a global fit and create an error statistic,
# then discard outliers that above some value of the 
# error statistic because the next level of sofistication
# will not perform much better, and have their own failure 
# modes
#
#####



mes = {}
mes["FeI top 5%"]  = FeI.measure_linecores(s6405_t5p)
mes["FeI segm"  ]  = FeI.measure_linecores(s6405_seg)
mes["Myst top 5%"] = myst.measure_linecores(s6405_t5p)
mes["Myst segm"  ] = myst.measure_linecores(s6405_seg)

mes["SiFe top 5%"]  = SiFe.measure_linecores(s6405_t5p)
mes["SiFe segm"  ]  = SiFe.measure_linecores(s6405_seg)

if False: #SiFe scatter
    pl.plot(mes["SiFe top 5%"][lin.lc.cont,:],mes["SiFe top 5%"][lin.lc.lbot,:],'bo',alpha=0.2);
    pl.title("Bottom, SiFe " + str(SiFe))
    pl.ylabel("Relative intensity of line bottom")
    pl.xlabel("Continuum value at line centre")
    pl.show()

    # Outliers
    # [332] [ 0.93543975] [ 640.71571066] Problem is that they have flat bottoms
    # [795] [ 1.10894355] [ 640.73377294] so that the narrow fit ends up poorly constrained in the
    #                                     spectral direction. Wonder how smoothing would work on that. 

#mes["SiFe top 5% by line"]  = SiFe.byline_measure_linecores(s6405_t5p)
#intr.select_linecore(mes["SiFe top 5% by line"])

#vis.show_line_and_corefit(SiFe,s6405_t5p.frames[0],331)

cv = 299792.458
line  = lins[4]
frame = s6405_seg.frames[0]
lmbd  = frame.group.lmbd

#####
#  Defining limits of fitting window
#####
width = len(line.idx)*0.16 # Fraction of points to be used
if width%2 == 0:
    width +=1
dwn = int(width - 1)/2; up = dwn+1
print width,dwn,up


#####
#  Making a global fit
#####

nrows = frame.data.shape[0]
guess   = frame.group.ref[line.idx].argmin()
bottom  = line.idx[guess +np.arange(-dwn,up)]
test    = line.idx[guess +np.arange(-(dwn+1),(up+1))]

fit     = pol.polyfit(lmbd[bottom],frame.data[:,bottom].T,2)
a,b,c   = fit[2,:],fit[1,:],fit[0,:]

lam_min = -b/(2*a)
lin_bot = pol.polyval(lam_min,fit,tensor=False)
pred    = pol.polyval(lmbd[test],fit)

vel = cv*(lam_min-line.cent)/line.cent
bot = lin_bot
con = frame.cont.val(lam_min) 
err = np.sqrt( np.mean( (frame.data[:,test]-pred)**2,axis=1) + 2*(lam_min-line.cent)**2 )


#####
# Making a line by line fit over entire frame 
#####

out   = np.zeros((4,nrows))
if False:
    for row in range(0,frame.data.shape[0]):
        cent    = frame.data[row,line.idx].argmin()

        if up > len(line.idx)-cent-1:
            bottom2  = line.idx[cent+np.arange(-dwn,up-2)]
            test2    = line.idx[cent+np.arange(-(dwn+1),(up-1))]
        else:   
            bottom2  = line.idx[cent+np.arange(-dwn,up)]
            test2    = line.idx[cent+np.arange(-(dwn+1),(up+1))]
        fit   = pol.polyfit(lmbd[bottom2],frame.data[row,bottom2],2)
        a,b,c = fit[2],fit[1],fit[0]
        lam_min = -b/(2*a)
        lin_bot = pol.polyval(lam_min,fit,tensor=False)
        prd    =  pol.polyval(lmbd[test2],fit)

        out[0,row] = lam_min
        out[1,row] = lin_bot
        out[2,row] = frame.cont.val(lam_min)[row]
        out[3,row] = np.sqrt( np.mean( (frame.data[row,test2]-prd)**2) + 2*(lam_min-line.cent)**2)

#####
# Making a line by line fit over parts of frame and printing diagnostics
#####

if True:
    for row in range(740,746):
        cent    = frame.data[row,line.idx].argmin()
        if up > len(line.idx)-cent-1:
            bottom2  = line.idx[cent+np.arange(-dwn,up-2)]
            test2    = line.idx[cent+np.arange(-(dwn+1),(up-1))]
        else:   
            bottom2  = line.idx[cent+np.arange(-dwn,up)]
            test2    = line.idx[cent+np.arange(-(dwn+1),(up+1))]
        fit   = pol.polyfit(lmbd[bottom2],frame.data[row,bottom2],2)
        a,b,c = fit[2],fit[1],fit[0]
        lam_min2 = -b/(2*a)
        lin_bot2 = pol.polyval(lam_min2,fit,tensor=False)
        prd    =  pol.polyval(lmbd[test2],fit)

        out[0,row] = lam_min2
        out[1,row] = lin_bot2
        out[2,row] = frame.cont.val(lam_min2)[row]
        out[3,row] = np.sqrt( np.mean( (frame.data[row,test2]-prd)**2) + 2*(lam_min2-line.cent)**2)

        print "Row {}".format(row)
        print "Global fit: cent = {}, bot = {}, RMSE = {}".format( lam_min[row], lin_bot[row],  err[row] )
        print "Line fit:   cent = {}, bot = {}, RMSE = {}".format(out[0,row], lin_bot2,  out[3,row])
        pl.step(lmbd[line.idx],frame.data[row,line.idx])
        pl.step(lmbd[bottom],frame.data[row,bottom])
        pl.plot(lmbd[line.idx][guess],frame.data[row,line.idx][guess],'o')
        pl.plot(lmbd[test2],prd,'r')
        pl.plot(lmbd[test],pred[row,:],'k')
        pl.show()

pl.plot(err)
pl.plot([0,800],(np.percentile(err,89))*np.ones(2))
pl.show()

pl.plot(out[3,:]);pl.plot([0,800],(err.mean()+2.5*err.std())*np.ones(2))
pl.show()

