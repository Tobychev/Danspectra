#encoding: utf8
import time
from wipfile import * 
from matplotlib import animation
import scipy.interpolate as intr
import scipy.signal as sg

def linevel(lam):
    c = 299792.458 # km/s
    return c*(lam - line.cent)/line.cent

s6405_t5p.normalize()
frm = s6405_t5p.frames[0]
line = CNq  

mes = myst.measure(s6405_t5p)
frmcon = mes[2][:,:798]

sdata = frm.data[np.argsort(frmcon),:].reshape(798,-1)

frmcon.sort()

lam  = frm.group.lmbd
lbox3 = np.convolve(lam[line.idx],np.ones(3)/3.,"valid")
lbox3 = linevel(lbox3)

lsav73 = sg.savgol_filter(lam[line.idx],7,3)
lsav73 = linevel(lsav73)


fig = pl.figure()
ax = pl.axes(
            xlim=(lbox3[-1], lbox3[0]),
#            xlim=(lsav73[-1], lsav73[0]),
#            xlim=(lam[line.idx[-1]], lam[line.idx[0]]),
            ylim=(0.914,1.03 ))
#            ylim=(-0.1,.1 ))

curve, = ax.step([], [], lw=2)
text   = ax.text(0.04, 0.92, '',fontsize=16, transform=ax.transAxes)
pause = False
def animinit():
    text.set_text('')
    curve.set_data([], [])
    return curve,

def onClick(event):
    global anim,pause
    pause = not pause
    if pause:
        anim.event_source.stop()
    else:
        anim.event_source.start()

def animate(i):
    data = sdata[i,line.idx]
    ys = np.convolve(data,np.ones(3)/3.,"valid"); xs = lbox3; method="rolling mean 3 points"
#    spli = intr.splrep(np.flipud(lam[line.idx]) ,np.flipud(data),s=0.007)  #Spline interpolation, smoothing hand-chosen, 0.0055 good
#    ys = intr.splev(lam[line.idx],spli); xs = lam[line.idx]; method="spline smoot s=0.007" 
#    ys = sg.wiener(data,mysize=3) - intr.splev(lam[line.idx],spli); xs = lam[line.idx] 
#    ys = sg.wiener(data,mysize=3); xs = lam[line.idx]; method="Wiener filter 3 points" 
#    ys = sg.savgol_filter(data[i,line.idx],7,3); xs = lsav73; method="Savgol filter p=7 n=3 "
#    ys = data; xs = lam[line.idx] 
    curve.set_data(xs, ys)
    text.set_text('Row {}, con {:.4f}, method {} '.format(i+1,frmcon[0,i],method) )
    return curve,text


fig.canvas.mpl_connect('button_press_event', onClick)
anim = animation.FuncAnimation(fig, animate, init_func=animinit,
                    frames=798, interval=325, blit=True)
#anim.save('frame_aS1_381_rollmean3.mp4', fps=3)
pl.show()
