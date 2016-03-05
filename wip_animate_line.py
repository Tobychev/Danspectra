#encoding: utf8
from wipfile import * 
from matplotlib import animation

import scipy.signal as sg

def linevel(lam):
    c = 299792.458 # km/s
    return c*(lam - line.cent)/line.cent

lam  = frm.group.lmbd
lbox3 = np.convolve(lam[line.idx],np.ones(17)/17.,"valid")
lbox3 = linevel(lbox3)

lsav73 = sg.savgol_filter(lam[line.idx],7,3)
lsav73 = linevel(lsav73)

fig = pl.figure()
ax = pl.axes(
            xlim=(lbox3[-1], lbox3[0]),
#            xlim=(lsav73[-1], lsav73[0]),
#            xlim=(lam[line.idx[-1]], lam[line.idx[0]]),
            ylim=(0.914,1.03 ))

curve, = ax.step([], [], lw=2)
text   = ax.text(0.04, 0.92, '',fontsize=16, transform=ax.transAxes)

def animinit():
    text.set_text('')
    curve.set_data([], [])
    return curve,

def animate(i):
    ys = np.convolve(frm.data[i,line.idx],np.ones(17)/17.,"valid"); xs = lbox3
#    ys = sg.savgol_filter(frm.data[i,line.idx],7,3); xs = lsav73
#    ys = frm.data[i,line.idx]; xs = lam[line.idx] 
    curve.set_data(xs, ys)
    text.set_text('Row {}'.format(i+1) )
    return curve,text


anim = animation.FuncAnimation(fig, animate, init_func=animinit,
                    frames=798, interval=100, blit=True)
#anim.save('frame_aS1_381.mp4', fps=3)
pl.show()
