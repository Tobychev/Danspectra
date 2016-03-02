#encoding: utf8
from wipfile import * 
from matplotlib import animation

def linevel(lam):
    c = 299792.458 # km/s
    return c*(lam - line.cent)/line.cent


lbox3 = np.convolve(lam[line.idx],np.ones(3)/3.,"valid")
lbox3 = linevel(lbox3)
fig = pl.figure()
ax = pl.axes(
            xlim=(lbox3[0], lbox3[-1]),
            ylim=(0.94,1.06 ))

curve, = ax.step([], [], lw=2)
text   = ax.text(0.04, 0.92, '',fontsize=16, transform=ax.transAxes)

def animinit():
    text.set_text('')
    curve.set_data([], [])
    return curve,

def animate(i):
    ys = np.convolve(frm.data[i,line.idx],np.ones(3)/3.,"valid")
    curve.set_data(lbox3, ys/ys.mean())
    text.set_text('Row {}'.format(i+1) )
    return curve,text


anim = animation.FuncAnimation(fig, animate, init_func=animinit,
                    frames=798, interval=350, blit=True)
#anim.save('frame_aS1_381.mp4', fps=3)
pl.show()
