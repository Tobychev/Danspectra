#Plot trying to show development of the unkonwn line in 3d, didn't work
if True:
    lbox3 = np.convolve(lam[line.idx],np.ones(3)/3.,"valid")
    rows = range(0,frm.data.shape[0])

    #Xs,Ys = np.meshgrid(lbox3,np.linspace(0,2000,num=798))
    #Zs = []
    for row in rows[0:400:40]:
    #    Zs.append( np.convolve(frm.data[row,line.idx],np.ones(3)/3.,"valid") ) ,zs=row,,
        ys = np.convolve(frm.data[row,line.idx],np.ones(3)/3.,"valid")
        ax.step(lbox3,ys,row/(400*1.),zdir='y')

    #ax.plot_wireframe(Xs, Ys, Zs, rstride=20, cstride=30)
    #ax.plot_surface(Xs, Ys, Zs, rstride=8, cstride=80, alpha=0.3)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim3d(640.5, 640.6)
    ax.set_ylim3d(0.8,1.2)
    ax.set_zlim3d(0.8,1.2)
    pl.show()
