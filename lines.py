import numpy as np
import scipy as sp
import danspec as dan
import interactive as intr

def make_pkwin_from_linegroup(lines):
    pkwin = [] 
    for itm in lines:
        pkwin.append(list(itm.win))

    return pkwin

def make_lines_from_wins(danspec,wins):
    lines = []
    for item in wins:
        lines.append(dan.line(danspec,item))

    return lines

def subdivide_line(line):
    mins =find_minima(line.ref)
    if len(mins) > 1:
        lines = []    
        mins  = intr.manual_delete_minima(mins,line)
        for li in mins:
            lines.append(dan.line( line.danspec, (line.idx[ li[0]] ,line.idx[li[2]]) ))

        return lines

    else:
        return line

def find_minima(curve):
    deriv = np.diff(curve)
    infl  = find_inflections(deriv)
    mins  = []
    for i,point in enumerate(infl):
        if is_minima(deriv,point):                      
            if i!= 0 and i < len(infl)-1:
                mins.append( (infl[i-1],point, infl[i+1]) )

    return mins

def find_inflections(deriv):
    mid =  np.where(np.diff(np.sign(deriv)))[0]+1
    return np.concatenate( ([0],mid,[len(deriv)]) )
    
def is_minima(curve,point):
    if curve[point-3:point].mean() < 0 and curve[point:point+3].mean() > 0:
        return True
    else:
        return False

def fit_linecores(line,xs,ys):
    guess = line.ref.argmin()
    bottom = line.idx[guess-4:guess+5]
    a,b,c = np.polyfit(xs[bottom],ys.T[bottom,:],2)
    return get_linecore(a,b,c) 

def fit_linecore(line,xs,ys):
    guess = line.ref.argmin()
    bottom = line.idx[guess-4:guess+5]
    a,b,c = np.polyfit(xs[bottom],ys[bottom],2)
    return get_linecore(a,b,c)    

def get_linecore(a,b,c):
    lam_min = -b/(2*a)
    lin_bot = np.polyval((a,b,c),lam_min)
    return lam_min,lin_bot
