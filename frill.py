from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
 
fig, ax = subplots()
ax.set_title( 'F(s)' )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'b-',markersize=1.)
ll, = plot([1], [0], 'b-',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axaC = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axrC = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
# SrMtxt = txt.Text(text = '1.0')

slaC =  Slider(axaC, r'$\varphi_s$', -pi, pi, valinit=0)
slrC =  Slider(axrC, r'$\rho_s$', 1e-7, 2, valinit=1)
slN =  Slider(axN, 'N', 1., 3, valinit=2)

def draw(N,aC,rC):         

    try:
        e = 1/0 # raise an error, hehe >:D
        from ethercalc import compute
        x_f, y_f, n_xy = compute.frill(N, aC, rC) 
        # print(frill.compute.__doc__)
        x = x_f[:n_xy-1]; y = y_f[:n_xy-1]
    except:
        x1 =[]
        y1 =[]
        x2 =[]
        y2 =[]
        x3 =[]
        y3 =[]
        x4 =[]
        y4 =[]
        ca = cos(aC)
        sa = sin(aC)
        for j in range(N): 
            rho = sqrt(abs(rC-1)) + (1-sqrt(abs(rC-1))) *j/N
            rho2 = rho*rho
            rho4 = rho2*rho2
            rs2 = rC*rC
            c2 = rs2 - 1 + (2 + rs2)*rho4 - rho4*rho4
            c2 = c2/rs2/rho2
            if abs(c2)<=2:
                c = rho*sqrt(c2 + 2)/2
                s = rho*sqrt(2 - c2)/2
                x1.append( ca*c - sa*s)
                y1.append( sa*c + ca*s)
                x3.append( sa*s - ca*c)
                y3.append(-sa*c - ca*s)
                x2.append( ca*c + sa*s)
                y2.append( sa*c - ca*s)
                x4.append(-sa*s - ca*c)
                y4.append( ca*s - sa*c)
            

        x = x1+[ca]+ x2[::-1]+[0 if rC ==1 else x1[0] if rC >1 else x3[0]]
        xx = x3+[-ca]+x4[::-1]+[0 if rC ==1 else x3[0] if rC >1 else x1[0]]
        y = y1+[sa]+y2[::-1]+[0 if rC ==1 else y1[0] if rC >1 else y3[0]]
        yy = y3+[-sa]+y4[::-1]+[0 if rC ==1 else y3[0] if rC >1 else y1[0]]
    ll.set_data(x,y)
    l.set_data(xx,yy)
    
    
    fig.canvas.draw_idle()

    
def update(val):
    aC = slaC.val
    rC = slrC.val
    N = int(exp(slN.val*log(1e1)))
    draw(N,aC,rC)

update(0)
slrC.on_changed(update)
slaC.on_changed(update)
slN.on_changed(update)
show()


