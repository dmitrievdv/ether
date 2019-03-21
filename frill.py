from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
 
fig, ax = subplots()
ax.set_title( 'F(s)' )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'r-',markersize=1.)
ll, = plot([1], [0], 'b.',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axaC = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axrC = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
# SrMtxt = txt.Text(text = '1.0')

slaC =  Slider(axaC, r'$\varphi_s$', -pi, pi, valinit=0)
slrC =  Slider(axrC, r'$\rho_s$', 0, 2, valinit=1)
slN =  Slider(axN, 'N', 1., 5, valinit=2)

def draw(N,aC,rC):         

    try:
        from ethercalc import compute
        x_f, y_f, n_xy = compute.frill(N, aC, rC) 
        # print(frill.compute.__doc__)
        x = x_f[:n_xy-1]; y = y_f[:n_xy-1]
    except:
        x =[]
        y =[]
        for j in range(N):
            sy = rC*(j*2/N-1+1/N)
            sx = sqrt(rC*rC-sy*sy)
            phi = arcsin(j*2/N-1+1/N)
            r = 1 - rC*rC
            q = 2*sx*sx- 2*sy*sy
            p = r-3
            P = poly1d([1,0,p,q,r])
            for c in P.r:
                if abs(c.imag)==0 and 0<=c.real<=1:
                    rho = sqrt(c.real)
                    x.append(rho*cos(aC+phi))
                    y.append(rho*sin(aC+phi))
                    x.append(-rho*cos(aC+phi))
                    y.append(-rho*sin(aC+phi))

    ll.set_data(x,y)
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


