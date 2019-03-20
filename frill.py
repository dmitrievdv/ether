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
slN =  Slider(axN, 'N', 1., 3, valinit=2)

def draw(N,aC,rC):

    sx = rC*cos(aC)
    sy = rC*sin(aC)

    x = []
    y = []

    for k in range(N):
        rho = k/N
        sa = 1+rho*rho
        sb = 1-rho*rho
        sc = sqrt(sa*sa-sb*sb)
        dif = 10
        pha = 0
        for j in range(N):
            phi = j*2*pi/N
            sa2 = sqrt((sc*cos(phi) - sx)**2 + (sc*sin(phi) - sy)**2)  + sqrt((sc*cos(phi) + sx)**2 + (sc*sin(phi) + sy)**2)
            if abs(sa2-2*sa)<1/4/N:
                x.append(rho*cos(phi))
                y.append(rho*sin(phi))

    ll.set_data(x,y)
    fig.canvas.draw_idle()

    
def update(val):
    aC = slaC.val
    rC = slrC.val
    N = int(exp(slN.val*log(1e1)))
    draw(N,aC,rC)

slrC.on_changed(update)
slaC.on_changed(update)
slN.on_changed(update)
show()

