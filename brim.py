from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt


fig, ax = subplots()
ax.set_title( 'B(t)' )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'r-',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

axaC = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axrC = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
# SrMtxt = txt.Text(text = '1.0')

slaC =  Slider(axaC, r'$\varphi_t$', -pi, pi, valinit=0)
slrC =  Slider(axrC, r'$\rho_t$', 0, 2, valinit=1)
slN =  Slider(axN, 'N', 1., 4, valinit=2.)

def draw(N,aC,rC):
    P = N+1
    piss = zeros((P,))
    puss = zeros((P,))
    C = rC*(cos(aC) + 1j*sin(aC))
    

    for k in range(P):
        argS = k/N*pi*2
        x = cos(argS)*(1+rC*rC/4)
        y = (1-rC*rC/4)*sin(argS)
        puss[k]=x*cos(aC)-y*sin(aC) 
        piss[k]=y*cos(aC)+x*sin(aC)

    l.set_data(puss,piss)
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

