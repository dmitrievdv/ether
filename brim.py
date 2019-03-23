from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt


fig, ax = subplots()
ax.set_title( 'B(t)' )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'r.',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

axaC = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axrC = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axM = axes([0.25, 0.07, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
axSrMbutt = axes([0.05, 0.025, 0.15, 0.04])
# SrMtxt = txt.Text(text = '1.0')

slaC =  Slider(axaC, r'$\varphi_t$', -pi, pi, valinit=0)
slrC =  Slider(axrC, r'$\rho_t$', 0, 4, valinit=1)
slM =  Slider(axM, r'$\mu$', 0, 2*pi, valinit=1)
slN =  Slider(axN, 'N', 1., 4, valinit=2.)
bSrM = Button(axSrMbutt, r'$\mu$ switch')

R=[]

def draw():
    aC = slaC.val
    rC = slrC.val
    N = int(exp(slN.val*log(1e1)))

    C = rC*(cos(aC) + 1j*sin(aC))
    if R :
        mu = slM.val
        k = 1
        for i in range(1,50):
            mup=mu*i/2/pi
            if abs(int(mup)-mup)*50<1:
                k = i
                break
        P = poly1d([1])
        for i in range(k):
            P = P*poly1d([1,C*sin(i*2*pi/k)])
        s = []
        Po = P[0]
        for i in range(N):
            P[0]=Po+exp(1j*2*pi/(i+1))
            s.extend(P.r)
        s = array(s)
        sc = max(abs(s))/2

        l.set_data(s.real/sc,s.imag/sc)
        ax.set_title('B(t), den: {:2}, scale: {:<.2e}'.format(k,sc) )

    else:

        P = N+1
        piss = zeros((P,))
        puss = zeros((P,))
        

        for k in range(P):
            argS = k/N*pi*2
            x = cos(argS)*(1+rC*rC/4)
            y = (1-rC*rC/4)*sin(argS)
            puss[k]=x*cos(aC)-y*sin(aC) 
            piss[k]=y*cos(aC)+x*sin(aC)

        l.set_data(puss,piss)
    fig.canvas.draw_idle()

    
def update(val):
    draw()


def inverse_rel(val):
    global R
    if R:
        R = False
    else:
        R = True 
    draw()


update(0)
bSrM.on_clicked(inverse_rel)
slrC.on_changed(update)
slaC.on_changed(update)
slN.on_changed(update)
slM.on_changed(update)
show()

