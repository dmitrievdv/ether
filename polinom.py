from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt

fig, ax = subplots()
ax.set_title( 'B(t)' )
# ax = axes[0]

subplots_adjust(left=0.25, bottom=0.25) 
brim, = plot([1], [0], 'b--', linewidth=0.5)
pol_brim, = plot([1], [0], 'r.',markersize=1.5)
ax.axis('scaled')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)

# axS = axes([0.25, 0.16, 0.65, 0.03], facecolor='gray')
axrt = axes([0.25, 0.15, 0.65, 0.03], facecolor='grey')
axpt = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.07, 0.65, 0.03], facecolor='grey')
axrN = axes([0.25, 0.03, 0.65, 0.03], facecolor='grey')

# slS =  Slider(axS, 'S', -pi, pi, valinit=pi/2)
slrt =  Slider(axrt, r'$t_r$', 0, 1, valinit=0.5)
slpt =  Slider(axpt, r'$t_{\phi}$', -pi, pi, valinit=0)
slN = Slider(axN, 'N', 0, 100, valinit=10)
slrN = Slider(axrN, 'one roots degree', 1, 20, valinit=10)

def draw(val):
    pt = slpt.val
    rt = slrt.val
    rN = int(floor(slrN.val))
    # argS = pi/2 - pt
    # z = rt**2
    t = rt*2*(cos(pt)+1j*sin(pt))
    # x = cos(argS)*(1+z)
    # y = (1-z)*sin(argS)
    # S = x*cos(pt)-y*sin(pt) + 1j*(y*cos(pt)+x*sin(pt))
    one_pol = poly1d([1.0])
    N = int(floor(slN.val))
    fund_roots = array([-t*sin(2*m/N*pi) 
                       for m in range(N)], 
                       dtype = 'complex128')
    pol = poly1d(fund_roots, r=True)
    one_roots = [1.0]
    for k in range(1,rN):
        for m in range(1,k):
            if(gcd(m,k) == 1):
                one_roots.append(exp(1.0j*2*m*pi/k))
    roots = []
    for one in one_roots:
        one_pol = poly1d([one])
        brim_pol = one_pol - pol
        roots.extend(brim_pol.r)
    # print(abs(pol(S)))
    # print(S)
    roots = array(roots)
    # print(abs(pol(roots)))
    brim_real = []; brim_imag = []
    P = 300
    for k in range(P+1):
        argS = k/P*pi*2
        x = cos(argS)*(1+rt*rt)
        y = (1-rt*rt)*sin(argS)
        brim_real.append(x*cos(pt)-y*sin(pt)) 
        brim_imag.append(y*cos(pt)+x*sin(pt))

    pol_brim.set_data(roots.real,roots.imag)
    brim.set_data(brim_real, brim_imag)
    # ll.set_data
    fig.canvas.draw_idle()


slrt.on_changed(draw)
slpt.on_changed(draw)
slrN.on_changed(draw)
slN.on_changed(draw)
draw(1)
show()
# axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
# axSrMbutt = axes([0.05, 0.025, 0.15, 0.04])
