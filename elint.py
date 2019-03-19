from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button

fig, ax = subplots()
ax.set_title( 'scale: {:<.2e}, sum: {:<.2e}'.format(1.0,0.0) )
subplots_adjust(left=0.15, bottom=0.20, right = 0.9) 
l, = plot([1], [0], 'k-')
zl, = plot([-pi,pi],[0,0],'b--')
ax.axis('scaled')
ax.set_xlim(-pi,pi)
ax.set_ylim(-2,2)

axaC = axes([0.15, 0.1, 0.7, 0.03], facecolor='gray')
axrC = axes([0.15, 0.05, 0.7, 0.03], facecolor='grey')

slaC =  Slider(axaC, 'x', 0, 1-1e-5, valinit=.5)
slrC =  Slider(axrC, 'z', 1e-5, 1-1e-5, valinit=.5)

N = 10000
def draw(z,x):
    th = linspace(0,pi*4,4*N)
    y = sqrt(((1-z**4)**2 -(1-z**2)**2*x*x )/(1+z**2)**2 )
    f = log((x+2*z*sin(th))**2+y*y)
    for i in range(N//2,2*N):
        if f[i+1]*f[i]<0: break
    sc = max(abs(f))/2
    l.set_data(th[i:i+2*N]-th[i]-pi,f[i:i+2*N]/sc)
    # ax.axis('scaled')
    # ax.set_xlim(-2,2)
    # ax.set_ylim(-2,2)
    ax.set_title( 'scale: {:<.2e}, integral: {:<.2e}'.format(sc,sum(f)*pi/N) )
    fig.canvas.draw_idle()

    
def update(val):
    z = slrC.val
    x = slaC.val*(1+z**2)
    draw(z,x)
update(0)
slrC.on_changed(update)
slaC.on_changed(update)
show()

