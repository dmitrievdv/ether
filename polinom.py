from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from fractions import gcd

fig, ax = subplots()
ax.set_title( 'B(t)' )
# ax = axes[0]

subplots_adjust(left=0.25, bottom=0.25) 
brim, = plot([1], [0], '.', color = 'royalblue', markersize=3)
pol_brim, = plot([1], [0], 'r.',markersize=1)
ax.axis('scaled')
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)

# axS = axes([0.25, 0.16, 0.65, 0.03], facecolor='gray')
axrt = axes([0.25, 0.15, 0.65, 0.03], facecolor='grey')
axpt = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.07, 0.65, 0.03], facecolor='grey')
axrN = axes([0.25, 0.03, 0.65, 0.03], facecolor='grey')
axSrMbutt = axes([0.05, 0.125, 0.15, 0.04])

# slS =  Slider(axS, 'S', -pi, pi, valinit=pi/2)
slrt =  Slider(axrt, r'$r$', 0, 3, valinit=0.5)
slpt =  Slider(axpt, r'$\phi$', -pi, pi, valinit=0)
slN = Slider(axN, 'N', 0, 20, valinit=4)
slrN = Slider(axrN, 'one roots degree', 1, 180, valinit=20)
bSrM = Button(axSrMbutt, r'$t/s$ switch')

R = 1

def draw(val):
    pt = slpt.val
    rt = slrt.val
    rN = int(floor(slrN.val))
    N = int(floor(slN.val))
    one_roots = [1.0]
    for k in range(1,rN):
        for m in range(1,k):
            if(gcd(m,k) == 1):
                one_roots.append(exp(1.0j*2*m*pi/k))
    if R:
        ax.set_title( 'B(t)' )
        t = rt*2*(cos(pt)+1j*sin(pt))
        fund_roots = array([-t*sin(2*m/N*pi) 
                           for m in range(N)], 
                           dtype = 'complex128')
        pol = poly1d(fund_roots, r=True)
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
        if rt<1:
            P = (rN*rN*N)//4
            for k in range(P+1):
                argS = k/P*pi*2
                x = cos(argS)*(1+rt*rt)
                y = (1-rt*rt)*sin(argS)
                brim_real.append(x*cos(pt)-y*sin(pt)) 
                brim_imag.append(y*cos(pt)+x*sin(pt))

        pol_brim.set_data(roots.real,roots.imag)
        brim.set_data(brim_real, brim_imag)
        # ll.set_data
    else:
        # rt = rt*2
        ax.set_title( 'F(s)' )
        s = rt*(cos(pt)+1j*sin(pt))
        
        pol = poly1d([1])
        for m in range(N):
            pol = pol*poly1d(array(
                [2*sin(m*2*pi/N),s], 
                dtype = 'complex128')
            )
        roots = []
        for one in one_roots:
            one_pol = poly1d([one])
            frill_pol = one_pol - pol
            roots.extend(frill_pol.r)
        # print(abs(pol(S)))
        # print(S)
        roots = array(roots)
        pol_brim.set_data(roots.real,roots.imag)
        # print(abs(pol(roots)))
        frill_real=[]
        frill_imag=[]
        if rt<2:
            P = (rN*rN*N)//16
            ca = cos(pt)
            sa = sin(pt)
            for j in range(P): 
                rho = sqrt(abs(rt-1)) + (1-sqrt(abs(rt-1))) *j/P
                rho2 = rho*rho
                rho4 = rho2*rho2
                rs2 = rt*rt
                c2 = rs2 - 1 + (2 + rs2)*rho4 - rho4*rho4
                c2 = c2/rs2/rho2
                if abs(c2)<=2:
                    c = rho*sqrt(c2 + 2)/2
                    s = rho*sqrt(2 - c2)/2
                    frill_real.append( ca*c - sa*s)
                    frill_imag.append( sa*c + ca*s)
                    frill_real.append( sa*s - ca*c)
                    frill_imag.append(-sa*c - ca*s)
                    frill_real.append( ca*c + sa*s)
                    frill_imag.append( sa*c - ca*s)
                    frill_real.append(-sa*s - ca*c)
                    frill_imag.append( ca*s - sa*c)
        brim.set_data(frill_real,frill_imag)
        

    fig.canvas.draw_idle()




def switch(val):
    global R
    if R:
        R = False
    else:
        R = True 
    draw(0)


draw(0)
bSrM.on_clicked(switch)
slrt.on_changed(draw)
slpt.on_changed(draw)
slrN.on_changed(draw)
slN.on_changed(draw)
draw(1)
show()
