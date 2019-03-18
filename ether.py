from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt

fig, ax = subplots()
ax.set_title( 'scale: {:<.2e}'.format(1.0) )
subplots_adjust(left=0.25, bottom=0.25)
l, = plot([1], [0], 'r.',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axS = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axM = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axC = axes([0.25, 0.07, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
axSrM = axes([0.05, 0.085, 0.15, 0.06])
axSrMbutt = axes([0.05, 0.045, 0.15, 0.04])
# SrMtxt = txt.Text(text = '1.0')

slS =  Slider(axS, 'S', -pi, pi, valinit=pi/2)
slM =  Slider(axM, 'M', -4, 4, valinit=1)
slC =  Slider(axC, 'C', 0., 1, valinit=.5)
slN =  Slider(axN, 'N', 0., 4, valinit=2.)
tbSrM = TextBox(axSrM, 'S/M')
bSrM = Button(axSrMbutt, 'inverse')
SrMset = False
SrMinv = False

def draw(M,N,S,C):
    rN =range(N)
    a = zeros((N,),dtype = complex)
    an = 1
    mabs = 1
    for n in rN:
        a[n] = an
        if abs(an)> mabs:
            mabs = abs(an)
        an = an *(S + 2*C*sin(n*M))
    l.set_data(a.real/mabs, a.imag/mabs)
    ax.set_title('scale: {:<.2e}'.format(mabs) )
    ax.axis('scaled')
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    fig.canvas.draw_idle()

def update_n(val):
    if(SrMset):
        update_rel(val)
    else:
        update(val)

    
def update(val):
    global SrMset
    if(SrMset and not SrMinv):
        slM.disconnect(0)
        slM.cnt = 0
        slM.on_changed(update)
    elif(SrMset):
        slS.disconnect(0)
        slS.cnt = 0
        slS.on_changed(update)
    argS = slS.val
    C = slC.val
    M = slM.val
    N = int(exp(slN.val*log(1e1)))
    S = cos(argS)*(1+C*C) + (1-C*C)*1j*sin(argS)
    draw(M,N,S,C)
    SrMset = False

def update_txt(val):
    global SrMset
    if(not SrMinv):
        slM.disconnect(0)
        slM.cnt = 0
        cid = slM.on_changed(update_rel)
    else:
        slS.disconnect(0)
        slS.cnt = 0
        slS.on_changed(update_rel)
    update_rel(val)
    SrMset = True

def update_rel(val):
    global SrMset
    C = slC.val
    if(not SrMinv):
        M = slM.val
        argS = float(tbSrM.text)*M
        SrMset = False
        slS.set_val(argS)
        SrMset = True
    else:
        argS = slS.val
        M = float(tbSrM.text)*argS
        SrMset = False
        slM.set_val(M)
        SrMset = True
    N = int(exp(slN.val*log(1e1)))
    S = cos(argS)*(1+C*C) + (1-C*C)*1j*sin(argS)
    draw(M,N,S,C)
    
def inverse_rel(val):
    global SrMinv
    SrMinv = not SrMinv
    tbSrM.label.set_visible(False)
    tbSrM.text_disp.set_visible(False)
    initial = tbSrM.text
    if(not SrMinv):
        tbSrM.__init__(axSrM, 'S/M', initial = initial)
        tbSrM.on_submit(update_txt)
        slS.disconnect(0)
        slS.cnt = 0
        slS.on_changed(update)
    else:
        tbSrM.__init__(axSrM, 'M/S', initial = initial)
        tbSrM.on_submit(update_txt)
        slM.disconnect(0)
        slM.cnt = 0
        slM.on_changed(update)
    try: 
        update_txt(val) 
    except:
        pass 

# print()
bSrM.on_clicked(inverse_rel)
slS.on_changed(update)
slN.on_changed(update_n)
slM.on_changed(update)
slC.on_changed(update_n)
tbSrM.on_submit(update_txt)
show()