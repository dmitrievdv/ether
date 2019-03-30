# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt

fig, ax = subplots()
ax.set_title( 'scale: {:<.2e}'.format(1.0) )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'r.',markersize=1.) 
lp, = plot([1], [0], 'b.',markersize=2.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axS = axes([0.25, 0.17, 0.65, 0.03], facecolor='gray')
axM = axes([0.25, 0.13, 0.65, 0.03], facecolor='grey')
axC = axes([0.25, 0.09, 0.65, 0.03], facecolor='grey')
axF = axes([0.25, 0.05, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.01, 0.65, 0.03], facecolor='green')
axSrM = axes([0.05, 0.065, 0.15, 0.06])
axSrMbutt = axes([0.05, 0.025, 0.15, 0.04])
axSpi = axes([0.05, 0.135, 0.15, 0.06])
axFpi = axes([0.05, 0.205, 0.15, 0.06])
axCpi = axes([0.05, 0.275, 0.15, 0.06])
axSvbutt = axes([0.05, 0.325, 0.15, 0.04])
# SrMtxt = txt.Text(text = '1.0')

slS =  Slider(axS, r'$\varphi$', -pi, pi, valinit=pi/2)
slM =  Slider(axM, r'$\mu$', -4, 4, valinit=1)
slC =  Slider(axC, r'$\rho_t$', 0., 1, valinit=.5)
slF =  Slider(axF, r'$\varphi_t$', 0, pi, valinit=pi/3)
slN =  Slider(axN, 'N', 0., 5, valinit=3.5)
tbSrM = TextBox(axSrM, r'$\varphi/\mu$')
bSrM = Button(axSrMbutt, 'inverse')
bSv = Button(axSvbutt, 'save')
tbSpi = TextBox(axSpi, r'$\pi/\varphi$')
tbFpi = TextBox(axFpi, r'$\pi/\varphi_t$')
tbCpi = TextBox(axCpi, r'$\rho_t$')
SrMset = False
SrMinv = False

def draw():
    N = int(exp(slN.val*log(1e1)))
    F = slF.val
    M = slM.val
    argS = slS.val -F
    z = slC.val**2
    C = slC.val*2*(cos(F)+1j*sin(F))
    x = cos(argS)*(1+z)
    y = (1-z)*sin(argS)
    S = x*cos(F)-y*sin(F) + 1j*(y*cos(F)+x*sin(F))
    try: 
        from ethercalc import compute
        a, mabs = compute.fundamental(S, C, M, N)
    except:
        rN =range(N)
        a = zeros((N,),dtype = complex)
        an = 1
        mabs = 1
        for n in rN:
            a[n] = an
            if abs(an)> mabs:
                mabs = abs(an)
            an = an *(S + C*sin(n*M))
    a = a/mabs
    l.set_data(a.real, a.imag)
    # lp.set_data(a.real[::4], a.imag[::4])
    ax.set_title('scale: {:<.2e}'.format(mabs) )
    # ax.axis('scaled')
    # ax.set_xlim(-1,1)
    # ax.set_ylim(-1,1)
    fig.canvas.draw_idle()
    return a*mabs

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
    draw()
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
    draw()
    
def inverse_rel(val):
    global SrMinv
    SrMinv = not SrMinv
    tbSrM.label.set_visible(False)
    tbSrM.text_disp.set_visible(False)
    initial = tbSrM.text
    if(not SrMinv):
        tbSrM.__init__(axSrM, r'$\varphi/\mu$', initial = initial)
        tbSrM.on_submit(update_txt)
        slS.disconnect(0)
        slS.cnt = 0
        slS.on_changed(update)
    else:
        tbSrM.__init__(axSrM, r'$\mu/\varphi$', initial = initial)
        tbSrM.on_submit(update_txt)
        slM.disconnect(0)
        slM.cnt = 0
        slM.on_changed(update)
    try: 
        update_txt(val) 
    except:
        pass 

def update_spi(val):
    slS.set_val(pi/float(tbSpi.text))
    update(val)

 
def update_fpi(val):
    slF.set_val(pi/float(tbFpi.text))
    update(val)

 
def update_cpi(val):
    slC.set_val(float(tbCpi.text))
    update(val)

def save(val):
    a = draw()
    rN = slS.val/(2*pi)
    for M in range(1,20):
        if abs(int(rN*M)-rN*M)<1e-3:
            break
    print(M)
    cuta = a[::M]
    cuta.tofile('petal.dat')
        

# print()
update_n(1)
bSrM.on_clicked(inverse_rel)
bSv.on_clicked(save)
slS.on_changed(update)
slN.on_changed(update_n)
slM.on_changed(update)
slC.on_changed(update_n)
slF.on_changed(update_n)
tbSrM.on_submit(update_txt)
tbSpi.on_submit(update_spi)
tbCpi.on_submit(update_cpi)
tbFpi.on_submit(update_fpi)
show()
