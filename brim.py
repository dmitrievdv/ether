# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from numpy.polynomial.legendre import leggauss
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from scipy.optimize import curve_fit as crvfit

fig, (sax, ax,) = subplots(1,2)
ax.set_title( 'scale: {:<.2e}'.format(1.0) )
sax.set_title( r'$\lambda(\sin(t))$' )
# subplots_adjust(left=0.25, bottom=0.25) 
ax.set_position([0.3,.25,.66,.66])
sax.set_position([0.05,.66,.22,.22])
# subplots_adjust(left=0.25, bottom=0.25) 
lp, = ax.plot([0], [0], 'b-',markersize=1.)
l, = ax.plot([1], [0], 'r.',markersize=2.) 
fl, = sax.plot([1], [0], 'k-',markersize=1.) 
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
sax.set_xlim(0,2*pi)
sax.set_xticks([0,pi/2,pi,1.5*pi,2*pi])
sax.set_xticklabels(['0',r'$\pi$',r'2$\pi$'])
sax.set_xticklabels(['0',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'2$\pi$'])

axS = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axP = axes([0.25, 0.1, 0.65, 0.03], facecolor='grey')
# axC = axes([0.25, 0.09, 0.65, 0.03], facecolor='grey')
# axF = axes([0.25, 0.05, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.05, 0.65, 0.03], facecolor='green')
axlf = axes([0.05, 0.5, 0.22, 0.06])
# axSrM = axes([0.05, 0.065, 0.15, 0.06])


slS =  Slider(axS, r'$\rho_t$', 0,1 , valinit=0.5)
slP =  Slider(axP, r'$\varphi_t$', -pi, pi, valinit=0)
# slC =  Slider(axC, r'$r$', 0, 2, valinit=1)
# slF =  Slider(axF, r'$\varphi$', -pi, pi, valinit=0)
slN =  Slider(axN, 'N', 0, 3, valinit=1.)
tblf = TextBox(axlf, r'$\lambda\,s\,:$')
# tbSrM = TextBox(axSrM, r'$\rho_{max}$')


# Func = lambda s: exp(s*s)*s/(s*s+1)
# Func = lambda s: s/(1-s*s+s*s*s*s)
Func = lambda s: s
rhom = 2



def draw(val):
    N = int(exp(slN.val*log(1e1)))
    C = slS.val*rhom#exp(slS.val)
    phis = linspace(0,pi/2,N)
    S = zeros((4*N,),dtype = complex)

    ns = linspace(0,2*pi,100)
    Fs = Func(sin(ns))
    lamax = amax(abs(Fs))
    fl.set_data(ns, Fs)
    sax.set_ylim(-lamax,lamax)
    sax.set_yticks([-1,0,1])


    def logrho1(ro,theta):
        return log((cf*ro+C*Func(sin(theta)))**2+sf*ro*sf*ro)


    for i in range(N):
        sv = 1e3
        dsv = sv/2
        Ng = 16
        Nl = 13
        phi = phis[i]
        cf = cos(phi)
        sf = sin(phi)
    
        gk, gw = leggauss(Ng)
        gk = (gk+1)/2
        lk, lw = linspace(0,2*pi,Nl+1,retstep=True)
        while (dsv/sv>1e-16 and sv>1e-16):
            s = 0
            for k in range(Nl):
                for j in range(Ng):
                    s+=logrho1(sv, lk[k]+lw*gk[j])*gw[j]
            
            if s>0:
                sv-=dsv
            else:
                sv+=dsv
            dsv/=2
        pS = sv*cf+1j*sv*sf
        cS = conj(pS)*exp(1j*slP.val)
        pS = pS*exp(1j*slP.val)
        S[i]=pS
        S[2*N+i]=-pS
        S[2*N-i-1]=-cS
        S[4*N-i-1]=cS

    # a = a/mabs
    l.set_data(S.real, S.imag)
    lp.set_data(S.real, S.imag)
  
    ax.axis('scaled')
    b = amax(abs(S))
    ax.set_xlim(-b,b)
    ax.set_ylim(-b,b)
    fig.canvas.draw_idle()
    return 0


def update_func(val):
    global Func
    Func = eval('lambda s : '+ tblf.text)
    draw(val)

# def update_rhom(val):
#     global rhom 
#     try:
#         rhom = float(tbSrM.text)
#     except:
#         rhom = 2
#     draw(val)

draw(1)
slS.on_changed(draw)
slP.on_changed(draw)
slN.on_changed(draw)

tblf.on_submit(update_func)
# tbSrM.on_submit(update_rhom)
# slC.on_changed(update)
# slF.on_changed(update)
show()
