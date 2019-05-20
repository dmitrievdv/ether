# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from scipy.optimize import curve_fit as crvfit
from numpy.polynomial.legendre import leggauss

fig, (sax, ax,) = subplots(1,2)
ax.set_title( 'scale: {:<.2e}'.format(1.0) )
sax.set_title( r'$\lambda\, s$' )
# subplots_adjust(left=0.25, bottom=0.25) 
ax.set_position([0.3,.3,.66,.6])
sax.set_position([0.0,.6,.33,.33])
# subplots_adjust(left=0.25, bottom=0.25) 
l, = ax.plot([1], [0], 'r.',markersize=2.) 
lp, = ax.plot([0], [0], 'b-',markersize=1.)
lc, = ax.plot([0], [0], 'g.',markersize=2.)
fp, = sax.plot([0,0], [0,0], 'g--',markersize=1.) 
fl, = sax.plot([1], [0], 'k-',markersize=1.) 
fc, = sax.plot([0], [0], 'r.',markersize=1.) 
ax.axis('scaled')
sax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
# sax.set_xlim(0,2*pi)
# sax.set_xticks([0,pi/2,pi,1.5*pi,2*pi])
# sax.set_xticklabels(['0',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'2$\pi$'])

axlf = axes([0.25, 0.21, 0.65, 0.04])
axX = axes([0.25, 0.17, 0.65, 0.025], facecolor='darkred')
axM = axes([0.25, 0.14, 0.65, 0.025], facecolor='darkred')
axP0 = axes([0.25, 0.105, 0.65, 0.025], facecolor='grey')
axP1 = axes([0.25, 0.075, 0.65, 0.025], facecolor='grey')
axP2 = axes([0.25, 0.045, 0.65, 0.025], facecolor='grey')
axN = axes([0.25, 0.01, 0.65, 0.025], facecolor='darkgreen')
axAP = axes([0.05, 0.06, 0.15, 0.05])
axXt = axes([0.05, 0.13, 0.05, 0.05])
axXm = axes([0.15, 0.13, 0.05, 0.05])
axP2t = axes([0.05, 0.20, 0.15, 0.05])
axP1t = axes([0.05, 0.27, 0.15, 0.05])
axP0t = axes([0.05, 0.34, 0.15, 0.05])
axMt = axes([0.05, 0.41, 0.15, 0.05])
axSPb = axes([0.05, 0.03, 0.15, 0.03]) 

slX =  Slider(axX, r'$\chi$', -pi,pi , valinit=-pi/2 ,valfmt = "%1.5f")
slM =  Slider(axM, r'$\nu$', -4, 4, valinit=-2,valfmt = "%1.5f")
slP0 =  Slider(axP0, r'$P0$', -1, 1, valinit=-.5,valfmt = "%1.5f")
slP1 =  Slider(axP1, r'$P1$', -1., 1, valinit=-.5,valfmt = "%1.5f")
slP2 =  Slider(axP2, r'$P2$', -1, 1, valinit=-.5,valfmt = "%1.5f")
slN =  Slider(axN, 'N', 3., 5, valinit=3.5)

bSP = Button(axSPb, 'save petal')

tbM = TextBox(axMt, r'$\nu$')
tbXt = TextBox(axXt, r'${\chi=\tau\times}$')
tbXm = TextBox(axXm, r'${+\nu\times}$')
tbP2 = TextBox(axP2t, r'$P2$')
tbP1 = TextBox(axP1t, r'$P1$')
tbP0 = TextBox(axP0t, r'$P0$')
tbAP = TextBox(axAP, r'$N_{ptls}$')
tblf = TextBox(axlf, r'$\lambda\,s\,:$')

# some functions :) 
# Func = lambda s: sin(s)# exp(1j*s)
# Func1 = lambda s: exp(sin(s)*sin(s))*sin(s)/(sin(s)*sin(s)+1)
# Func2 = lambda s: sin(s)/(1-sin(s)*sin(s)+sin(s)*sin(s)*sin(s)*sin(s))
# Func3 = lambda s: sin(s)/(.5+sin(s)*sin(s))
# Func4 = lambda s: (1-2*sin(s)*sin(s))/(2+sin(s)) 
Ellipse = lambda a, s: (a*a+1)*cos(s) + (1-a*a)*1j*sin(s)
Quellipse = lambda a, s: 4j*a**2*cos(s)**2 + (sqrt(4-2*a*a)-sqrt(2)*a)**2*sin(s)**2
# Func = lambda s: Ellipse(1,s)/(1-sin(s)*sin(s)+sin(s)*sin(s)*sin(s)*sin(s)) 


def petal(a):
    N = a.shape[-1]
    na = a - sum(a)/N
    sp = fft.fft(na)[1:N//2]
    pf = fft.fftfreq(N)[1:N//2]
    asp = abs(sp)

    profileL = lambda x, A, M, S: A/(1 + abs(x - M)/S)
    profileN = lambda x, A, M, S: A*exp(-(x - M)**2/S**2) 
    profile =  lambda x, A, M, S: A*(exp(-(x - M)**2/S**2) + 1/(1 + abs(x - M)/S))
    mix = []
    masp = amax(asp)
    for i in range(len(asp)):
        if asp[i]> masp/180:
            mix.append((i,asp[i]))

    width = pf[0]
    Ps = []
    for  mi,mx in mix:
        casp = asp[mi-5:mi+5]
        cupf = pf[mi-5:mi+5]
        p = crvfit(profile, cupf, casp, [mx, pf[mi], width] )[0]
        cupft=linspace(cupf[0],cupf[-1],200)

        Ps.append(abs(1/p[1]))

    #     plot(cupft,profile(cupft,p[0],p[1],p[2]),'b-')
    #     plot(cupf,casp,'r.')
    # plot(pf,asp,'-')
    # show()


    # print(Ps)
    SD = 2
    # print(N,SD,1/width)
    while width>N**(-1.4):
        liph = int(sqrt(1/width))
        # print(liph,width,N**(-.85))
        rL = 0
        rP = 0
        for P in Ps:
            iph = sort(array([ (i,(1.*i)%P) for i in range(liph)], dtype=[('ind', '<i4'),('phase', '<f8')] ),order='phase')
            L = 0
            pik = 0
            l = []
            w = []
            r = []
            for j in range(len(iph)):
                i2, ph2 = iph[j-2]
                i1, ph1 = iph[j-1]
                i0, ph0 = iph[j-0]
                v10 = a[i1]-a[i0]
                v21 = a[i2]-a[i1]

                r.append(abs(a[i1]))
                av21 = abs(v21)
                L += av21
                l.append(L)
                av10 = abs(v10)
                da = (v21/av21/av21 + v10/av10/av10)/(1/av21+1/av10) 
                phase = log(da)
                phase = (phase.imag+pi)%(2*pi)-pi+pik*2*pi
                dp = phase - w[-1] if w else 0 
                if dp>4:
                    pik-=1
                    phase -= 2*pi
                if dp<-4:
                    pik+=1
                    phase += 2*pi
                w.append(phase)
            if 1/L>rL:
                rP = 1/P
                rL = 1/L
                dwidth = 0
                # print(P,L,'nP0')
                ona = a[iph['ind']]
                oma = array(w)
                lie = array(l)
                rad = array(r)
            elif 1/L==rL:
                dwidth = 1/P - rP 
                # print(P,L,'nP', dwidth)
            else:
                pass
                # print(P,L,'nP-')
        rP += dwidth/2 
        # print('P',1/rP,'L',1/rL)
        # print('P',1/rP,'L',sqrt(width)/rL)
        width /= SD
        Ps = [1/(rP+i*width) for i in range(-2*SD,2*SD+1)]

    return lie,oma,ona,rad,1/rP



Ng = 16
Nl = 160
gk, gw = leggauss(Ng)
gk = (gk+1)/2
lk, lw = linspace(0,2*pi,Nl+1, retstep=True)
    

def draw():
    N = int(exp(slN.val*log(1e1)))
    Nnet = int(exp(slN.val*log(1e1)/4))
    Nbin = int(slN.val*3.321928)+7
    argS = slX.val 
    nu = slM.val
    P0 = slP0.val
    P1 = slP1.val #- slX.val 
    P2 = slP2.val #- slX.val 
    if True: #tries 
        try :
            Funct = eval('lambda s, P0, P1, P2: '+ tblf.text)
            Func = lambda s: Funct(s,P0,P1,P2)
        except:
            # Func = lambda s: cos(s)*(1+P1)+1j*(1-P1)*sin(s)+exp(1j*pi*P2+tan(P0*pi/2))
            
            # Func = lambda s: cos(s)**2*(1+P1)+1j*(1-P1)*sin(s)**2+exp(1j*pi*P2+tan(P0*pi/2))
            
            # Func =  lambda s: (1-P1*P1)*exp(1j*s)/(1+P1*cos(s))+exp(1j*pi*P2*3+tan(P0*pi/2))
            
            Func =  lambda s: exp(1j*(s+P1*pi+tan(P0*pi/2)*cos(s)))+2/(1+P2)

            # Func = lambda s: (cos(P2*pi)*(1+P1**2)+1j*(1-P1**2)*sin(P2*pi)+P1*2*sin(s))*exp(-1j*P2*pi+cos(s)*tan(P0*pi/2))

    net = zeros((Nnet+1,Nnet+1),dtype=int)
    net2 = zeros((Nnet*4+1,Nnet*4+1),dtype=int)
  
    
    s = 0
    c = 0
    cp = angle(Func(0.0))
    pik = 0 
    for k in range(Nl):
        for j in range(Ng):
            s+=log(abs(Func(lk[k]+lw*gk[j])))*gw[j]
            cn = angle(Func(lk[k]+lw*gk[j]))
            if abs(cn-cp)>pi:
                pik += sign(cp-cn)*2*pi
                print(cn-cp,pik)
            cp = cn 
            c += (cn+pik)*gw[j]
    
    # s = s+1j*c
    print(slX.val,slM.val,slP0.val,slP1.val,slP2.val,c*lw/4/pi)
    R = exp(1j*argS-(s+1j*c)*lw/4/pi) 

    try: 
        1/0
        from ethercalc import compute
        a, mabs, mibs = compute.fundamental(S, C, nu, N)
        print('fort')
    except:
        rN =range(N)
        a = zeros((N,),dtype = complex)
        an = 1
        mabs = 1
        mibs = 1
        for n in rN:
            a[n] = an
            if abs(an)> mabs:
                mabs = abs(an)
            if abs(an)< mibs:
                mibs = abs(an)
            an = an *(R*Func(n*nu))

    try:
        M = int(tbAP.text)
        lie,oma,ona,rad,T  = petal(a[::M])
        lp.set_data(ona.real/mabs,ona.imag/mabs)
        L = lie [-1]

    except:
        lp.set_data([0],[0])
        M = 0
        L = 0
        T = 0


    ns = linspace(0,2*pi,1000)
    Fs = R*Func(ns)
    lamax = amax(abs(Fs))
    fl.set_data(Fs.real, Fs.imag)
    fp.set_data(sin(ns), cos(ns))
    sax.set_ylim(-lamax,lamax)
    sax.set_yticks([])
    sax.set_xlim(-lamax,lamax)
    sax.set_xticks([])

    l.set_data(a.real/mabs, a.imag/mabs)
    dim = 1 if M<1.5 else 2
    ax.set_title('outer radius: {:<.2e}; inner radius: {:<.2e}\n dim: {} petals: {}, of length {:<.3f} and period {:<.2f}'.format(mabs,mibs,dim,M,L,T) )

    fig.canvas.draw_idle()
    return a*mabs


    
def update(val):

    draw()

def update_chi(val):
    chi = 0  
    try :
        slX.set_val(2*pi*eval(tbXt.text))
        chi = slX.val
    except:
        pass 
    try :
        slX.set_val(chi+slM.val*eval(tbXm.text))
    except:
        pass
    update(val)


def update_nu(val):
    try :
        slM.set_val(eval(tbM.text))
    except:
        pass
    update(val)

def update_p0(val):
    try :
        slP0.set_val(eval(tbP0.text))
    except:
        pass
    update(val)
 
def update_p1(val):
    try :
        slP1.set_val(eval(tbP1.text))
    except:
        pass
    update(val) 
 
def update_p2(val):
    try :
        slP2.set_val(eval(tbP2.text))
    except:
        pass
    update(val)


def update_func(val):

    update(val)

def save(val):
    print('petal saved')
    a = draw()
    try: 
        M =  int(tbAP.text)
    except:
        M = 1
    cuta = a[::M]
    cuta.tofile('petal.dat')
        

update(1)
slN.on_changed(update)

slP0.on_changed(update)
slP1.on_changed(update)
slP2.on_changed(update)

slX.on_changed(update)
slM.on_changed(update_chi)

tbXt.on_submit(update_chi)
tbXm.on_submit(update_chi)

tbM.on_submit(update_nu)

tbP0.on_submit(update_p0)
tbP1.on_submit(update_p1)
tbP2.on_submit(update_p2)

tblf.on_submit(update_func)

tbAP.on_submit(update)

bSP.on_clicked(save)

show()
