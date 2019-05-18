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
sax.set_position([0.0,.5,.33,.33])
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
axS = axes([0.25, 0.17, 0.65, 0.03], facecolor='gray')
axM = axes([0.25, 0.13, 0.65, 0.03], facecolor='grey')
axC = axes([0.25, 0.09, 0.65, 0.03], facecolor='grey')
axF = axes([0.25, 0.05, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.01, 0.65, 0.03], facecolor='green')
axAP = axes([0.05, 0.065, 0.15, 0.06])
axSrM = axes([0.05, 0.345, 0.15, 0.06])
axSpi = axes([0.05, 0.135, 0.15, 0.06])
axFpi = axes([0.05, 0.205, 0.15, 0.06])
axCpi = axes([0.05, 0.275, 0.15, 0.06])
axSvbutt = axes([0.05, 0.025, 0.15, 0.04]) 
# axSrMbutt = axes([0.05, 0.325, 0.15, 0.04])
# SrMtxt = txt.Text(text = '1.0')
                              # 1.5559549046215868 
slS =  Slider(axS, r'$\chi$', -pi,pi , valinit=0 )
slM =  Slider(axM, r'$\mu$', -4, 4, valinit=0)
slC =  Slider(axC, r'$P1$', -1., 1, valinit=0)
slF =  Slider(axF, r'$P2$', -1, 1, valinit=0)
slN =  Slider(axN, 'N', 3., 5, valinit=3.5)
tbSrM = TextBox(axSrM, r'$\mu$')
bSv = Button(axSvbutt, 'save petal')
tbSpi = TextBox(axSpi, r'${\chi}$')
tbFpi = TextBox(axFpi, r'$P2$')
tbCpi = TextBox(axCpi, r'$P1$')
tbAP = TextBox(axAP, r'$N_{ptls}$')
tblf = TextBox(axlf, r'$\lambda\,s\,:$')
# bSrM = Button(axSrMbutt, 'inverse')
# SrMset = False
# SrMinv = False

# some functions :) 
# Func = lambda s: sin(s)# exp(1j*s)
# Func1 = lambda s: exp(sin(s)*sin(s))*sin(s)/(sin(s)*sin(s)+1)
# Func2 = lambda s: sin(s)/(1-sin(s)*sin(s)+sin(s)*sin(s)*sin(s)*sin(s))
# Func3 = lambda s: sin(s)/(.5+sin(s)*sin(s))
# Func4 = lambda s: (1-2*sin(s)*sin(s))/(2+sin(s)) 
Ellipse = lambda a, s: (a*a+1)*cos(s) + (1-a*a)*1j*sin(s)
Quellipse = lambda a, s: 4j*a**2*cos(s)**2 + (sqrt(4-2*a*a)-sqrt(2)*a)**2*sin(s)**2
# Func = lambda s: Ellipse(1,s)/(1-sin(s)*sin(s)+sin(s)*sin(s)*sin(s)*sin(s)) 
# P1 = 1
# P2 = 0 
# Func = lambda s: P1+exp(1j*P2)*sin(s)


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
    argS = slS.val 
    Mu = exp(slM.val)
    P1 = slC.val #- slS.val 
    P2 = slF.val #- slS.val 
    try :
        Funct = eval('lambda s, P1, P2: '+ tblf.text)
        Func = lambda s: Funct(s,P1,P2)
    except:
        Func = lambda s: P1*2+exp(1j*P2*pi)*sin(s)

    net = zeros((Nnet+1,Nnet+1),dtype=int)
    net2 = zeros((Nnet*4+1,Nnet*4+1),dtype=int)
  
    # z = slC.val**2
    # C = slC.val*2*(cos(F)+1j*sin(F))
    # x = cos(argS)*(1+z)
    # y = (1-z)*sin(argS)
    # S = x*cos(F)-y*sin(F) + 1j*(y*cos(F)+x*sin(F))
    # cf = cos(argS)
    # sf = sin(argS)
    # C = 1-exp(-slC.val*20)
    def logrho1(ro,theta):
        return log( abs( Func(theta)) )

    s = 0
    for k in range(Nl):
        for j in range(Ng):
            s+=logrho1(1., lk[k]+lw*gk[j])*gw[j]
    # s = s*lw/2
    # sv = exp(-s*lw/2/pi)
    S = exp(1j*argS-s*lw/4/pi) 
    try:
        print(2*pi/argS,2*pi/slF.val,slC.val,slM.val,sv,sv*(1-C),sv*C)
    except:
        pass

    try: 
        1/0
        from ethercalc import compute
        a, mabs, mibs = compute.fundamental(S, C, Mu, N)
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
            an = an *(S*Func(n*Mu))

    # for an in a:
    #     # print()
    #     net[ int((an.real/mabs+1)*Nnet/2), int((an.imag/mabs+1)*Nnet/2) ] = 1
    #     net2[ int((an.real/mabs+1)*Nnet*2), int((an.imag/mabs+1)*Nnet*2) ] = 1
    # dim = log(sum(net2)/sum(net))/log(4)
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
    Fs = S*Func(ns)
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

def update_spi(val):
    try :
        slS.set_val(eval(tbSpi.text))
    except:
        pass
    update(val)

 
def update_fpi(val):
    try :
        slF.set_val(float(tbFpi.text))
    except:
        pass
    update(val)

def update_srm(val):
    try :
        slM.set_val(float(tbSrM.text))
    except:
        pass
    update(val)

 
def update_cpi(val):
    try :
        slC.set_val(float(tbCpi.text))
    except:
        pass
    update(val)

def update_func(val):

    update(val)

def save(val):
    print('petal saved')
    a = draw()
    try: 
        M =  eval(tbAP.text)
    except:
        M = 1
    cuta = a[::M]
    cuta.tofile('petal.dat')
        

# print() ro 0.1757649
update(1)
slS.on_changed(update)
slN.on_changed(update)
slM.on_changed(update)
slC.on_changed(update)
slF.on_changed(update)
tbSrM.on_submit(update_srm)
tbSpi.on_submit(update_spi)
tbCpi.on_submit(update_cpi)
tbAP.on_submit(update)
tblf.on_submit(update_func)
bSv.on_clicked(save)
tbFpi.on_submit(update_fpi)
# bSrM.on_clicked(inverse_rel)

show()
