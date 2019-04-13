# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from scipy.optimize import curve_fit as crvfit
from numpy.polynomial.legendre import leggauss

fig, ax = subplots()
ax.set_title( 'scale: {:<.2e}'.format(1.0) )
subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'r.',markersize=1.) 
lp, = plot([0], [0], 'b-',markersize=2.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axS = axes([0.25, 0.17, 0.65, 0.03], facecolor='gray')
axM = axes([0.25, 0.13, 0.65, 0.03], facecolor='grey')
axC = axes([0.25, 0.09, 0.65, 0.03], facecolor='grey')
axF = axes([0.25, 0.05, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.01, 0.65, 0.03], facecolor='green')
axSrM = axes([0.05, 0.065, 0.15, 0.06])
axSpi = axes([0.05, 0.135, 0.15, 0.06])
axFpi = axes([0.05, 0.205, 0.15, 0.06])
axCpi = axes([0.05, 0.275, 0.15, 0.06])
axAP = axes([0.05, 0.345, 0.15, 0.06])
axlf = axes([0.05, 0.445, 0.2, 0.06])
axSvbutt = axes([0.05, 0.025, 0.15, 0.04]) 
# axSrMbutt = axes([0.05, 0.325, 0.15, 0.04])
# SrMtxt = txt.Text(text = '1.0')
                              # 1.5559549046215868 
slS =  Slider(axS, r'$\chi$', -pi,pi , valinit=1.5559548642369434 )
slM =  Slider(axM, r'$\mu$', -4, 4, valinit=1)
slC =  Slider(axC, r'$\rho_t$', 0., 1, valinit=.7)
slF =  Slider(axF, r'$\varphi_t$', -pi, pi, valinit=pi/2)
slN =  Slider(axN, 'N', 3., 5, valinit=3.5)
tbSrM = TextBox(axSrM, r'$\mu$')
bSv = Button(axSvbutt, 'save petal')
tbSpi = TextBox(axSpi, r'$\frac{2\pi}{\chi}$')
tbFpi = TextBox(axFpi, r'$\frac{2\pi}{\varphi_t}$')
tbCpi = TextBox(axCpi, r'$\rho_t$')
tbAP = TextBox(axAP, r'$N_{ptls}$')
tblf = TextBox(axlf, r'$\lambda\,s\,:$')
# bSrM = Button(axSrMbutt, 'inverse')
# SrMset = False
# SrMinv = False

# some functions :) 
# Func = lambda s: exp(s*s)*s/(s*s+1)
# Func = lambda s: s/(1-s*s+s*s*s*s)
# Func = lambda s: s/(.5+s*s)
# Func = lambda s: (1-2*s*s)/(2+s) 
Func = lambda s: s

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
        if asp[i]> masp/3:
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

    return lie,oma,ona,rad



Ng = 16
Nl = 16
gk, gw = leggauss(Ng)
gk = (gk+1)/2
lk, lw = linspace(0,2*pi,Nl+1,retstep=True)
    

def draw():
    N = int(exp(slN.val*log(1e1)))
    Nbin = int(slN.val*3.321928)+7
    F = slF.val - slS.val 
    M = slM.val
    argS = slS.val 
    # tbSpi.set_text(slS.val/pi)

    # z = slC.val**2
    # C = slC.val*2*(cos(F)+1j*sin(F))
    # x = cos(argS)*(1+z)
    # y = (1-z)*sin(argS)
    # S = x*cos(F)-y*sin(F) + 1j*(y*cos(F)+x*sin(F))
    sv = 1e3
    dsv = sv/2
    cf = cos(argS)
    sf = sin(argS)
    C = slC.val*2
    def logrho1(ro,theta):
        return log((cf*ro+C*Func(sin(theta)))**2+sf*ro*sf*ro)

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
    
    

    S = sv*exp(1j*F)*exp(1j*argS)
    C = C*exp(1j*F)
    print(argS,F,sv)

    try: 
        1/0
        from ethercalc import compute
        a, mabs, mibs = compute.fundamental(S, C, M, N)
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
            an = an *(S + C*Func(sin(n*M)))
    rN = slF.val/(2*pi)
    MM = 99
    for M in range(1,MM+1):
        if abs(round(rN*M)-rN*M)<1/MM/MM:
            break
    try:
        M = int(tbAP.text)
        lie,oma,ona,rad  = petal(a[::M])
        lp.set_data(ona.real,ona.imag)
        L = lie [-1]

    except:
        lp.set_data([0],[0])
        M = 0
        L = 0

    # a = a/mabs
    l.set_data(a.real, a.imag)
    # print(slS.val,mibs,mabs     )
    # lp.set_data(a.real[::4], a.imag[::4])
    ax.set_title('outer radius: {:<.2e}; inner radius: {:<.2e}\n petals: {}, of length {:<.3f}'.format(mabs,mibs,M,L) )

    ax.axis('scaled')
    ax.set_xlim(-mabs,mabs)
    ax.set_ylim(-mabs,mabs)
    # ax.set_xlim(-1,1)
    # ax.set_ylim(-1,1)
    fig.canvas.draw_idle()
    return a*mabs

# def update_n(val):
#     if(SrMset):
#         update_rel(val)
#     else:
#         update(val)

    
def update(val):
    # global SrMset
    # if(SrMset and not SrMinv):
    #     slM.disconnect(0)
    #     slM.cnt = 0
    #     slM.on_changed(update)
    # elif(SrMset):
    #     slS.disconnect(0)
    #     slS.cnt = 0
    #     slS.on_changed(update)
    draw()
    # SrMset = False

# def update_txt(val):
#     global SrMset
#     if(not SrMinv):
#         slM.disconnect(0)       
#         slM.cnt = 0
#         cid = slM.on_changed(update_rel)
#     else:
#         slS.disconnect(0)
#         slS.cnt = 0
#         slS.on_changed(update_rel)
#     update_rel(val)
#     SrMset = True

# def update_rel(val):
#     global SrMset
#     if(not SrMinv):
#         M = slM.val
#         argS = float(tbSrM.text)*M
#         SrMset = False
#         slS.set_val(argS)
#         SrMset = True
#     else:
#         argS = slS.val
#         M = float(tbSrM.text)*argS
#         SrMset = False
#         slM.set_val(M)
#         SrMset = True
#     draw()
    
# def inverse_rel(val):
#     global SrMinv
#     SrMinv = not SrMinv
#     tbSrM.label.set_visible(False)
#     tbSrM.text_disp.set_visible(False)
#     initial = tbSrM.text
#     if(not SrMinv):
#         tbSrM.__init__(axSrM, r'$\varphi/\mu$', initial = initial)
#         tbSrM.on_submit(update_txt)
#         slS.disconnect(0)
#         slS.cnt = 0
#         slS.on_changed(update)
#     else:
#         tbSrM.__init__(axSrM, r'$\mu/\varphi$', initial = initial)
#         tbSrM.on_submit(update_txt)
#         slM.disconnect(0)
#         slM.cnt = 0
#         slM.on_changed(update)
#     try: 
#         update_txt(val) 
#     except:
#         pass 

def update_spi(val):
    slS.set_val(2*pi/float(tbSpi.text))
    update(val)

 
def update_fpi(val):
    slF.set_val(2*pi/float(tbFpi.text))
    update(val)

def update_srm(val):
    slM.set_val(float(tbSrM.text))
    update(val)

 
def update_cpi(val):
    slC.set_val(float(tbCpi.text))
    update(val)

def update_func(val):
    global Func
    Func = eval('lambda s : '+ tblf.text)
    update(val)

def save(val):
    print('ololol')
    a = draw()
    M = int(tbAP.text)
    cuta = a[::M]
    cuta.tofile('petal.dat')
        

# print()
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
