# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from scipy.optimize import curve_fit as crvfit


def petal(a):
    N = a.shape[-1]
    na = a - sum(a)/N
    sp = fft.fft(na)[1:N//2]
    pf = fft.fftfreq(N)[1:N//2]
    asp = abs(sp)

    profile = lambda x, A, M, S: A/(1 + (x - M)**2/S**2)
    # profile = lambda x, A, M, S: A*exp(-(x - M)**2/S**2) 
    mix = []
    masp = amax(asp)
    for i in range(len(asp)):
        if asp[i]> masp/2:
            mix.append((i,asp[i]))

    Ps = []
    for  mi,mx in mix:
        casp = asp[mi-5:mi+5]
        cupf = pf[mi-5:mi+5]
        p = crvfit(profile, cupf, casp, [mx, pf[mi], pf[mi+1]-pf[mi-1]] )[0]
        cupft=linspace(cupf[0],cupf[-1],200)
        Ps.append(abs(1/p[1]))

    #     plot(cupft,profile(cupft,p[0],p[1],p[2]),'b-')
    #     plot(cupf,casp,'r.')
    # plot(pf,asp,'-')
    # show()


    # print(Ps)
    rL = 0
    liph = int(sqrt(N)*1.5)
    for P in Ps:
        iph = sort(array([ (i,(1.*i)%P) for i in range(liph)], dtype=[('ind', '<i4'),('phase', '<f8')] ),order='phase')

        L = 0
        pik = 0
        l = []
        w = []
        # r = []
        for j in range(0 ,len(iph)):
            l.append(L)
            i, ph = iph[j-1]
            i1, ph1 = iph[j]
            v = a[i1]-a[i]
            # r.append(abs(a[i1]))
            av = abs(v)
            L += abs(v)
            phase = log(v/av)
            phase = (phase.imag+pi)%(2*pi)-pi+pik*2*pi
            dp = phase - w[-1] if w else 0 
            if dp>4:
                pik-=1
                phase -= 2*pi
            if dp<-4:
                pik+=1
                phase += 2*pi
            w.append(phase)
        if 1>rL*L:
            rL = 1/L
            ona = a[iph['ind']]
            oma = array(w)
            lie = array(l)
            # rad = array(r)
    return lie,oma,ona
lie,oma,ona = petal(fromfile('petal.dat',dtype = complex))
L = lie[-1]

liph = len(lie)

figure(1)
title("Length = "+str(L))
plot(ona.real,ona.imag,'-')
plot([0],[0],'r.', markersize=10.)
for j in range(0,liph,liph//10):
    plot(ona[j].real, ona[j].imag,'.')
axis('scaled')

figure(2)
plot(lie, oma)
for j in range(0,liph,liph//10):
    plot(lie[j], oma[j],'.')

# figure(3)
# plot(lie, rad)
# for j in range(0,liph,liph//10):
#     plot(lie[j], rad[j],'.')

show()
