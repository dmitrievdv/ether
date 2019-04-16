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
        # Ps = [rP*(1+i*width/rP/20) for i in range(-30*SD,30*SD+1)]

    return lie,oma,ona,rad


lie,oma,ona, rad = petal(fromfile('petal.dat',dtype = complex))
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
yr = range(int(amin(oma)*2/pi),int(amax(oma)*2/pi)+1)
yl = []
for yt in yr:
    if yt == 0 :
        yl.append('0')
    elif yt%2:
        yl.append(str(yt)+r'$\frac{\pi}{2}$')
    else:
        yl.append(str(yt//2)+r'$\pi$')

yticks(array(yr)*pi/2,array(yl))


for j in range(0,liph,liph//10):
    plot(lie[j], oma[j],'.')

# figure(3)
# plot(lie, rad)
# for j in range(0,liph,liph//10):
#     plot(lie[j], rad[j],'.')

show()
