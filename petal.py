# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from scipy.optimize import curve_fit as crvfit


a = fromfile('petal.dat',dtype = complex)
N = a.shape[-1]
na = a - sum(a)/N
sp = fft.fft(na)[:N//2]
pf = fft.fftfreq(N)[:N//2]
asp = abs(sp)

def profile(x, A, M, S):
    # return A * exp(- (x-M)**2 / S**2)
    return A/(1+(x-M)**2/S**2)

def fp(asp):
    mi = 0  
    mx = 0
    for i in range(len(asp)):
        if asp[i]> mx:
            mx = asp[i]
            mi = i
    casp = asp[mi-5:mi+5]
    cupf = pf[mi-5:mi+5]

    p = crvfit(profile, cupf, casp, [mx, pf[mi], pf[mi+1]-pf[mi-1]] )[0]
    cupft=linspace(cupf[0],cupf[-1],200)
    plot(cupft,profile(cupft,p[0],p[1],p[2]),'b-')
    plot(cupf,casp,'r.')
    plot(pf,sp.real,'-')
    plot(pf,sp.imag,'-')
    show()
    return 1/p[1]




P= fp(asp)
print(P)
# P = 1/0.29772
# print(P)
iph = sort(array([ (i,(1.*i)%P) for i in range(int(3*sqrt(N)))], dtype=[('ind', '<i4'),('phase', '<f8')] ),order='phase')
# print(inds)

ona = a[iph['ind']]


L = 0
w = []
l = []
pik = 0
for j in range(-2,len(iph)):
    l.append(L)
    i, ph = iph[j-1]
    i1, ph1 = iph[j]
    v = a[i1]-a[i]
    av = abs(v)
    L += abs(v)
    phase = log(v/av)
    phase = phase.imag%(2*pi)+pik*2*pi
    dp = phase - w[-1] if w else 0 
    if dp>4:
        pik-=1
        phase -= 2*pi
    if dp<-4:
        pik+=1
        phase += 2*pi
    w.append(phase)

w = array(w)

figure(1)
title("Length = "+str(L))
plot(ona.real,ona.imag,'-')
for j in range(-len(iph)//20,len(iph),len(iph)//10):
    plot(ona[j].real, ona[j].imag,'.')
axis('scaled')

figure(2)
plot(l, w)
for j in range(-len(iph)//20,len(iph),len(iph)//10):
    plot(l[j], w[j],'.')
show()
# 
