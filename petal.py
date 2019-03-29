# from ethercalc import compute
from matplotlib.pyplot import *
from numpy import *
from random import random as rnd
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt


a = fromfile('petal.dat')
sp = fft.fft(a)
print(a.shape)
freq = fft.fftfreq(a.shape[-1])
plot(freq, sp.real, freq, sp.imag)
show()
