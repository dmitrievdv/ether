import numpy as np 
import scipy.integrate as integrate
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.text as txt
from matplotlib.pyplot import *

fig, ax = subplots()

subplots_adjust(left=0.25, bottom=0.25) 
l, = plot([1], [0], 'b-',markersize=1.)
# ll, = plot([1], [0], 'b-',markersize=1.)
ax.axis('scaled')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

axaC = axes([0.25, 0.15, 0.65, 0.03], facecolor='gray')
axrC = axes([0.25, 0.11, 0.65, 0.03], facecolor='grey')
axN = axes([0.25, 0.07, 0.65, 0.03], facecolor='grey')
axtet = axes([0.25, 0.03, 0.65, 0.03], facecolor='green')
# SrMtxt = txt.Text-np.pi, np.pi = '1.0')

slrt =  Slider(axaC, r'$\rho_t$', 0, 1, valinit=0.5)
slpt =  Slider(axrC, r'$\varphi_t$', -np.pi, np.pi, valinit=0)
slps =  Slider(axN, r'$\varphi_s$', -np.pi, np.pi, valinit=0)
sltt =  Slider(axtet, r'$\theta_k/\pi$', 0, 20, valinit=2)

def draw(rhot,phit,phis,tetk):         
	rhos = (1/(np.cos(phis - phit)**2/(1+rhot**2)**2+np.sin(phis - phit)**2/(1-rhot**2)**2))**0.5
	S = rhos*np.exp(1.0j*phis)
	C = 2*rhot*np.exp(1.0j*phit)
	# F = 1
	# M = 1
	# argS = np.pi/2 -F
	# z = slC.val**2
	# C = slC.val*2*(cos(F)+1j*sin(F))
	# x = cos(argS)*(1+z)
	# y = (1-z)*sin(argS)
	# S = x*cos(F)-y*sin(F) + 1j*(y*cos(F)+x*sin(F))
	f = lambda tet: np.log(S+C*np.sin(tet))
	# thetk = 8*np.pi
	n = 5000
	m = 0
	ints = np.zeros(n, dtype='complex')
	thets = np.linspace(0, tetk, n, endpoint = True)
	for i,thet in enumerate(thets):
		m = f(thet)*tetk/n + m
		ints[i] = np.exp(m)
	# ints = np.exp(ints)
	x = ints.real
	y = ints.imag
	ro = np.sqrt(x**2 + y**2)
	rom = np.amax(ro)
	ax.set_title( 'Integral, scale: ' + str(rom) )
	l.set_data(x/rom,y/rom)
	fig.canvas.draw_idle()

def update(val):
	rhot = slrt.val
	phit = slpt.val
	phis = slps.val
	tetk = np.pi*sltt.val
	# N = int(exp(slN.val*log(1e1)))
	draw(rhot,phit,phis, tetk)

update(0)
slpt.on_changed(update)
slps.on_changed(update)
slrt.on_changed(update)
sltt.on_changed(update)
show()