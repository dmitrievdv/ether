import sympy as sp  
import numpy as np
from scipy.misc import comb
from scipy.special import chebyu
import matplotlib.pyplot as plt

n_fin = 5
int_val_n = np.zeros((n_fin,n_fin+1), dtype = complex)

u, a_symb, b_symb, m_symb, k_symb = sp.symbols('u a_symb b_symb m k')
expr = u**(2*k_symb+2)/((1+u**2)**(m_symb+2)*((a_symb+b_symb)+(a_symb-b_symb)*u**2))

a=0.75
b=1*sp.I

for m in range(n_fin):
	for kk in range(m+1):
		k = m-kk
		expr_mk = expr.subs(m_symb, m).subs(k_symb, k)
		integral = sp.integrate(expr_mk, (u, 0, sp.oo))
		# integral_ab = integral.subs(a, 1.25).subs(b, 1)
		int_val_n[m,k] = sp.N(integral, subs = {a_symb:a, b_symb:b}) 
		print(k)
	print('m = ', m, ' : ', int_val_n[m,:m+1])

for n1 in range(n_fin):
	n = n1+1
	u_n1 = chebyu(n-1)

	sum_k = np.zeros(n, dtype = complex)


	for m in range(n):
		for kk in range(m+1):
			k = m-kk
			value = int_val_n[m,k]
			sum_k[m] += value*(-1)**k*comb(m,k)
			
	# print(sum_k)
	# print(np.flip(sum_k)*u_n1)
	print('n = ', n, ': f_n = ', np.sum(np.flip(sum_k)*u_n1)/n*16/np.pi*complex(b))



# expr_ab = expr.subs(a, 1.25)
# expr_ab = expr_ab.subs(b, 1)
# expr_mk = expr.subs(m, 3)
# expr_mk = expr_mk.subs(k, 3)

# integral = sp.integrate(expr_mk, (u, 0, sp.oo))
# integral_ab = integral.subs(a, 1.25).subs(b, 1)
# value = integral.evalf(subs = {a:1.25, b:1}) 
# print(value)