module Ether

using Polynomials
using Plots

function chebsecondkind(n::Int)
	U₀ = Poly([1])
	U₁ = Poly([0,2])
	x = Poly([0,1])
	if n == 0
		return U₀
	elseif n == 1
		return U₁
	end
	Uₙ₋₁ = U₁
	Uₙ₋₂ = U₀
	Uₙ = Uₙ₋₁
	for i=2:n
	    Uₙ = 2*x*Uₙ₋₁ - Uₙ₋₂
	    Uₙ₋₂ = Uₙ₋₁
	    Uₙ₋₁ = Uₙ
	end
	return Uₙ
end

function mppowerreduction(m::Int, p::Int, Ap, Bp, Cp, c, d)
	C = -1/(2*(c-d)*(p+1))*Cp
	A = c*(Ap-Bp)*(m-1) 
	B = d*(Ap-Bp)*(m-1) - 2*(Bp*c-Ap*d)*(p+1)
	m -= 2
	p += 1
	return (m, p, A, B, C)
end

function ppowerreduction(p::Int, Ap, Bp, Cp, c, d)
	C = 1/(2*(c-d)*(p+1))*Cp
	A = c*(Ap-Bp) + 2*Ap*(c-d)*(p+1)
	B = d*(Ap-Bp)*(1 + 2*(p+1))
	p += 1
	return (p, A, B, C)
end

function mkintegral(m::Int, k::Int, a, b)
	A = 1; B = 0
	c = a+b; d = a-b
	C = 1
	n = 2*k+2
	p = -(m+2)
	# print(n, " ",  p, " ", A, " ", B, " ", C, "\n")
	# print(c, " ",  d, "\n")
	while n > 0
	    n, p, A, B, C = mppowerreduction(n, p, A, B, C, c, d)
	    # print(n, p, A, B, C, "\n")
	end
	while p < -1
		p, A, B, C = ppowerreduction(p, A, B, C, c, d)
		# print(n, p, A, B, C, "\n")
	end
	integral = π*(A/√c + B/√d)/(2*(√c + √d))
	return C*integral
end 

function ksum(m::Int, a, b)
	s = 0
	for k=0:m
		# print(m, k, a, b, " ")
		# print(mkintegral(m, k, a, b), "\n")
		s += mkintegral(m, k, a, b)*(-1)^k*binomial(m,k)
	end
	return s
end

function fₙ(n::Int, a, b)
	s = 0
	Uₙ₋₁ = chebsecondkind(n-1)
	for m=0:(n-1)
		s += ksum(m, a, b)*Uₙ₋₁[m]
	end
	s = s*16*b/n/π
	return s
end

function expfₙ(n::Int, a, b)
	if n == 0
		return 0
	elseif n > 0
		return 1/2*fₙ(n, a, b)
	elseif n < 0 
		return 1/2*fₙ(-n, a, b)
	end
end


function calc_f(m::Int, a, b)
	function f(x)
		s = 0
		for n=-m:m
			s += expfₙ(n, a, b)*exp(-1.0im*n*x)
		end
		return s
	end 
	return f
end

function calc_g(m::Int, a, b)
	function g(x)
		s = 0
		for n=1:m
			# print(expfₙ(n, a, b)*(exp(-1.0im*n*x) - sin(n*x + n/2.0)/sin(n/2.0))*exp(-1.0im*n*x), "\n")
			s += expfₙ(n, a, b)*(exp(-1.0im*n*x) - sin(n*x + n/2.0)/sin(n/2.0))*exp(-1.0im*n*x)
			s += expfₙ(n, a, b)*(exp(1.0im*n*x) - sin(n*x + n/2.0)/sin(n/2.0))*exp(1.0im*n*x)
			# print(n, " ", expfₙ(n, a, b), " ", s, "\n")
		end
		return s
	end 
	return g
end

function calc_petal(t, φ; n = 10)
	areal = cos(φ) * (1+t^2)
	aimag = sin(φ) * (1-t^2)
	a = areal + aimag*1.0im
	# print(a)
	return calc_g(n, a, 2*t)
end

function plot_petal(t, φ; nx = 400)
	function petal(n); calc_petal(t, φ; n = n); end
	x = [0:π/nx:π;]
	pyplot(size = (700, 500))
	plot(real.(exp.(petal(10).(x))), imag.(exp.(petal(10).(x))), aspect_ratio="equal", label = "10")
	plot!(real.(exp.(petal(20).(x))), imag.(exp.(petal(20).(x))), aspect_ratio="equal", label = "20")
	plot!(real.(exp.(petal(30).(x))), imag.(exp.(petal(30).(x))), aspect_ratio="equal", label = "30")
end

# function get_next(an :: Number, n :: Int)
# 	return an + sin(n)
# end


end