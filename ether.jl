module Ether

using Polynomials

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
# function get_next(an :: Number, n :: Int)
# 	return an + sin(n)
# end


end