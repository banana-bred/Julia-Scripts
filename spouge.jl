## Generate Γ expansion coefficients for a given numerical precision

# -- set precision
setprec(bits :: Integer)= (setprecision(BigFloat, bits); nothing)

# -- Spouge error bound
spouge_bound(a :: Integer) = 1/sqrt(2BigFloat(pi)) * a^(BigFloat(0.5) - a) * exp(BigFloat(a))

"""
    spouge_coeffs(a :: Int, bits = 256)

Compute Spouge coefficients c0..c(a-1) at given precision (bits)
c0 = sqrt(2π)
ck = (-1)^{k-1}/(k-1)! * exp(a-k)*(a-k)^(k-1/2)/ sqrt(2π), k=1..a-1
"""
function spouge_coeffs(a :: Int, bits :: Int=256)
  @assert a ≥ 3
  setprec(bits)
  c = Vector{BigFloat}(undef, a)
  twopi = 2BigFloat(pi)
  sqrt2pi= sqrt(twopi)
  c[begin + 0]= sqrt(twopi) # c0
  for k in 1:(a-1)
    kk = BigFloat(k)
    ak = BigFloat(a) - kk
    term = (-one(BigFloat))^(k-1) / factorial(big(k-1))
    term *= exp(ak)
    term *= ak^(kk - BigFloat(0.5))
    term /= sqrt2pi
    c[begin + k] = term
  end
  return c
end

" Pick the minimal a for at least `digits` decimal digits"
function a4digits(digits :: Integer)
  targ = big(10.0)^(-digits)
  a =3
  while spouge_bound(a) > targ
    a+=1
  end
  return a
end
