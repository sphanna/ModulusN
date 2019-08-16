#Code written by Scott Hanna 6/27/2019
#Last Update: 8/14/2019
#Contact: shanna7@jhu.edu
#
#The Modulus type is here defined as the residue class of a modulo n where
#[a] = a + nZ := {a +nz | z ∈ Z}
#the quotient set with respect to Rn, is often denoted by Z/nZ.
#Thus, Z/nZ := {[0], [1], . . . , [n − 1]}.
#
#The set of all classes modulo n is a field only in the case that n is prime.
#i.e. it only guarantees an inverse element if n is prime.  If n is not prime
#it is possible that a particular value, a,has an inverse (that is a*a^-1 = 1),
#but the set is still not considered a field.
module ModulusN

using LinearAlgebra, Random
import Base: convert,promote_rule,+,-,*,/,zero,one,abs,abs2,^,<,<=,>,>=,==,isequal,inv,isnan
import LinearAlgebra: conj
import Random: rand

export getVal,getMod,Modulus,hasInverse

struct Modulus{N} <: Number
    val::Real
    function Modulus{N}(val) where N
        @assert N isa Real
        @assert N != 0
        new{N}(mod(val,N))
    end
    Modulus{N}(x::Modulus{N}) where N = Modulus{N}(x.val)
end

convert(::Modulus{N}, x::M) where {N,M<:Real} = Modulus{N}(x)
promote_rule(::Type{Modulus{N}}, ::Type{M}) where {N,M<:Real} = Modulus{N}
promote_rule(::Type{M}, ::Type{Modulus{N}}) where {N,M<:Real} = Modulus{N}

getVal(x::Modulus{N}) where N = return x.val
getMod(x::Modulus{N}) where N = return N

+(x::Modulus{N},y::Modulus{N}) where N = Modulus{N}(x.val + y.val)
-(x::Modulus{N},y::Modulus{N}) where N = Modulus{N}(x.val - y.val)
*(x::Modulus{N},y::Modulus{N}) where N = Modulus{N}(x.val * y.val)
/(x::Modulus{N},y::Modulus{N}) where N = inv(y) * x
-(x::Modulus{N}) where N = Modulus{N}(-x.val)
zero(x::Modulus{N}) where N = Modulus{N}(0)
one(x::Modulus{N}) where N = Modulus{N}(1)
abs(x::Modulus{N}) where N = Modulus{N}(abs(x.val))
^(x::Modulus{N}, y::Int) where N = Modulus{N}(x.val^y)
abs2(x::Modulus{N}) where N = abs(x)^2

<(x::Modulus{N}, y::Modulus{N}) where N = x.val < y.val
<=(x::Modulus{N}, y::Modulus{N}) where N = x.val <= y.val
>(x::Modulus{N}, y::Modulus{N}) where N = x.val > y.val
>=(x::Modulus{N}, y::Modulus{N}) where N = x.val >= y.val
==(x::Modulus{N}, y::Modulus{M}) where {N,M} = x.val == y.val && N==M
isequal(x::Modulus{N}, y::Modulus{M}) where {N,M} = x.val == y.val && N==M

isnan(x::Modulus{N}) where N = isnan(x.val)
conj(x::Modulus{N}) where N = abs(x)

function inv(x::Modulus{N}) where N
    divs,v = gcdx(x.val,N)
    if divs != 1 || N == 1
        error("No Inverse Exists")
    end
    return Modulus{N}(v)
end

function hasInverse(x::Modulus{N}) where N
    return gcd(x.val,N)==1 && N != 1
end

rand(rng::AbstractRNG, ::Random.SamplerType{Modulus{N}}) where N = Modulus{N}(rand(rng,Int))

end #end module ModulusN
