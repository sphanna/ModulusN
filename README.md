# ModulusN
Julia support for Modulus type

The Modulus type is here defined as the residue class of a modulo n where [a] = a + nR := {a +nz | z ∈ R}.  The integer quotient set with respect to Rn, is often denoted by Z/nZ.  Thus, Z/nZ := {[0], [1], . . . , [n − 1]}.

The set of all classes modulo n is a field only in the case that n is prime, i.e. it only guarantees an inverse element if n is prime.  If n is not prime it is possible that a particular integer value has an inverse, but the set is still not considered a field (and may not work with linear algebra operations).

Modulus real values are also supported, however inverse elements and linear algebra functionality will not work with non-integer modulus types.

## Example Usage
Declaration:
Declare a Modulus type with the format: Modulus{N}(a) where N is the modifier and a is the value.  This will create a new type of Modulus{N} and the associated stored value for that type will be mod(a,N).
```julia
julia> using ModulusN
julia> Modulus{4}(10)
Modulus{4}(2)
```

Simple Math:
```julia
julia> a = Modulus{4}(10)
Modulus{4}(2)

julia> b = Modulus{4}(7)
Modulus{4}(3)

julia> a+b
Modulus{4}(1)

julia> a*b
Modulus{4}(2)
```

The module provides support for the following simple operations and conditionals: +,-,*,/,^,abs,abs2,<,>,<=,>=,==

In addition support for finding an inverse element is also implemented:
```julia
julia> a = Modulus{7}(4)
Modulus{7}(4)

julia> inv(a)
Modulus{7}(2)
```

Note that if N is not prime then there is no gaurantee that an inverse exists.
```julia
julia> b = Modulus{6}(4)
Modulus{6}(4)

julia> inv(b)
ERROR: No inverse exists.
```

Automatic promotion and conversion is implemented for real value operations.  Real values are converted to the modulus of the type operated against.
```julia
julia> a = Modulus{5}(9)
Modulus{6}(4)

julia> a + 3
Modulus{5}(2)

julia> a / 3
Modulus{5}(3)
```

Additional functionality (namely, defining a one and zero element) is provided to allow creation of Matrices of Modulus types that work well with the LinearAlgebra package:
```julia
julia> using LinearAlgebra
julia> A = rand(Modulus{5},3,3)
3×3 Array{Modulus{5},2}:
 Modulus{5}(2)  Modulus{5}(3)  Modulus{5}(2)
 Modulus{5}(3)  Modulus{5}(0)  Modulus{5}(3)
 Modulus{5}(4)  Modulus{5}(1)  Modulus{5}(4)

julia> B = rand(Modulus{5},3,3)
3×3 Array{Modulus{5},2}:
 Modulus{5}(0)  Modulus{5}(2)  Modulus{5}(4)
 Modulus{5}(2)  Modulus{5}(4)  Modulus{5}(4)
 Modulus{5}(1)  Modulus{5}(1)  Modulus{5}(2)

julia> A*B
3×3 Array{Modulus{5},2}:
 Modulus{5}(3)  Modulus{5}(3)  Modulus{5}(4)
 Modulus{5}(3)  Modulus{5}(4)  Modulus{5}(3)
 Modulus{5}(1)  Modulus{5}(1)  Modulus{5}(3)

julia> det(B)
Modulus{5}(2)

julia> L,U = lu(B)
LU{Modulus{5},Array{Modulus{5},2}}
L factor:
3×3 Array{Modulus{5},2}:
 Modulus{5}(1)  Modulus{5}(0)  Modulus{5}(0)
 Modulus{5}(3)  Modulus{5}(1)  Modulus{5}(0)
 Modulus{5}(0)  Modulus{5}(3)  Modulus{5}(1)
U factor:
3×3 Array{Modulus{5},2}:
 Modulus{5}(2)  Modulus{5}(4)  Modulus{5}(4)
 Modulus{5}(0)  Modulus{5}(4)  Modulus{5}(0)
 Modulus{5}(0)  Modulus{5}(0)  Modulus{5}(4)
 ```
