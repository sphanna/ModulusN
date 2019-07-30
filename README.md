# ModulusN
Julia support for Modulus type

The Modulus type is here defined as the residue class of a modulo n where
[a] = a + nZ := {a +nz | z ∈ Z}
the quotient set with respect to Rn, is often denoted by Z/nZ.
Thus, Z/nZ := {[0], [1], . . . , [n − 1]}.

The set of all classes modulo n is a field only in the case that n is prime.
i.e. it only guarantees an inverse element if n is prime.  If n is not prime
it is possible that a particular value, a,has an inverse (that is a*a^-1 = 1),
but the set is still not considered a field.
