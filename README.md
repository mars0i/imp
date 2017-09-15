imp
===

Imprecise probability and population genetics in [OCaml](http://ocaml.org/). Experimental
work in progress.

Some of this will be reorganized in the future.  Some of it was
written when I was learning OCaml.  (I am still learning OCaml.)

### What's here

#### src/lib/models

#### wrightfisher.ml: simultaneous Wright-Fisher models, PDF plot generation:

For modeling finite sets of transition probabilities generated by
Wright-Fisher models.

#### setchains.ml: implementations of algorithms in chapter 2 of *Markov Set-Chains*, by Darald J. Hartfiel, Springer 1998:

These model continuous sets of stochastic transition matrices that fall
within a matrix "interval", i.e. all stochastic matrices such that each
element `x` is s.t. `l <= x <= h`, where `l` and `h` are corresponding
elements of a low matrix `L` and a high matrix `H`.
You can also start from an initial probability "interval", i.e. all
stochastic vectors such that each element `x` is s.t. `l <= x <= h`,
where `l` and `h` are correspoding elements of low and high vectors.
Note that the low and high vectors are not typically stochastic vectors,
nor the low and high matrices typically transition matrices.
(Stochastic vectors are row vectors and it's each row of a matrix that
sums to 1, so multiplication of a vector and a matrix typically
happens with the vector on the left.  I use the
[Owl](https://github.com/ryanrhymes/owl) library; both vectors
and matrices are Owl matrices.)

#### src/lib/utils

#### genl.ml: general-purpose utilities    
#### prob.ml: probability-related utilities


#### src/misc
miscellaneous code, potentially useful now or in the future,
including:

#### setchain_egs.ml: definitions for specific examples in Hartfiel (see above)


#### src/bin

#### wrightfisherPDFs.ml: source for command line executable that generates simultaneous Wright-Fisher model PDFs

#### setchaintest.ml: miscellaneous tests for setchains.ml
