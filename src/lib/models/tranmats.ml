(** Tranmats: Lists, Lazy Lists, and other structures containing
    transition matrices or vectors produced by them.  Agnostic about
    the source of the matrices.  For example, there's nothing here about
    fitnesses. *)

module Mat = Owl.Mat
module L = Batteries.List
module T = Tdists

(** The general goal here is to create a "distlist", which is a 
    LazyList of Lists Owl row vector matrices representing probability
    distributions over possible frequencies of alleles (or other organism
    types0 in a population. *)

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

let always_true _ = true

(** Given a list of transition matrices, mats, and a list of
    probability distributions, dists, returns a new list of
    probability distributions produced by multiplying all
    distributions by all matrices. *)
let next_dists tranmats dists =
   L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                   dists)

(** Given a list of transition matrices and a list of initial distributions
    (often one distribution with all weight on one frequency),
    return a lazy list of lists of probability distributions.
    Note this function does not not drop the first element. That way, the 
    number of dists in the nth distlist = (length tranmats)**n for one 
    initial distribution.  e.g. with
    two transition matrices and one initial distribution,
       T.at distlists 1
    will produce 2**1 = 2 dists.  Or with more initial distributions, the
    number of dists at n is (length init_dists) * (length tranmats)**n . *)
let make_distlists_from_mats tranmats init_dists =
  T.seq init_dists (next_dists tranmats) always_true
