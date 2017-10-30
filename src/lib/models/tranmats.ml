(** Tranmats: Lists, Lazy Lists, and other structures containing
    transition matrices or vectors produced by them.  Agnostic about
    the source of the matrices.  For example, there's nothing here about
    fitnesses. *)

module Mat = Owl.Mat
module L = Batteries.List
module LL = Batteries.LazyList

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

let always_true x = true

type dists_at_t = {t : int; dists : Mat.mat list}

(** Return pair of old distribution vector (i.e. the second argument [dist]
    of this function) and a distribution vector (i.e. the product of the two 
    arguments).  For use with [Batteries.LazyList.from_loop] *)
let deprecated_next_dist tranmat dist = 
  (dist, dist *@ tranmat)

(** Returns a LazyList of probability distributions that are the result of 
    a stationary Markov chain with initial distribution (usually with 1 at 
    one entry and zero elsewhere, representing that the population has that 
    initial frequency) and a transition matrix. The initial distribution is 
    included in the list.*)
let deprecated_make_dists tranmat init_dist =
  LL.from_loop init_dist (deprecated_next_dist tranmat)

(** Given a list of transition matrices, mats, and a list of
    probability distributions, dists, returns a new list of
    probability distributions produced by multiplying all
    distributions by all matrices. *)
let old_next_dists tranmats t_dists =
  let {t; dists} = t_dists in
  (dists,
   {t     = t + 1;
    dists = L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                            dists)})

let new_next_dists tranmats t_dists =
  let {t; dists} = t_dists in
  {t     = t + 1;
   dists = L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                           dists)}

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
       LL.at distlists 1
    will produce 2**1 = 2 dists.  Or with more initial distributions, the
    number of dists at n is (length init_dists) * (length tranmats)**n . *)

let make_distlists_from_mats tranmats init_dists =
  LL.seq init_dists (next_dists tranmats) always_true

let new_make_distlists_from_mats tranmats init_dists =
  LL.seq init_dists (new_next_dists tranmats) always_true

let old_make_distlists_from_mats tranmats init_dists =
  LL.from_loop init_dists (old_next_dists tranmats)

