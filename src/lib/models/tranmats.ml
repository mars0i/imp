(** Tranmats: Lists, Lazy Lists, and other structures containing
    transition matrices or vectors produced by them.  Agnostic about
    the source of the matrices.  For example, there's nothing here about
    fitnesses. *)

module Mat = Owl.Mat
module L = Batteries.List
module LL = Batteries.LazyList

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
       LL.at distlists 1
    will produce 2**1 = 2 dists.  Or with more initial distributions, the
    number of dists at n is (length init_dists) * (length tranmats)**n . *)
let make_distlists_from_mats tranmats init_dists =
  LL.seq init_dists (next_dists tranmats) always_true

type tdists = {t : int ; dists : Mat.mat list}
(** accessor, constructor functions: *)
let t tds = tds.t
let dists tds = tds.dists
let make_tdists t dists = {t ; dists}

(* Transform lazy lists of dists to/from lazy lists of tdists: *)
let ints_from n = LL.seq n ((+) 1) always_true

let add_ts ?(first_tick=0) dists_llist =
  LL.map2 make_tdists (ints_from first_tick) dists_llist

let remove_ts tdists_llist = LL.map dists tdists_llist

(** [tdists_sublist start_t finish_t tdists_llist] returns a lazy list that's
    a sublist of [tdists_llist], from the element with [t]=[start_t] to the 
    element with [t]=[finish_t], inclusive.  Note that if the list is infinite
    and there is no element with [t]=[start_t] or with [t]=[finish_t], the 
    function will run forever, or until the system is overloaded. *)
let tdists_sublist start_t finish_t tdists_llist =
  LL.take_while (fun tds -> tds.t < finish_t)
                (LL.drop_while (fun tds -> tds.t < start_t)
		               tdists_llist)

(* Alias for [tdists_sublist *)
let subtdsl = tdists_sublist 
