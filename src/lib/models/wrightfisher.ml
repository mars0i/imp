(** Wrightfisher:
    Calculate Wright-Fisher transition probabilities with selection and
    create plots of probabilities of population frequencies.

    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

module Mat = Owl.Mat
module Prob = Owl.Stats.Pdf
module L = Batteries.List
module LL = Batteries.LazyList
module U = Utils.Genl

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

let make_init_dist allele_popsize a1count =
  let m = Mat.zeros 1 (allele_popsize + 1) in
  Mat.set m 0 a1count 1.0;
  m

type fitnesses = {w11 : float; w12 : float; w22 : float}

(** 1.59 in Ewens *)
let weight_i {w11; w12; w22} allele_popsize freq  =
  let i, i' = float freq, float (allele_popsize - freq) in
  let a_hom = w11 *. i *. i in
  let het   = w12 *. i *. i' in
  let b_hom = w22 *. i' *. i' in
  (a_hom +. het) /. (a_hom +. 2. *. het +. b_hom)

(** Wright-Fisher transition probability from frequency prev_freq (row index)
    to frequency next_freq (column index). *)
let prob_ij fitns allele_popsize prev_freq next_freq =
  let p = weight_i fitns allele_popsize prev_freq in
  Prob.binomial next_freq p allele_popsize

(** Make a transition matrix from fitnesses *)
let make_tranmat allele_popsize fitns =
  (* prob_ij with an extra ignored argument: *)
  let prob_ij_ fitns allele_popsize prev_freq next_freq _ =
    prob_ij fitns allele_popsize prev_freq next_freq
  in
  let dim = allele_popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (prob_ij_ fitns allele_popsize) m

(** Return pair of old distribution vector (i.e. the second argument) and
    distribution vector (i.e. the product of the two arguments), for use 
    as an element in a Batteries.LazyList *)
let next_dist tranmat dist = 
  (dist, dist *@ tranmat)

(** Returns a LazyList of probability distributions that are the result of 
    a stationary Markov chain with initial distribution (usually with 1 at 
    one entry and zero elsewhere, representing that the population has that 
    initial frequency) and a transition matrix. The initial distribution is 
    included in the list.*)
let make_dists tranmat init_dist =
  LL.from_loop init_dist (next_dist tranmat)

(** Given a list of transition matrices, mats, and a list of
    probability distributions, dists, returns a new list of
    probability distributions produced by multiplying all
    distributions by all matrices. *)
let next_dists tranmats dists =
  (dists, L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                          dists))

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
  LL.from_loop init_dists (next_dists tranmats)

(** Like make_distlists_from_mats, but uses basic parameters to generate the 
    transition matrices and initial distributions that are arguments to 
    make_distslists.  Arg 1, size is the number of alleles, i.e. 2N; arg 2,
    init_freqs, is a list of all initial frequencies for the population 
    (usually there is only one, so the list will have only one element); 
    arg 3, fitn_list is a list of fitness structures. *)
let make_distlists size init_freqs fitn_list =
  let init_dists = L.map (make_init_dist size) init_freqs in
  let tranmats = L.map (make_tranmat size) fitn_list in
  make_distlists_from_mats tranmats init_dists

(** Given a list of float fitness values, which should be in the order
       w11, w12, w22, w11, w12, w22, ...
    eat them in groups of three, using each three to create a
    fitness record and return a list of these records in order. 
    This can be used for commandline processing. *)
let group_fitns fitn_float_list =
  let rec loop l acc =
    match l with
    | [] -> acc
    | w11::w12::w22::tl -> loop tl ({w11=w11; w12=w12; w22=w22}::acc)
    | _ -> raise (Failure "Missing/extra fitness(es)")
  in L.rev (loop fitn_float_list [])
