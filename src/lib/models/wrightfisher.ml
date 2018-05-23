(** Wrightfisher:
    Calculate Wright-Fisher transition probabilities with selection and
    create plots of probabilities of population frequencies.

    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

module L = Batteries.List
module Mat = Owl.Mat
module Stats = Owl.Stats

module TM = Tranmats
module U = Utils.Genl

(** One goal here is to create a "distslist", which is a lazy list of Lists 
    Owl row vector matrices representing probability distributions over 
    possible frequencies of alleles (or other organism types0 in a population. 
    The basic distlist creation functions are in tranmats.ml. *)

let make_init_dist allele_popsize a1count =
  let m = Mat.zeros 1 (allele_popsize + 1) in
  Mat.set m 0 a1count 1.0;
  m

type fitnesses = {w11 : float; w12 : float; w22 : float}

(** Return a new fitnesses record in which w11 and w12 have changed places. *)
let swap_fitns {w11; w12; w22} = {w11=w22; w12=w12; w22=w11}

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
  Stats.binomial_pdf next_freq ~n:allele_popsize ~p:p

(** [make_tranmat allele_popsize fitns] makes a transition matrix for 
    a population with size [allele_popsize/2] from fitnesses [fitns]. 
    Note that the dimension of the square transition matrix will be
    [allele_popsize + 1] since 0 and [allele_popsize] are allowed
    frequencies. *)
let make_tranmat allele_popsize fitns =
  let dim = allele_popsize + 1 in  (* frequencies from zero to N *)
  Mat.init_2d dim dim (prob_ij fitns allele_popsize)

(** Like make_distlists_from_mats, but uses basic parameters to generate the 
    transition matrices and initial distributions that are arguments to 
    make_distslists.  Arg 1, size is the number of alleles, i.e. 2N; arg 2,
    init_freqs, is a list of all initial frequencies for the population 
    (usually there is only one, so the list will have only one element); 
    arg 3, fitn_list is a list of fitness structures. Note that the first
    element will simply be init_freqs. *)
let make_distlists size init_freqs fitn_list =
  let init_dists = L.map (make_init_dist size) init_freqs in
  let tranmats = L.map (make_tranmat size) fitn_list in
  TM.make_distlists_from_mats tranmats init_dists

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
