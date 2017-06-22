
(** Wrightfisher:
    Calculate Wright-Fisher transition probabilities with selection and
    create plots of probabilities of population frequencies.

    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

module Mat = Owl.Mat
module Math = Owl.Maths (* note British->US translation *)
module Pl = Owl.Plot
module L = Batteries.List
module LL = Batteries.LazyList

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)


(** Return a lazy list that's a sublist of the argument, from element start 
    (zero-based) to element finish, inclusive. *)
let sub_lazy_list start finish ll =
  LL.take (finish - start + 1) (LL.drop start ll)

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
  let wt = weight_i fitns allele_popsize prev_freq in
  let other_wt = 1. -. wt in
  let j = float next_freq in
  let j' = float (allele_popsize - next_freq) in
  let comb = Math.combination_float allele_popsize next_freq in
  comb  *.  wt**j  *.  other_wt**j'

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

(** Convenience function: Takes elements from start to finish, inclusive, 
    from a LazyList, and convert the result to a List. *)
let take_to_list start finish ll = 
  LL.to_list (sub_lazy_list start finish ll)

(** Return second dimension of a matrix or vector. *)
let length m = snd (Mat.shape m)

(* Make a series of n 2D plot pdfs from dists using basename. *)
let make_2D_pdfs basename dists n =
  let dist_length = length (LL.at dists 0) in
  let xs = Mat.sequential 1 dist_length in (* vector of x-axis indices *)
  let make_pdf i dist =
    let filename = basename ^ (Printf.sprintf "%03d" i) ^ ".pdf" in
    let h = Pl.create filename in 
    Pl.set_yrange h 0.0 0.25;
    Pl.scatter ~h xs dist; 
    Pl.output h
  in LL.iteri make_pdf (LL.take n dists)

(** Return a triple containing x-coord, y-coord, and z-coord matrices.
    dist_list is a list of row vectors representing prob dists that will
    be concatenated into z coords.  Note that the x and y coord matrices
    will have the same shape, which will be transposed/rotated wrt the z
    coord matrix that results. That's what Owl.Plot.{mesh,surf} need. *)
let make_coords dist_list =
  let (_, width) = Mat.shape (L.hd dist_list) in
  let height = L.length dist_list in
  let xs = Mat.repeat ~axis:0 (Mat.sequential 1 height) width in
  let ys = Mat.repeat ~axis:1 (Mat.sequential width 1) height in
  let zs = L.reduce Mat.concat_vertical dist_list in
  (xs, ys, zs)


(** Given a list of transition matrices, mats, and a list of
    probability distributions, dists, returns a new list of
    probability distributions produced by multiplying all
    distributions by all matrices. *)
let next_dists tranmats dists =
  (dists, L.concat (L.map (fun dist -> L.map (Mat.dot dist) tranmats)
                          dists))

(** Given a list of transition matrices and a list of initial distributions
    (often one distribution with all weight on one frequency) 
    Note this function does not not drop the first element. That way, the 
    number of dists in the nth distlist = (length tranmats)**n for one 
    initial distribution.  e.g. with
    two transition matrices and one initial distribution,
       LL.at distlists 1
    will produce 2**1 = 2 dists.  Or with more initial distributions, the
    number of dists at n is (length init_dists) * (length tranmats)**n . *)
let make_distlists_from_mats tranmats init_dists =
  LL.from_loop init_dists (next_dists tranmats)

(** Like make_distlists_from_mats, but uses basic parameters to generate the transition
    matrices and initial distributions that are arguments to make_distslists.
    Arg 1, size is the number of alleles, i.e. 2N; arg 2, init_freqs, is a 
    list of all initial frequencies for the population (usually there is only
    one, so the list will have only one element); arg 3, fitn_list is a list
    of fitness structures. *)
let make_distlists size init_freqs fitn_list =
  let init_dists = L.map (make_init_dist size) init_freqs in
  let tranmats = L.map (make_tranmat size) fitn_list in
  make_distlists_from_mats tranmats init_dists

(* Given a list of probability distribution vectors, tries to sort them
   so that similar lists are close in the order.  
   Uses Utils.difference_compare. *)
let sort_dists dists = L.sort Utils.difference_compare dists

(** Make a series of n 3D plot pdfs from distlists using basename.
    Example:
    let distlists = make_distlists 500 [200] 
                  [{w11=1.0; w12=0.8; w22=0.7}; {w11=1.0; w12=0.3; w22=0.7}];;
    make_3D_pdfs "distsN=500init=200w11=1w22=0.7w12=0.8or0.3gen" distlists 9;;
 *)
let make_3D_pdfs ?(altitude=50.) ?(azimuth=75.) basename distlists start_gen last_gen =
  let make_pdf i dists =  (* i = t-1; dists = prob dists at t *)
    let gen = i + start_gen in
    let filename = basename ^ (Printf.sprintf "%03d" gen) ^ ".pdf" in 
    let xs, ys, zs = make_coords (sort_dists dists) in
    let h = Pl.create filename in
      Pl.set_background_color h 255 255 255;
      Pl.set_foreground_color h 150 150 150; (* grid lines *)
      Pl.set_ylabel h "frequency of A allele";
      Pl.set_xlabel h "possible distributions";
      Pl.set_zlabel h "probability";
      Pl.set_altitude h altitude;
      Pl.set_azimuth h azimuth;
      Pl.mesh ~h xs ys zs;
      Pl.output h;
      Printf.printf "%s\n%!" filename
  in LL.iteri make_pdf (sub_lazy_list start_gen last_gen distlists)

  (*
      match altitude with Some a -> Pl.set_altitude h a | None -> ();
      match azimuth  with Some a -> Pl.set_azimuth  h a | None -> ();
  *)

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
