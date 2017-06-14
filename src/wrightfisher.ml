
(** Calculate Wright-Fisher transition probabilities with selection.
    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

(* Example:
 * module LL = Batteries.LazyList;;
 * open Owl;;
 * let xs = Mat.sequential 1 1001;;
 * let s = make_init_dist 1000 500;;
 * let t = make_tranmat 1000 (0.3, 0.3, 0.4);;
 * let h = Plot.create "yo.pdf" in Plot.scatter ~h xs (LL.at dists 10); Plot.output h;;
 *)

module Mat = Owl.Mat
module Math = Owl.Maths (* note British->US translation *)
module Plot = Owl.Plot
module L = Batteries.List
module LL = Batteries.LazyList

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

let combination_float = Gsl.Sf.choose

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
  let comb = combination_float allele_popsize next_freq in
  comb  *.  wt ** j  *.  other_wt ** j'

(** prob_ij with an extra ignored argument; can be used mapi to
    initialize a matrix. *)
let prob_ijf fitns allele_popsize prev_freq next_freq _ =
  prob_ij fitns allele_popsize prev_freq next_freq

let make_tranmat allele_popsize fitns =
  let dim = allele_popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (prob_ijf fitns allele_popsize) m

let next_dist tranmat dist = 
  (dist, dist *@ tranmat)

(** Returns a LazyList of probability distributions that are the result of 
    a stationary Markov chain with initial distribution (usually with 1 at 
    one entry and zero elsewhere, representing that the population has that 
    initial frequency) and a transition matrix. The initial distribution is 
    not included in the list.*)
let make_dists tranmat init_dist =
  LL.drop 1 (LL.from_loop init_dist (next_dist tranmat))

(** take n elements from a LazyList and convert the result to a List. *)
let take_to_list n ll = 
  LL.to_list (LL.take n ll)

let length m = snd (M.shape m)

(* Make a series of n plot pdfs from dists using basename. *)
let make_pdfs basename dists n =
  let dist_length = length (LL.at dists 0) in
  let xs = Mat.sequential 1 dist_length in (* vector of x-axis indices *)
  let make_pdf i dist =
    let filename = basename ^ (Printf.sprintf "%03d" i) ^ ".pdf" in
    let h = Plot.create filename in 
    Plot.set_yrange h 0.0 0.25;
    Plot.scatter ~h xs dist; 
    Plot.output h
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
  (dists, L.concat (L.map (fun dist -> L.map (M.dot dist) tranmats)
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
let make_dist_lists tranmats init_dists =
  LL.from_loop init_dists (next_dists tranmats)


let sort_dists dists = L.sort Utils.difference_compare dists


(*
 
Note this is a small pop below, with big effects from initial frequencies.
Better to use at least N=100, or more depending on intensity of selection.

let t0 = make_tranmat 50 {w11=0.7; w12=0.8; w22=1.0};;
let t1 = make_tranmat 50 {w11=0.7; w12=0.3; w22=1.0};;

let s = make_init_dist 50 25;;
let t0 = make_tranmat 50 {w11=1.0; w12=0.8; w22=0.7};;
let t1 = make_tranmat 50 {w11=0.7; w12=0.8; w22=1.0};;
let distlists = make_dist_lists [t0; t1] [s];;
let xs, ys, zs = make_coords (LL.at distlists 2) in let h = Plot.create "yo.pdf" in Plot.mesh ~h xs ys zs; Plot.output h;;

i.e.:
let xs, ys, zs = make_coords (LL.at distlists 2) in
let h = Plot.create "yo.pdf" in
   Plot.mesh ~h xs ys zs;
  Plot.output h

or:
let n = 3 in let xs, ys, zs = make_coords (LL.at distlists2 n) in let h = Plot.create "yo.pdf" in Plot.mesh ~h xs ys zs; Plot.output h;;
i.e.:
let gen = 3 in
let xs, ys, zs = make_coords (LL.at distlists2 gen) in
let h = Plot.create "yo.pdf" in
  Plot.mesh ~h xs ys zs;
  Plot.output h;;

*)
  
