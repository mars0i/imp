
(** Calculate Wright-Fisher transition probabilities with selection.
    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

(* Example:
 * module LL = Batteries.LazyList;;
 * open Owl;;
 * let xs = Mat.sequential 1 1001;;
 * let s = make_init_state 1000 500;;
 * let t = make_tranmat 1000 (0.3, 0.3, 0.4);;
 * let h = Plot.create "yo.pdf" in Plot.scatter ~h xs (LL.at states 10); Plot.output h;;
 *)

module Mat = Owl.Mat
module Math = Owl.Maths (* note British->US translation *)
module LL = Batteries.LazyList

let ( *@ ) = Mat.( *@ )  (* = dot: matrix multiplication *)

(* TODO prob_ij inefficiently calls weight_i multiple times with same args,
 * i.e.  within each row. *)

(* TODO Also there's no need to calculate the combinations again for each
 * row in tranmat; they're the same in each row. *)

(* TODO Replace with Owl.Maths.combination_float when I get newer version. *)
(* Owl.Math.combination just wraps the following in a conversion to int.
 * This produces odd results for larger coefficients because of OCaml's unsafe,
 * limited-precision integers.  I need a float in the end anyway. *)
let combination_float n k = 
  Gsl.Sf.choose n k

let make_init_state allele_popsize num_alleles =
  let m = Mat.zeros 1 (allele_popsize + 1) in
  Mat.set m 0 num_alleles 1.0;
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

let next_state tranmat state = 
  (state, state *@ tranmat)

let make_states tranmat init_state =
  LL.from_loop init_state (next_state tranmat)

let length m = snd (M.shape m)

(* Make a series of n plot pdfs from states using basename. *)
let make_pdfs basename states n =
  let state_length = length (LL.at states 0) in
  let xs = Mat.sequential 1 state_length in (* vector of x-axis indices *)
  let f i state =
    let filename = basename ^ (string_of_int i) ^ ".pdf" in
    let h = Plot.create filename in 
      Plot.scatter ~h xs state; 
      Plot.output h
  in LL.iteri f (LL.take n states)
