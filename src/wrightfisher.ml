
(** Calculate Wright-Fisher transition probabilities with selection.
    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

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

(** 1.59 *)
let weight_i allele_popsize fitnesses freq =
  let i, i' = float freq, float (allele_popsize - freq) in
  let w11, w12, w22 = fitnesses in
  let a_hom, het, b_hom = w11 *. i *. i,
                          w12 *. i *. i',
                          w22 *. i' *. i' in
  (a_hom +. het) /. (a_hom +. 2. *. het +. b_hom)

(** Wright-Fisher transition probability from frequency = i to frequency = j *)
let prob_ij allele_popsize fitnesses prev_freq next_freq =
  let wt = weight_i allele_popsize fitnesses prev_freq in
  let other_wt = 1. -. wt in
  let j = float next_freq in
  let j' = float (allele_popsize - next_freq) in
  let comb = combination_float allele_popsize next_freq in
  comb  *.  wt ** j  *.  other_wt ** j'

  (** prob_ij with an extra ignored argument; can be used mapi to
   * initialize a matrix. *)
let prob_ijf allele_popsize fitnesses prev_freq next_freq _ =
  prob_ij allele_popsize fitnesses prev_freq next_freq

let make_tranmat allele_popsize fitnesses =
  let dim = allele_popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (prob_ijf allele_popsize fitnesses) m

let next_state tranmat state = 
  (state, state *@ tranmat)

let make_states tranmat init_state =
  LL.from_loop init_state (next_state tranmat)
