
(** Calculate Wright-Fisher transition probabilities with selection.
    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

module Mat = Owl.Mat
module Math = Owl.Maths (* note British->US translation *)

(** 1.59 *)
let weight_i popsize fitnesses freq =
  let i, i' = float freq, float ((2 * popsize) - freq) in
  let w11, w12, w22 = fitnesses in
  let a_hom, het, b_hom = w11 *. i *. i,
                          w12 *. i *. i',
                          w22 *. i' *. i' in
  (a_hom +. het) /. (a_hom +. 2. *. het +. b_hom)

(* Owl.Math.combination just wraps the following in a conversion to int.
 * This produces odd results for larger coefficients because of OCaml's unsafe,
 * limited-precision integers.  I need a float in the end anyway. *)
let float_comb n k = 
  Gsl.Sf.choose n k

(* TODO prob_ij inefficiently calls weight_i multiple times with same args,
 * i.e.  within each row. *)

(* TODO Also there's no need to calculate the combinations again for each
 * row in tranmat; they're the same in each row. *)

(** Wright-Fisher transition probability from frequency = i to frequency = j *)
let prob_ij popsize fitnesses prev_freq next_freq =
  let alleles = 2 * popsize in
  let wt = weight_i popsize fitnesses prev_freq in
  let other_wt = 1. -. wt in
  let j = float next_freq in
  let j' = float (alleles - next_freq) in
  let comb = float_comb alleles next_freq in
  comb  *.  wt ** j  *.  other_wt ** j'

  (** prob_ij with an extra ignored argument; can be used mapi to
   * initialize a matrix. *)
let prob_ijf popsize fitnesses prev_freq next_freq _ =
  prob_ij popsize fitnesses prev_freq next_freq

let tranmat popsize fitnesses =
  let dim = 2 * popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (prob_ijf popsize fitnesses) m
