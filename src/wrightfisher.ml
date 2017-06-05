
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

(* TODO seems to work for small pops but not large.  Is there a problem
 * with Math.combination?? 
 *
 * Yeah, maybe ...
 *
 * # combination 100 25;;
 * - : int = 0
 * # (factorial 100) /. ( (factorial 25)  *. (factorial 75));;
 * - : float = 2.42519269720337125e+23
 *
 * This comes from Gnu Scientific Lib btw.
 *
 * But there maybe OCaml-specific problems:
 * # combination 100 84;;
 * - : int = 1345860629046814720
 * # combination 100 83;;
 * - : int = -2573237163917575168
 *
 * It's wrapping around!
 * Maybe I need Clojure for this ....
 *
 * Is there a different algorithm I can use?  Even if you do it more
 * efficiently than with factorials The top divided by one on the bottom
 * Is just a shorter multiplication, there might still be problems
 *
 * Perhaps the problem is partly in the conversion to int from the
 * GSL's floats:
 *
 * # Gsl.Sf.choose 100 25;;
 * - : float = 2.42519269720337091e+23
 * # int_of_float (Gsl.Sf.choose 100 25);;
 * - : int = 0
 * The latter is the def of comination in Owl.
 *
 * *)

(* Owl.Math.combination just wraps the following in a conversion to int.
 * This doesn't work for larger coefficients, which messes up the tran probs.
 * I need a float in the end anyway, so use this instead. *)
let float_comb n k = 
  Gsl.Sf.choose n k

(* TODO inefficiently calls weight_i multiple times with same args. *)
(** Wright-Fisher transition probability from frequency = i to frequency = j *)
let prob_ij popsize fitnesses prev_freq next_freq =
  let alleles = 2 * popsize in
  let wt = weight_i popsize fitnesses prev_freq in
  let other_wt = 1. -. wt in
  let j = float next_freq in
  let j' = float (alleles - next_freq) in
  let comb = float_comb alleles next_freq in
  comb  *.  wt ** j  *.  other_wt ** j'

(*
(** prob_ij with an extra argument that's ignored.  For use with mapi. *)
let prob_ijf popsize fitnesses prev_freq next_freq ignored_float =
  prob_ij popsize fitnesses prev_freq next_freq
*)

let tranmat popsize fitnesses =
  let dim = 2 * popsize + 1 in
  let m = Mat.empty dim dim  in
  Mat.mapi (fun row col _ -> prob_ij popsize fitnesses row col) m
