
(** Calculate Wright-Fisher transition probabilities with selection.
    Based on Ewens _Mathematical Population Genetics I, 2nd ed, 
    equations 1.58, 1.59, and 1.25, though similar formulas can be 
    found in many places. *)

(* TODO create transition matrix from pij calls. *)

module type WFParams = sig
  val popsize : int
  val fitnesses : float * float * float
end

module WF (Params : WFParams) = struct
  let popsize = Params.popsize
  let fitnesses = Params.fitnesses

  let weight_i popsize fitnesses freq =
    let i, i' = float freq, float ((2 * popsize) - freq) in
    let w11, w12, w22 = fitnesses in
    let a_hom, het, b_hom = w11 *. i *. i,
                            w12 *. i *. i',
                            w22 *. i' *. i' in
    (a_hom +. het) /. (a_hom +. 2. *. het +. b_hom)

  (** Wright-Fisher transition probability from frequency = i to frequency = j *)
  let pij popsize fitnesses prev_freq next_freq =
    let alleles = 2 * popsize in
    let wt = weight_i popsize fitnesses prev_freq in
    let other_wt = 1. -. wt in
    let j = float next_freq in
    let j' = float (alleles - next_freq) in
    let comb = float (Owl.Maths.combination alleles next_freq) in
    comb  *.  wt ** j  *.  other_wt ** j'
end
