
(* NOT equation (3) in Skulj ?? *)
let flowers max_n min_alg max_alg =
  let flower (idxs, lower) = 
    let complement_upper = (L.assoc (list_complement max_n idxs) max_alg)
    in (idxs, (max lower (1. -. complement_upper)))
  in L.map flower min_alg

(* NOT equation (4) in Skulj ?? *)
let fuppers max_n min_alg max_alg =
  let fupper (idxs, upper) = 
    let complement_lower = (L.assoc (list_complement max_n idxs) min_alg)
    in (idxs, (min upper (1. -. complement_lower)))
  in L.map fupper max_alg


let rec subtract_list xs ys =
  match xs with 
  | [] -> []
  | x::more_xs -> match ys with
    | [] -> xs
    | y::more_ys -> 
      if x = y 
      then subtract_list more_xs more_ys
      else x::(subtract_list more_xs ys)

(* Note on next_intsets: An alternative would be to skip the append 
 * and just have each element contain the new additions to the power 
 * set.  Then to get a power set, you have to concat all sets up until 
 * and including that one.  This is more space efficient but means 
 * re-appending again if you want another powerset. *)

(* Note on algebra_probs: Could use Sets instead of lists as keys, but they don't 
 * display their contents by default, which makes playing around at the repl 
 * inconvenient. *)

(* Perhaps it's a bit fragile, but since I'm representing sets of
 * atoms as lists of integers in decreasing order, we can use that
 * to compute complements.  Since this is for exploratory
 * experimentation, I prefer this representation to using Batteries or Core
 * Sets, which don't display their contents by default.
*)


(** OLD DEFINITIONS OF pri_f_field_lowers/uppers  **)

(** OLD Calculate a lower, L value from a probability interval.
    Equation (3) p. 1317 in Skulj 2009 (funny capitalization evokes the paper) *)
let pri_F_field_Lower  omega_max  atom_mins  atom_maxs  subset_idxs = 
  let mins_sum = prob_sum atom_mins subset_idxs in
  let maxs_comp_sum = prob_sum atom_maxs (list_complement omega_max subset_idxs) in
  max mins_sum (1. -. maxs_comp_sum)

(** OLD Calculate an upper, U value from a probability interval.
    Equation (4) p. 1317 in Skulj 2009 *)
let pri_F_field_Upper  omega_max  atom_mins  atom_maxs  subset_idxs = 
  let maxs_sum = prob_sum atom_maxs subset_idxs in
  let mins_comp_sum = prob_sum atom_mins (list_complement omega_max subset_idxs) in
  min maxs_sum (1. -. mins_comp_sum)

(** OLD Calculate L values for all members of the algebra and return an
    (atoms, L-value) alist *)
let pri_F_field_Lowers omega_max atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets omega_max in
  let lower idx_set =
    (idx_set, pri_F_field_Lower omega_max atom_mins atom_maxs idx_set) in
  L.map lower algebra_idx_sets

(** OLD Calculate U values for all members of the algebra and return an
    (atoms, U-value) alist *)
let pri_F_field_Uppers omega_max atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets omega_max in
  let upper idx_set = 
    (idx_set, pri_F_field_Upper omega_max atom_mins atom_maxs idx_set) in
  L.map upper algebra_idx_sets

(** Given *two* algebra_probs lists, return a similar alist in which
    values are the maximum of the two corresponding probabilities. *)
let algebra_maxs = algebra_extrema max
