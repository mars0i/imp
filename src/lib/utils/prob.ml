(** Tools for exploring imprecise Markov processes *)

(* For names mentioned, see references section at end of file. *)

module M  = Owl.Mat
module A  = Batteries.Array
module L  = Batteries.List


(*********** misc. convenience definitions **********)

(* Pull out infixes I might want so I don't have to open Owl.Mat and shadow Pervasives: *)
let ( *@ ) = M.( *@ )  (* = dot: matrix multiplication *)
let (/$) = M.(/$) (* divide matrix by float *)



(*********** imprecise probability vectors and matrices **********)

(** Given a list of matrices return one with the same dimensions
    whose elements are the minima of corresponding locations. *)
let min_elts matlist = L.reduce M.min2 matlist

(** Given a list of matrices return one of the same dimensions
    whose elements are the maxima of corresponding locations. *)
let max_elts matlist = L.reduce M.max2 matlist



(*********** imprecise probabilities over algebras **********)

(** Given *two* algebra_probs alists, return a similar alist in which values are 
   the minimum/maximum/etc (according to relat) of the two corresponding probs. *)
let algebra_extrema relat alg_probs1 alg_probs2:((int list * float) list) =  (* compiler needs a little help with types here *)
  let make_entry (event, prob1) (_, prob2) =
    (event, relat prob1 prob2) in (* first elts s/b same *)
  L.map2 make_entry alg_probs1 alg_probs2

(** Given a list of *multiple* algebra_probs, return an algebra_prob-like
    list with minima of all probs for each set of indexes. *)
let min_algebra_probs alg_probs_list =
  let min_combine combo alg = algebra_extrema min combo alg in
  L.reduce min_combine alg_probs_list

(** Given a list of *multiple* algebra_probs, return an algebra_prob-like
    list with maxima of all probs for each set of indexes. *)
let max_algebra_probs alg_probs_list =
  let max_combine combo alg = algebra_extrema max combo alg in
  L.reduce max_combine alg_probs_list

(** Given lists of atoms, value pairs from pri_f_field_lowers and
    pri_f_field_uppers, return a list of pairs that combine the lower
    and upper values into intervals represented as 2-element lists. *)
let pri_f_field_intervals lowers uppers =
  let add_elt (event, lower_prob) (_, upper_prob) elts =   (* events s/b same in lower and upper *)
    (event, (lower_prob, upper_prob)) :: elts
  in L.fold_right2 add_elt lowers uppers []
