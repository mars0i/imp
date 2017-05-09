(** Tools for exploring imprecise Markov processes *)

(* For names mentioned, see references section at end of file. *)

module M  = Owl.Mat
module A  = Batteries.Array
module L  = Batteries.List
module LL = Batteries.LazyList
module S  = Batteries.String


(*********** misc. convenience definitions **********)

(* Pull out infixes I might want so I don't have to open Owl.Mat and shadow Pervasives: *)
let ( *@ ) = M.( *@ )  (* = dot: matrix multiplication *)
let (/$) = M.(/$) (* divide matrix by float *)

(** length of a row vector or number of columns of a matrix *)
let size m = snd (M.shape m)


(*********** algebras of indexes representing atoms **********)

(** Create a power set of integers from a smaller power set.
  *
  * Given a sequence of sequences representing the power set of non-negative
  * integers up through n-1, returns a pair consisting of the sequence of 
  * sequences for up through n-1 and one up through n.  (This kind of pair 
  * is what LazyList.from_loop wants.)  Assumes that the first element of 
  * the first element of the power set sequence pset passed in is n-1. *)
let next_intsets pset =
  let n = 1 + L.hd (L.hd pset) in  (* Get next integer; previous one must be first in the first element. *)
  let addl_sets = (L.map (fun xs -> n :: xs) pset) in
  (pset, addl_sets @ pset)

(** Generate a list of subsequent integer power sets. 
 *
 *  Return a lazy list of subsequent power sets of integers from 0 to n. 
 *  They can be retreived using e.g., to get the power set of integers
 *  up to 5: LazyList.at intsets 4 .  Note each power set is in the form
 *  of a regular list; only the top level list is lazy. *)
let make_intsets () = LL.from_loop [[0]; []] next_intsets

(** A lazy list of integer power sets. *)
let algebra_sets = make_intsets ()

(** Set difference for ordered lists of integers.
  *
  * Given lists xs and ys that are both ordered in the same way (e.g. 
  * monotonically ordered integers), return a list containing all
  * elements of xs not in ys (with order preserved).  The elements
  * in ys must be a (perhaps improper) subset of xs.
  * This version by RichN at https://codereview.stackexchange.com/a/162407/61384 *)
let rec subtract_list xs ys =
  match xs, ys with 
  | [], _ -> []
  | _, [] -> xs
  | x::xs', y::ys' when x = y -> subtract_list xs' ys'
  | x::xs', _ -> x::(subtract_list xs' ys)

(** Set complement for ordered lists of integers.
  * i.e. return the atoms representing the negation of the original set.
  *
  * Given a maximum element omega_max and a subset of a domain of atoms 
  * represented by list of integers in decreasing order, return a list
  * of the integers, between 0 and omega_max inclusive, that are not in
  * the subset.  The returned list will also be in decreasing order.
  * Note that omega_max is one less than the size of the domain. *)
let list_complement omega_max subset =
  let omega = L.range omega_max `Downto 0 in
  subtract_list omega subset


(*********** imprecise probability vectors and matrices **********)

(** Given a list of matrices return one with the same dimensions
    whose elements are the minima of corresponding locations. *)
let min_elts matlist = L.reduce M.min2 matlist

(** Given a list of matrices return one of the same dimensions
    whose elements are the maxima of corresponding locations. *)
let max_elts matlist = L.reduce M.max2 matlist


(*********** probabilities over algebras **********)

(** Given a vector of atom probabilities and a list of indexes representing
    atoms, returns sum of probabilities for the set containing those atoms. *)
let prob_sum probs atom_idxs =
  let add_prob sum idx = sum +. M.get probs 0 idx  (* Owl.Mat.get rather than .{i,j} to get type right *)
  in L.fold_left add_prob 0. atom_idxs

(** Given a probability vector, returns an alist of pairs for all
    elements in the algebra of sets based on the atoms represented by
    elements in the vector.  Each pair contains a list of indexes 
    representing atoms in a set followed by the probability of that set. *)
let algebra_probs probs = 
  let i = (size probs) - 1 in
  let idx_sets = LL.at algebra_sets i in
  let make_entry event = (event, prob_sum probs event) in
  L.map make_entry idx_sets 


(*********** imprecise probabilities over algebras **********)

(** Given *two* algebra_probs alists, return a similar alist in which values are 
 * the minimum/maximum/etc (according to relat) of the two corresponding probs. *)
let algebra_extrema relat alg_probs1 alg_probs2:((int list * float) list) =  (* compiler needs a little help with types here *)
  let make_entry (event, prob1) (_, prob2) =
    (event, relat prob1 prob2) in (* first elts s/b same *)
  L.map2 make_entry alg_probs1 alg_probs2

(** Given a list of *multiple* algebra_probs, return an algebra_prob-like
    list with minima of all probs for each set of indexes. *)
let min_algebra_elts alg_probs_list =
  let min_combine combo alg = algebra_extrema min combo alg in
  L.reduce min_combine alg_probs_list

(** Given a list of *multiple* algebra_probs, return an algebra_prob-like
    list with maxima of all probs for each set of indexes. *)
let max_algebra_elts alg_probs_list =
  let max_combine combo alg = algebra_extrema max combo alg in
  L.reduce max_combine alg_probs_list

(** Return (1 - the sum of probs of extreme probs) for a set of atoms
    [for (3) and (4) in Skulj] *)
let invert_prob_sum omega_max atom_extrema subset_idxs = 
  1. -. prob_sum atom_extrema (list_complement omega_max subset_idxs)

(** Map prob_sum over each possible combination of atoms. *)
let simple_sums omega_max atom_extrema =
  L.map (prob_sum atom_extrema) (LL.at algebra_sets omega_max)

(** Map invert_prob_sum over each possible combination of atoms. *)
let inverted_sums omega_max atom_extrema =
  L.map (invert_prob_sum omega_max atom_extrema) (LL.at algebra_sets omega_max)

(** Calculate L values for all members of the algebra and return an
    (atoms, L-value) alist [(3) in Skulj] *)
let pri_f_field_lowers omega_max atom_mins atom_maxs =
  let mins = simple_sums omega_max atom_mins in
  let inverted_maxs = inverted_sums omega_max atom_maxs in
  let minmins = L.map2 max mins inverted_maxs in
  L.combine (LL.at algebra_sets omega_max) minmins

(** Calculate U values for all members of the algebra and return an
    (atoms, U-value) alist [(4) in Skulj] *)
let pri_f_field_uppers omega_max atom_mins atom_maxs =
  let maxs = simple_sums omega_max atom_maxs in
  let inverted_mins = inverted_sums omega_max atom_mins in
  let maxmaxs = L.map2 min maxs inverted_mins in
  L.combine (LL.at algebra_sets omega_max) maxmaxs

(** Given lists of atoms, value pairs from pri_f_field_lowers and
  * pri_f_field_uppers, return a list of pairs that combine the lower
  * and upper values into intervals represented as 2-element lists. *)
let pri_f_field_intervals lowers uppers =
  let add_elt (event, lower_prob) (_, upper_prob) elts =   (* events s/b same in lower and upper *)
    (event, [lower_prob; upper_prob]) :: elts
  in L.fold_right2 add_elt lowers uppers []


(*********** Ways to make matrices **********)

(** Given a matrix such as a row or column vector, return a version 
    normalized so that elements sum to 1.0. *)
let normalize_vec m =
  let tot = M.sum m in
  m /$ tot;;

(* The next definitions produce stochastic row vectors and stochastic
   matrices with rows that each sum to 1.  We multiply with vectors
   on the left and matrices on the right.  To get the other scheme,
   transpose these vectors and matrices. *)

(** Make a random uniform stochastic row vector of length len. *)
let unif_stoch_vec dim = normalize_vec (M.uniform 1 dim)

(** Make a random uniform stochastic matrix of size dim x dim.
    It's each row that has a sum of 1, so this is to be the right 
    multiplicand with a row vector as a left multiplicand. *)
let unif_stoch_mat dim =
  let rows = A.init dim (fun _ -> unif_stoch_vec dim) in
  M.of_rows rows

(** Make a square matrix initializing it with function f, which is
    passed indexes, e.g: mat_from_fn 4 (fun _ -> Owl.Stats.Rnd.uniform ())  *)
let mat_from_fn dim f = M.map f (M.empty dim dim)


(*********** Ways to process matrices **********)

(** Apply f to all combinations of elements in xs and ys using fold_right.  
    Preserves order. *)
let cross_apply f xs ys =
  L.fold_right (fun x xacc -> L.fold_right
                   (fun y yacc -> (f x y)::yacc)
                   ys xacc)
    xs [];;

(** Given two lists of matrices, return a list containing the products of all 
    combinations of one matrix from the first list and the other from the 
    second list. *)
let mult_mats = cross_apply M.dot


(*********** Strings for printing **********)

let string_of_t_list string_of_t l = 
  "[" ^ String.concat "; " (L.map string_of_t l) ^ "]"

let string_of_int_list = string_of_t_list string_of_int

let string_of_float_list = string_of_t_list string_of_float

(* Convert one algebra entry, a record containing a list of indexes
 * and a probibility, into a string. *)
let string_of_alg_prob alg_prob =
  let (k, v) = alg_prob in
  "(" ^ string_of_int_list k ^ ", " ^ string_of_float v ^ ")"

let string_of_alg_interval alg_interval =
  let (k, v) = alg_interval in
  "(" ^ string_of_int_list k ^  ", " ^ string_of_float_list v ^ ")"

(* Convert a list of indexes, probability entries into a string. *)
let string_of_alg_probs = string_of_t_list string_of_alg_prob

(* Convert a list of indexes, probability interval entries into a string. *)
let string_of_alg_intervals = string_of_t_list string_of_alg_interval


(*********** Convenience function for testing **********)

let generate_system omega_size num_dists =
  let omega_max = omega_size - 1 in
  (* make num_dists atomic dists *)
  let ps = L.init num_dists (fun _ -> unif_stoch_vec omega_size) in
  (* probabilities for algebras for each of the num_dists distributions *)
  let algs = L.map algebra_probs ps in (* alists mapping atom lists to probs *)
  (* min and max values of atomic probs across all num_dists distributions *)
  let mins = min_elts ps in
  let maxs = max_elts ps in
  (* min and max values of probs for each member of the algebra *)
  let min_alg = min_algebra_elts algs in
  let max_alg = max_algebra_elts algs in
  (* prob values for each member of the algebra computed using (3) in Skulj 
   * The first two are min'ed to produce the third. *)
  let f_mins = simple_sums omega_max mins in
  let f_inverted_maxs = inverted_sums omega_max maxs in
  let f_lowers = pri_f_field_lowers omega_max mins maxs in
  (* prob values for each member of the algebra computed using (4) in Skulj
   * The first two are max'ed to produce the third. *)
  let f_maxs = simple_sums omega_max maxs in
  let f_inverted_mins = inverted_sums omega_max mins in
  let f_uppers = pri_f_field_uppers omega_max mins maxs in
  (* interval enteries constructed from f_lowers, f_uppers *)
  let f_intervals = pri_f_field_intervals f_lowers f_uppers in
  (* return all of the above: *)
  ((ps, mins, maxs),
   (algs, min_alg, max_alg),
   (f_lowers, f_mins, f_inverted_maxs),
   (f_uppers, f_maxs, f_inverted_mins),
   f_intervals)




(*********** References **********)

(* Damjan Škulj, "Discrete time Markov chains with interval probabilities",
 * _International Journal of Approximate Reasoning_, Volume 50, Issue 8, 
 * 2009, Pages 1314-1329. *)
