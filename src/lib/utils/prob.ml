(** Tools for exploring imprecise Markov processes *)

(* For names mentioned, see references section at end of file. *)

module M  = Owl.Mat
module A  = Batteries.Array
module L  = Batteries.List


(*********** misc. convenience definitions **********)

(* Pull out infixes I might want so I don't have to open Owl.Mat and shadow Pervasives: *)
let ( *@ ) = M.( *@ )  (* = dot: matrix multiplication *)
let (/$) = M.(/$) (* divide matrix by float *)


(*********** algebras of indexes representing atoms **********)

(** Set difference for ordered lists of integers.
   
    Given lists xs and ys that are both ordered in the same way (e.g. 
    monotonically ordered integers), return a list containing all
    elements of xs not in ys (with order preserved).  The elements
    in ys must be a (perhaps improper) subset of xs.
    This version by RichN at https://codereview.stackexchange.com/a/162407/61384 *)
let rec subtract_list xs ys =
  match xs, ys with 
  | [], _ -> []
  | _, [] -> xs
  | x::xs', y::ys' when x = y -> subtract_list xs' ys'
  | x::xs', _ -> x::(subtract_list xs' ys)

(** Set complement for ordered lists of integers.
    i.e. return the atoms representing the negation of the original set.
   
    Given a maximum element omega_max and a subset of a domain of atoms 
    represented by list of integers in decreasing order, return a list
    of the integers, between 0 and omega_max inclusive, that are not in
    the subset.  The returned list will also be in decreasing order.
    Note that omega_max is one less than the size of the domain. *)
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

(** Return (1 - the sum of probs of extreme probs) for a set of atoms
    See (3) and (4) in Skulj. *)
let invert_prob_sum omega_max atom_extrema subset_idxs = 
  1. -. prob_sum atom_extrema (list_complement omega_max subset_idxs)

(** Given lists of atoms, value pairs from pri_f_field_lowers and
    pri_f_field_uppers, return a list of pairs that combine the lower
    and upper values into intervals represented as 2-element lists. *)
let pri_f_field_intervals lowers uppers =
  let add_elt (event, lower_prob) (_, upper_prob) elts =   (* events s/b same in lower and upper *)
    (event, (lower_prob, upper_prob)) :: elts
  in L.fold_right2 add_elt lowers uppers []


(*********** Ways to create matrices **********)

(** Given a matrix such as a row or column vector, return a version 
    normalized so that elements sum to 1.0. *)
let normalize_vec m =
  let tot = M.sum' m in
  M.div_scalar m tot;;

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

