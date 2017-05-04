(** Tools for exploring imprecise Markov processes *)

(* RWO says you shouldn't do this, but Clojure shows it works well. *)
module M  = Owl.Mat
module L  = Batteries.List
module LL = Batteries.LazyList
module A  = Batteries.Array
(* module P = Pervasives *) (* to override defs from Owl.Mat *)

(* Pull out infixes I might want so I don't have to open Owl.Mat and shadow Pervasives: *)
let ($@) = Owl.Mat.($@)  (* = dot: matrix multiplication *)
let (+@) = Owl.Mat.(+@) (* matrix addition *)
let (-@) = Owl.Mat.(-@) (* matrix subtraction *)
let (/$) = Owl.Mat.(/$) (* divide matrix by float *)
let ($/) = Owl.Mat.($/) (* divide float by matrix *)
let ( *@ ) = Owl.Mat.( *@ )  (* = mul: elementwise multiplication *)

(*********** Ways to process matrices **********)

(** Apply f to all combinations of elements in xs and ys using fold_right.  
    Preserves order. *)
let cross_apply f xs ys =
  L.fold_right (fun x xacc -> L.fold_right
                   (fun y yacc -> (f x y)::yacc)
                   ys xacc)
    xs [];;

(** Given two lists of matrices, return a list containing the products of all 
    combinations of one matrix from the first list and the other from the second  list. *)
let mult_mats = cross_apply M.dot

(*********** Constructing imprecise probability vectors, algebras, and matrices **********)


(** Given a list of matrices return one with the same dimensions
    whose elements are the minima of corresponding locations. *)
let min_elts matlist = L.reduce M.min2 matlist

(** Given a list of matrices return one of the same dimensions
    whose elements are the maxima of corresponding locations. *)
let max_elts matlist = L.reduce M.max2 matlist

(** length of a row vector or number of columns of a matrix *)
let size m = snd (M.shape m)

(** Given a sequence of sequences representing the power set of non-negative
    integers up through n-1, returns a pair consisting of the sequence of 
    sequences for up through n-1 and one up through n.  (This kind of pair 
    is what LazyList.from_loop wants.)  Assumes that the first element of 
    the first element of the power set sequence pset passed in is n-1. *)
let next_intsets pset =
  let n = 1 + L.hd (L.hd pset) in  (* Get next integer; previous one must be first in the first element. *)
  let addl_sets = (L.map (fun xs -> n :: xs) pset) in
  (pset, addl_sets @ pset)

(** Return a lazy list of subsequent power sets of integers from 0 to n. 
    They can be retreived using e.g., to get the power set of integers
    up to 5: LazyList.at intsets 4 .  Note each power set is in the form
    of a regular list; only the top level list is lazy. *)
let make_intsets () = LL.from_loop [[0]; []] next_intsets
let algebra_sets = make_intsets ()

(** Given a vector of probabilities and a list of indexes representing
    atoms, returns sum of probabilities for the set containing those atoms. *)
let prob_sum probs atom_idxs =
  let add_another prob idx = prob +. M.get probs 0 idx   (* Owl.Mat.get rather than .{i,j} to get type right *)
  in L.fold_left add_another 0. atom_idxs

(** Given a probability vector, returns an alist of pairs for all
    elements in the algebra of sets based on the atoms represented by
    elements in the vector.  Each pair contains a list of indexes 
    representing atoms in a set followed by the probability of that set. *)
let algebra_probs probs = 
  let num_atoms = size probs in
  let idx_sets = LL.at algebra_sets (num_atoms - 1) in
  let idx_prob_entry idxs = (idxs, prob_sum probs idxs) in
  L.map idx_prob_entry idx_sets 

(** Given two algebra_probs alists, return a similar alist in which values are 
 * the minimum/maximum/etc (according to relat) of the two corresponding probs. *)
let algebra_extrema relat alg_probs1 alg_probs2 =
  let make_entry = 
    (fun set_prob1 set_prob2 ->
       let (atoms, prob1) = set_prob1 in
       let (_    , prob2) = set_prob2 in  (* first elts s/b same *)
       (atoms, relat prob1 prob2))
  in L.map2 make_entry alg_probs1 alg_probs2

(** Given two algebra_probs lists, return a similar list in which
    values are the minimum of the two corresponding probabilities. *)
let algebra_mins = algebra_extrema min

(** Given two algebra_probs lists, return a similar list in which
    values are the maximum of the two corresponding probabilities. *)
let algebra_maxs = algebra_extrema max

(** Given a list of algebra_probs, return an algebra_prob-like list
    with minima of all probs for each set of indexes. *)
let min_algebra_elts alg_probs_list =
  L.reduce (fun combo alg -> algebra_mins combo alg) alg_probs_list

(** Given a list of algebra_probs, return an algebra_prob-like list
    with maxima of all probs for each set of indexes. *)
let max_algebra_elts alg_probs_list =
  L.reduce (fun combo alg -> algebra_maxs combo alg) alg_probs_list

(** Given lists xs and ys that are both ordered in the same way (e.g. 
    monotonically ordered integers), return a list containing all
    elements of xs not in ys (with order preserved).  The elements
    in ys must be a (perhaps improper) subset of xs. *)
(* This version is by RichN at https://codereview.stackexchange.com/a/162407/61384
 * as a code review answer to my original version. *)
let rec subtract_list xs ys =
  match xs, ys with 
  | [], _ -> []
  | _, [] -> xs
  | x::xs', y::ys' when x = y -> subtract_list xs' ys'
  | x::xs', _ -> x::(subtract_list xs' ys)

(** Given a maximum element max_n and a subset of a domain of atoms 
    represented by list of integers in decreasing order, return a list
    of the integers, between 0 and max_n inclusive, that are not in
    the subset.  The returned list will also be in decreasing order.
    Note that max_n is one less than the size of the domain. *)
let list_complement max_n subset =
  let omega = L.range max_n `Downto 0 in
  subtract_list omega subset

(* let algebra_complements omega idx_sets = L.map (fun s -> subtract_list omega s) idx_sets *)


(** Calculate a lower, L value from a probability interval.
    Equation (3) p. 1317 in Skulj 2009 (funny capitalization evokes the paper) *)
let pri_F_field_Lower  max_n  atom_mins  atom_maxs  subset_idxs = 
  let mins_sum = prob_sum atom_mins subset_idxs in
  let maxs_comp_sum = prob_sum atom_maxs (list_complement max_n subset_idxs) in
  max mins_sum (1. -. maxs_comp_sum)

(** Calculate an upper, U value from a probability interval.
    Equation (4) p. 1317 in Skulj 2009 *)
let pri_F_field_Upper  max_n  atom_mins  atom_maxs  subset_idxs = 
  let maxs_sum = prob_sum atom_maxs subset_idxs in
  let mins_comp_sum = prob_sum atom_mins (list_complement max_n subset_idxs) in
  min maxs_sum (1. -. mins_comp_sum)

(** Calculate L values for all members of the algebra and return an
    (atoms, L-value) alist *)
let pri_F_field_Lowers max_n atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets max_n in
  let lower idx_set =
    (idx_set, pri_F_field_Lower max_n atom_mins atom_maxs idx_set) in
  L.map lower algebra_idx_sets

(** Calculate U values for all members of the algebra and return an
    (atoms, U-value) alist *)
let pri_F_field_Uppers max_n atom_mins atom_maxs =
  let algebra_idx_sets = LL.at algebra_sets max_n in
  let upper idx_set = 
    (idx_set, pri_F_field_Upper max_n atom_mins atom_maxs idx_set) in
  L.map upper algebra_idx_sets



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
let unif_stoch_vec len = normalize_vec (M.uniform 1 len)

(** Make a random uniform stochastic matrix of size dim x dim.
    It's each row that has a sum of 1, so this is to be the right 
    multiplicand with a row vector as a left multiplicand. *)
let unif_stoch_mat dim =
  let rows = A.init dim (fun _ -> unif_stoch_vec dim) in
  M.of_rows rows

(** Make a square matrix initializing it with function f, which is
    passed indexes, e.g: mat_from_fn 4 (fun _ -> Owl.Stats.Rnd.uniform ())  *)
let mat_from_fn dim f = M.map f (M.empty dim dim)


(*********** test data ***********)

let size = 4
let ps = L.init 10 (fun _ -> unif_stoch_vec size)
let mins = min_elts ps
let maxs = max_elts ps
let f_lowers = pri_F_field_Lowers (size - 1) mins maxs
let f_uppers = pri_F_field_Uppers (size - 1) mins maxs
let algs = L.map algebra_probs ps
let min_alg = min_algebra_elts algs
let max_alg = max_algebra_elts algs
