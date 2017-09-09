(** This file contains versions of functions before they were
    optimized further.  The intent of these functions is clearer than
    the optimized versions; studying these functions can make it
    easier to understand the optimized versions. *)

(******************************************************)
(** Alternate versions of hilo_mult copied from git commit 4940ea6:
    Given [recombine_lo] or [recombine_hi], the original tight interval bounds
    P and Q, and either the previous tight component lo or hi bound (as
    appropriate), return the next lo or hi tight component bound.
    NOTE:
      If recomb is recombine_lo, the arguments should be P, Q, and the 
      previous lo matrix.  
      If recomb is recombine_hi, the arguments should be (notice!) Q, P,
      and the previous hi matrix. *)

(** This version is close to the earliest full working version. 
 *  The mapping to Hartfiel's algorithm descriptions is clearest. *)
let hilo_mult recomb p_mat q_mat prev_bound_mat = 
  let (m, n) = M.shape p_mat in
  let prev_cols = M.to_cols prev_bound_mat in
  let p_rows = M.to_rows p_mat in
  let q_rows = M.to_rows q_mat in
  let new_bound_mat = M.empty m n in
  for j = 0 to n - 1 do    (* j indexes columns in L and in p and q, as on pp. 50f *)
    for i = 0 to m - 1 do  (* and i indexes rows in p, q *)
      let prev_col = prev_cols.(j) in 
      let bar_row = recomb p_rows.(i) q_rows.(i) prev_col in
      M.(set new_bound_mat i j (get (bar_row *@ prev_col) 0 0)) (* result of multiplication is 1x1 *)
    done
  done;
  new_bound_mat

(** This version is closest to the Parmap-basedd def I used in the end. *)
let hilo_mult  recomb p_mat q_mat prev_bound_mat = 
  let (m, n) = M.shape p_mat in
  let len = m * n in
  let bounds_array = A.init len (calc_bound_val_from_flat_idx recomb p_mat q_mat prev_bound_mat m)
  in
  M.of_array bounds_array m n

(******************************************************)
(** Alternate version of recombine copied from git commit 33c1218 *)

(** Return the sum of all values in matrix except the one at i j. *)
let sum_except mat i j = 
  (M.sum mat) -. (M.get mat i j)

(** Given a relation (>=), a column l vec and two tight row vecs p and q s.t. 
    p<=q, return a stochastic row vec ("p bar") with high values from q where l
    is low and low values from p where l is high.  Or pass (<=), l, and tight
    row vecs s.t. p >= q to return a stoch row vec ("q bar") with low values 
    from p where l is low.  Note that the latter swaps the normal meanings of 
    p and q in Hartfiel, i.e. here the arguments should be (<=), l, q, p
    according to the normal senses of p and q. *)
let recombine_old relation p q lh =
  let pbar = M.clone p in
  let rec find_crossover idxs =
    match idxs with
    | i::idxs' -> 
        let qi = M.get q 0 i in
        let sum_rest = sum_except pbar 0 i in (* pbar begins <= 1 if p<=q, or >= 1 if p, q swapped *) (* TODO Can this be made more efficient?? *)
        if relation (qi +. sum_rest) 1.
        then M.set pbar 0 i (1. -. sum_rest) (* return--last iter put it over/under *)
        else (M.set pbar 0 i qi;             (* still <= 1, or >=1; try next one *)
          find_crossover idxs') 
    | [] -> raise (Failure "bad vectors") (* this should never happen *)
  in 
  find_crossover (idx_sort lh);
  pbar
