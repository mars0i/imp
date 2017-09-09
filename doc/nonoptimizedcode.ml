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
