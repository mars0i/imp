(** pops.ml *)

module L = Batteries.List
module M = Owl.Mat

(** Returns a row vector representing the state of a population with two 
    competing types at t0.  Each element represents a possible absolute
    frequency for the focal type A.  For example, if 1.0 is the third element
    and all other elements are zero (which they should be), that means that
    there are 2 A's in a population of size n.  The result is a row vector 
    of length pop_size *)
let init_freq_dist pop_size freq =
  let v = M.zeros 1 (pop_size + 1) in
  M.set v 0 freq 1.0;
  v


let make_tran_mat step_probs =
  let dim = 1 + L.length step_probs in
  let pop_size = dim - 1 in
  let m = M.zeros dim dim in
  (* first and last rows: *)
  M.set m 0 0 1.0;
  M.set m pop_size pop_size 1.0;
  (* middle rows: *)
  for row = 1 to pop_size - 1 do
    (* step_probs has the length of a row.
       we can to shift to correct place, multiply by freqs,
       and sum whatevers extra at the end of the line. *)
    
  done

