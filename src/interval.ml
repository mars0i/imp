
module M  = Owl.Mat
module I = Core.Interval;;

(* TODO is this right? *)
(** Return true iff pred is true for all corresponding elements of
    matrices m1 and m2. Short-circuits on the first false. *)
let forall2 pred m1 m2 =
  let rows, cols as dims = M.shape m1 in
  if dims <> M.shape m2 then raise (Failure "matrices have different shapes");
  let rec loop_cols i j =
    if j >= cols then true else
      pred (M.get m1 i j) (M.get m2 i j) && loop_cols i (j + 1)
  in let rec loop_rows i j =
    if i >= rows then true else
      loop_cols i j && loop_rows (i + 1) j
  in loop_rows 0 0

(** A compare function for matrices that returns zero of all elements of
    both matrices are equal, and if not returns -1 only if the matrices are 
    not equal and all elements of m1 are less than or equal to corresponding 
    elements of m2, and if not returns 1 only if all elements of m1 are 
    greater than or equal to corresponding elements of m2, and raises an
    exception of none of these relationships hold. (The normal compare
    function short-circuits on the first non-equal elements *)
let compare m1 m2 =
  if m1 = m2 then 0 else
    if forall2 (<=) m1 m2 then -1 else 
      if forall2 (>=) m1 m2 then 1 
      else raise (Failure "incomparable matrices")
