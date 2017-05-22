module M  = Owl.Mat

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

(* Fold f over matrices m1 and m2 starting with initial value init. *)
let fold2 f init m1 m2 =
  let rows, cols as dims = M.shape m1 in
  if dims <> M.shape m2 then raise (Failure "matrices have different shapes")
  ;
  let last_col = cols - 1 in
  let apply_f acc i j = 
    f acc (M.get m1 i j) (M.get m2 i j)
  in
  let rec loop acc i j =
    if i < rows
    then loop (apply_f acc i j) (i + 1) j
    else if j < last_col         (* don't start on next col if at final col *)
         then loop acc 0 (j + 1) (* start over on next col *)
         else acc
  in loop init 0 0