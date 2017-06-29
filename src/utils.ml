module M  = Owl.Mat
module L = Batteries.List

(** Returns a memoizing version of function f of one argument.
    By Andrej Bauer: https://stackoverflow.com/a/14503530/1455243
    Caveats: 
    Only works for one arg.
    Not designed for recursive functions.  See above URL if you want that.
    Uses an association list, so if you have a lot of different results,
    could be inefficient.  Consider replacing with a Map or HashTable
    in that case. *)
let memo f =
  let m = ref [] in
    fun x -> try L.assoc x !m with 
              Not_found -> let y = f x in
              m := (x, y) :: !m ;
              y

(** Return true iff pred is true for all corresponding elements of
    matrices m1 and m2. Short-circuits on the first false. *)
let forall2 pred m1 m2 =
  let rows, cols as dims = M.shape m1 in
  if dims <> M.shape m2 then failwith "matrices have different shapes"
  ;
  let rec loop_cols i j =
    if j >= cols then true else
      pred (M.get m1 i j) (M.get m2 i j) && loop_cols i (j + 1)
  in let rec loop_rows i j =
    if i >= rows then true else
      loop_cols i j && loop_rows (i + 1) j
  in loop_rows 0 0

(** Fold f over matrices m1 and m2 starting with initial value init: 
    Folds f through all corresponding pairs of elements of matrices m1 
    and m2 by repeatedly applying f acc element_from_m1 element_from_m2,
    where acc is the result of previous applications.  init is the
    initial value for acc. *)
let fold2 f init m1 m2 =
  let rows, cols as dims = M.shape m1 in
  if dims <> M.shape m2 then failwith "matrices have different shapes"
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

(** Fold f over matrices m1 and m2 starting with initial value init, 
    short-circuiting if stop_val is encountered:
    Folds f through all corresponding pairs of elements of matrices m1 
    and m2 by repeatedly applying f acc element_from_m1 element_from_m2,
    where acc is the result of previous applications.  init is the
    initial value for acc.  If f ever returns stop_val, it will be
    returned immediately. *)
let short_circuit_fold2 stop_val f init m1 m2 =
  let rows, cols as dims = M.shape m1 in
  if dims <> M.shape m2 then failwith "matrices have different shapes"
  ;
  let last_col = cols - 1 in
  let apply_f acc i j = 
    f acc (M.get m1 i j) (M.get m2 i j)
  in
  let rec loop acc i j =
    if acc = stop_val then stop_val
    else if i < rows
    then loop (apply_f acc i j) (i + 1) j
    else if j < last_col         (* don't start on next col if at final col *)
         then loop acc 0 (j + 1) (* start over on next col *)
         else acc
  in loop init 0 0

(* NOTE The compare function below is complicated by the fact that it returns 0 for 
 * equivalent matrices.  However, if it's only used interval-creation,
 * a pair of equal matrices create an Empty interval, at least in
 * Jane Street-style Interval modules.  So you might as well return 1
 * for those.  -1 is the only return value that matters. *)


(** A compare function for matrices that returns zero if all elements of
    both matrices are equal, and if not returns -1 only if all elements 
    of m1 are less than or equal to corresponding elements of m2; otherwise
    returns 1, indicating that at least one element in m1 is greater than 
    the corresponding element in m2. *)
let biased_compare m1 m2 =
  let f acc e1 e2 =
    if e1 > e2 then 1     (* don't need to check for acc = 1 since fold exits first *)
    else match acc with   (* at this point we know that e1 <= e2 *)
         | -1 -> -1       (* all previous pairs were <= *)
         |  0 -> if e1 = e2 then 0 else -1
         |  _ -> failwith "bug: acc is not -1, 0, or 1" (* avoid match warning *)
  in short_circuit_fold2 1 f 0 m1 m2


(** A compare function for matrices that returns zero if all elements of
    both matrices are equal, -1 if each element of the first is less
    than or equal to the corresponding element of the second, or 1 if
    each element of the first is greater than or equal to each element
    of the second.  Raises an exception otherwise.  Note that the size
    of the differences between the values has no effect. *)
let standard_compare m1 m2 =
  let f acc e1 e2 =
    match acc with
         |  0 -> if e1 < e2 then -1 else
                 if e1 > e2 then 1 else 0
         | -1 -> if e1 <= e2 then -1 else failwith "incomparable"
         |  1 -> if e1 >= e2 then  1 else failwith "incomparable"
         |  _ -> failwith "bug: acc is not -1, 0, or 1" (* avoid match warning *)
  in fold2 f 0 m1 m2

(* Owl.Mat.signum, which converts each element into its sign, i.e.
 * -1, 0, or 1, might also be useful below. *)

(** Subtract mat2 from mat1 and sum the result.  A sort of poor person's
    non-normalized integral of the difference between the matrices
    (which might be vectors). *)
let sumdist mat1 mat2 = M.(sum (mat1 - mat2))

(** A compare function for matrices that determines whether the summed
    differences between corresponding matrix elements is negative, zero, 
    or positive.  This differs from standard_compare, which ignores sizes of
    differences.  With this compare function, by contrast, a large
    difference on one value can override smaller differences on other
    values. *)
let difference_compare m1 m2 =
  let dif = sumdist m1 m2 in
  if dif > 0.0 then 1
  else if dif < 0.0 then -1
  else 0

(** Given a matrix return a narrower matrix in which each every_nth element
    in each row is present.  The intervening elements are ignored.  *)
let subsample_in_rows every old_mat =
  if every <= 1 then old_mat
  else let (height, width) = M.shape old_mat in
       let new_width = width / every in (* WHAT ABOUT if doesn't divide evenly? *)
       let new_mat = M.empty height new_width in
       for i = 0 to height - 1 do
         for j = 0 to new_width - 1 do
           new_mat.{i,j} <- old_mat.{i, j*every}
         done
       done;
       new_mat

let insert_after n new_elt l = 
  if n = -1 then new_elt::l else
  L.fold_righti 
    (fun i elt acc -> if n = i then elt::new_elt::acc else elt::acc)
    l []

let insert_before n new_elt l = 
  if n = L.length l then l @ [new_elt]
  else L.fold_righti 
    (fun i elt acc -> if n = i then new_elt::elt::acc
                      else elt::acc)
    l []
