
module M  = Owl.Mat
module I = Core.Interval;;

(* Something like this is what one wants, I think:
module Mat_interval = Interval.Make(struct 
  type t = Owl.Mat.mat
  let compare = Comparemat.compare
end);;
*)

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
(* TODO This is inefficient.  It potentially loops through all pairs
 * of elements three times.  
 * A new strategy:
 * Look at this element.  If it's =, then keep looking for anything.
 * If this element is not zero, then if it's <, all subsequent pairs
 * must be < or =; else raise an exception.
 * If this element is >, all subsequent elements must be > or =; else
 * raise an exception.
 * If we are still in the = state when we get done, then return 0.
 * If we are in the < or > state, then return -1 or 1, respectively.
 * So maybe have two functions or branchs:
 * The equality function/branch, which is open-minded, but that will
 * choose the other function/branch if it encounters a non-equal pair.
 * To this other fn/branch is passed a relation, <= or >=, which is
 * uses repeatedly. *)

(* TODO? An alternative would be to return 1 rather than raising an
 * exception when the matrices are incomparable.  This is an odd choice,
 * conceptually, but since when 1, i.e. when m1 > m2, we have no interval--
 * i.e. the interval creation function will return None, we could just
 * return None in the incomparable case, too, by returning 1.  i.e. in
 * general 1 essentially means that you'll get a non-interval, so why not
 * just piggy-back on this behavior? *)
(* i.e. how about this function (which is still inefficient, though less so)? *)
let compare' m1 m2 =
  if m1 = m2 then 0 else
    if forall2 (<=) m1 m2 then -1 else 1
(* or even: *)
let compare'' m1 m2 =
  if forall2 (<=) m1 m2 then -1 else 1
(* That misrepresents when m1 = m1.  Is that a problem? *)
