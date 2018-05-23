(** General-purpose utility functions for matrices, sequences, etc. *)

module M  = Owl.Mat
module A = Batteries.Array
module L = Batteries.List
module F = Core.Float

let always_true _ = true

let (%) f g = (fun x -> f (g x))

(********************************************)
(** Timing functions *)

(** Measure execution time of function of one argument.  Note that if used 
    with partial application on a function that expects multiple arguments,
    you might just get the time needed to return the first function, which
    would probably not be what you wanted. *)
let time0 f () =
    let cpu_time, wall_time = Sys.time(), Unix.gettimeofday() in
    let result = f in
    Printf.printf "cpu: %fs, wall: %fs\n%!" (Sys.time() -. cpu_time) (Unix.gettimeofday() -. wall_time);
    result

let time1 f x =
    let cpu_time, wall_time = Sys.time(), Unix.gettimeofday() in
    let result = f x in
    Printf.printf "cpu: %fs, wall: %fs\n%!" (Sys.time() -. cpu_time) (Unix.gettimeofday() -. wall_time);
    result

let time2 f x y =
    let cpu_time, wall_time = Sys.time(), Unix.gettimeofday() in
    let result = f x y in
    Printf.printf "cpu: %fs, wall: %fs\n%!" (Sys.time() -. cpu_time) (Unix.gettimeofday() -. wall_time);
    result

let time3 f x y z =
    let cpu_time, wall_time = Sys.time(), Unix.gettimeofday() in
    let result = f x y z in
    Printf.printf "cpu: %fs, wall: %fs\n%!" (Sys.time() -. cpu_time) (Unix.gettimeofday() -. wall_time);
    result

(********************************************)
(** Memo functions *)

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

(********************************************)
(** Sequence manipulators, type converters, etc. *)
let list_range ?(step=1) start stop = 
  if start = stop then [start]
  else let adjust_by, ineq =
    if stop > start then (+), (>) else (-), (<)
  in let rec aux curr stop' =
     if ineq curr stop' then []
     else curr::(aux (adjust_by curr step) stop')
  in aux start stop

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

(********************************************)
(** Matrix comparison functions *)

(* NOTE The compare function below is complicated by the fact that it returns
 * 0 for equivalent matrices.  However, if it's only used for 
 * interval-creation, a pair of equal matrices create an Empty interval, at 
 * least in Jane Street-style Interval modules.  So you might as well return 1
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

let make_compare diff m1 m2 = 
  let dif = diff m1 m2 in
  if dif > 0.0 then 1
  else if dif < 0.0 then -1
  else 0

(** Subtract mat2 from mat1 and sum the result.  A sort of poor person's
    non-normalized integral of the difference between the matrices
    (which might be vectors). *)
let sumdiff mat1 mat2 = M.(sum' (mat1 - mat2))

(** A compare function for matrices that determines whether the summed
    differences between corresponding matrix elements is negative, zero, 
    or positive.  This differs from standard_compare, which ignores sizes of
    differences.  With this compare function, by contrast, a large
    difference on one value can override smaller differences on other
    values. *)
let difference_compare = make_compare sumdiff

let absdiff mat1 mat2 = M.(sum' (abs (mat1 - mat2)))
let absdiff_compare = make_compare absdiff

(** L2 distance between mat1 and mat2: Subtract corresponding elements,
    square the results, sum those, and take the square root of the sum. *)
let l2diff mat1 mat2 = M.(l2norm' (mat1 - mat2))

(** A compare function for matrices based on the L2 distance. *)
let l2_compare = make_compare l2diff


(********************************************)
(** Matrix manipulation functions *)

(** Given a matrix, return a narrower matrix in which each [every] nth element
    in each row is present.  The intervening elements are ignored.  *)
let subsample_in_rows skip old_mat =
  if skip <= 1 then old_mat
  else let (height, width) = M.shape old_mat in
       let new_width = width / skip in (* WHAT ABOUT if doesn't divide evenly? *)
       let new_mat = M.empty height new_width in
       for i = 0 to height - 1 do
         for j = 0 to new_width - 1 do
           M.set new_mat i j (M.get old_mat i (j*skip))
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

let mapmap f outer = L.map (fun inner -> L.map f inner) outer

(********************************************)
(** Ways to process multiple matrices *)

(** Apply f to all combinations of elements, i.e. to the Cartesian product of
    all elements in xs and ys using fold_right.  Preserves order. *)
let cross_apply_lists f xs ys =
  L.fold_right (fun x xacc -> L.fold_right
                   (fun y yacc -> (f x y)::yacc)
                   ys xacc)
    xs [];;

(** Given two lists of matrices, return a list containing the products of all 
    combinations of one matrix from the first list and the other from the 
    second list. *)
let mult_mat_lists = cross_apply_lists M.dot

(** Apply f to all combinations of elements, i.e. to the Cartesian product of
    all elements in xs and ys using fold_right.  Preserves order. *)
let cross_apply_arrays f xs ys =
  A.fold_right (fun x xacc -> A.fold_right
                   (fun y yacc -> (f x y)::yacc)
                   ys xacc)
    xs [];;

(** Given two lists of matrices, return a list containing the products of all 
    combinations of one matrix from the first list and the other from the 
    second list. *)
let mult_mat_arrays = cross_apply_arrays M.dot


(********************************************)
(** Functions to generate matrices from lists *)

(** Converts a 1xN matrix, i.e. row vector, into a list. *)
let vec_to_list = Batteries.(A.to_list % M.to_array)

let list_to_vec l =
  let a = A.of_list l in
  M.of_array a 1 (A.length a)

(** Converts an MxN matrix into a list of M lists of length N. *)
let mat_to_lists m =
  let open Batteries in
  A.to_list (A.map (A.to_list % M.to_array) (M.to_rows m))

(** Create a 1xN vector from a list of floats of length N. *)
let vec_from_list ?(col=false) l =
  match col with
  | false -> M.of_array (A.of_list l) 1 (L.length l) 
  | true -> M.of_array (A.of_list l) (L.length l)  1

(** Create a 1xN vector from a list of integers of length N. *)
let vec_from_int_list ?(col=false) l =
  vec_from_list ~col (L.map float l)

(** Throw an exception if the length of list l is not equal to cols. *)
let check_length cols l = 
  if cols <> L.length l
  then raise (Failure "input rows have different lengths")

(* This CAN BE REVISED to use of_array. *)
(** Create an MxN matrix from a list of M lists of N floats. 
    Throws an exception if the internal lists have different lengths. *)
let mat_from_lists ls =
  match ls with
  | [] -> M.empty 0 0
  | l' :: ls' -> 
      let rows, cols = L.length ls, L.length l' in
      L.iter (check_length cols) ls';
      let mat = M.empty rows cols in
      let fill_row i l = L.iteri (fun j e -> M.set mat i j e) l in
      L.iteri fill_row ls;
      mat

(** Create an MxN matrix from a list of M lists of N integers. 
    Throws an exception if the internal lists have different lengths. *)
let mat_from_int_lists ls =
  mat_from_lists (L.map (fun l -> L.map float l) ls)


(********************************************)
(** Numerical utility functions *)

(** Calculate a pair of matrix indexes from an index into a vector and
    a row width for the matrix.  i.e. if we laid out a matrix, one row 
    after another in vector form, idx would be an index into it, and width
    woud be the row width of the original matrix.  The corresponding index
    pair is recovered by this function. *)
let flat_idx_to_rowcol width idx =
  let row = idx / width in
  let col = idx mod width in
  row, col

(** [rowcol_to_flat_idx width row col] calculates a 1D index into an array or
    ndarray in which matrix rows are laid out sequentially. *)
let rowcol_to_flat_idx width row col = (row * width) + col

(** Rounds second arg to number of decimal digits specified by first arg. *)
let roundto digits x =
  let scale = F.int_pow 10. digits in
  (F.round (scale *. x)) /. scale

let int_is_positive = ((<) 0);;
let int_is_nonnegative = ((<=) 0);;
let int_is_negative = ((>) 0);;
let int_is_nonpositive = ((>=) 0);;

let float_is_positive = ((<) 0.);;
let float_is_nonnegative = ((<=) 0.);;
let float_is_negative = ((>) 0.);;
let float_is_nonpositive = ((>=) 0.);;
