
(* module M = Owl.Mat *)

(* Tips
  * In po_examples.ml and hasse.ml, the keys are the things that get
  * partially ordered.  These models just unit for *all* of the values.
  * Note that type el is the *key* type.
  * 
  * When you make an empty pomap, it has type 'a pomap, where 'a is the
  * value type.  When you add an element to it and return a new pomap
  * object, it will have the type of whatever kind of value you used.
  *)


(* notes *)
#require "pomap";;
module P = Pomap_impl.Make(struct type el = int type ord = Unknown | Lower | Equal | Greater let compare x y = if x < 0 || y < 0 then Unknown else if x = y then Equal else if x < y then Lower else Greater end);;
let pm = P.empty;;
pm;;
P.add 2 3 pm;;
pm;;
P.fold (fun n acc -> (P.get_key n, P.get_el n)::acc) pm [] ;;
P.topo_fold (fun n acc -> (P.get_key n, P.get_el n)::acc) pm [] ;;
P.find 2 pm;;
let pm = P.add 2 3 pm;;
P.find 2 pm;;
open P;;
P.find 2 pm;;
P.get_el (P.find 2 pm);;
let (x,y) = P.find 2 pm;;
let (x,y) = P.find 2 pm in P.get_el y;;

type ord_type = Unknown | Lower | Equal | Greater

let simple_compare x y = 
  if x = y then Equal
  else if x < y then Lower
  else if x > y then Greater 
  else Unknown

let efficient_compare m1 m2 =
  let f acc e1 e2 =
    match acc with
    | Equal -> if e1 < e2 then Lower else 
               if e1 > e2 then Greater else Equal
    | Lower ->   if e1 <= e2 then Lower else Unknown
    | Greater -> if e1 >= e2 then Greater else Unknown
    | Unknown -> Unknown (* should never happen given short circuit *)
  in Utils.short_circuit_fold2 Unknown f Equal m1 m2


module MatPO = struct 
  type el = Owl.Mat.mat
  type ord = ord_type
  let compare = efficient_compare
end

module P = Pomap_impl.Make(MatPO)

let add_mat m pm = P.add m () pm

let pm = 
  P.singleton (mat_from_int_lists [[4; 4];[4; 4]]) () |>
  add_mat (mat_from_int_lists [[3;3];[3;3]]) |>
  add_mat (mat_from_int_lists [[2;3];[3;3]]) |>
  add_mat (mat_from_int_lists [[2;2];[3;3]]) |>
  add_mat (mat_from_int_lists [[2;2];[3;1]]) |>
  add_mat (mat_from_int_lists [[2;2];[3;2]]) |>
  add_mat (mat_from_int_lists [[2;1];[3;2]]) |>
  add_mat (mat_from_int_lists [[1;1];[3;2]]) |>
  add_mat (mat_from_int_lists [[3;3];[1;3]]) |>
  add_mat (mat_from_int_lists [[3;3];[3;5]]) |>
  add_mat (mat_from_int_lists [[3;3];[4;5]]) |>
  add_mat (mat_from_int_lists [[4;3];[4;5]]) |>
  add_mat (mat_from_int_lists [[4;4];[4;5]]) |>
  add_mat (mat_from_int_lists [[4;6];[4;5]]) |>
  add_mat (mat_from_int_lists [[4;6];[4;6]]) |>
  add_mat (mat_from_int_lists [[4;6];[4;7]]) |>
  add_mat (mat_from_int_lists [[4;6];[7;7]]) |>
  add_mat (mat_from_int_lists [[8;6];[7;7]])

