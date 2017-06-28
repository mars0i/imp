(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998 *)

module L = Batteries.List
module A = Batteries.Array
module M = Owl.Mat

(** Given two lists of length n, return a list containing the 2^n lists
    containing each combination of elements from p and q at the same indices.
    e.g.
       vertices [1; 2; 3] [100; 200; 300]
    returns
       [[1; 2; 3];   [1; 2; 300];   [1; 200; 3];   [1; 200; 300];
        [100; 2; 3]; [100; 2; 300]; [100; 200; 3]; [100; 200; 300]] *)
let rec vertices p q =
  match p, q with
  | [], [] -> []
  | hp::[], hq::[] -> [[hp];[hq]]
  | hp::tp, hq::tq -> 
      let tailverts = vertices tp tq in
      L.concat [L.map (L.cons hp) tailverts; L.map (L.cons hq) tailverts]
  | _, _ -> raise (Failure "lists are not the same length")

(* No doubt there's some slightly better way to do this. *)
let mat_to_lists m = L.map A.to_list (A.to_list (A.map M.to_array (M.to_rows m)))

(** Can also be used with vectors. *)
let mat_vertices m1 m2 =
  let m1_list, m2_list = mat_to_lists m1, mat_to_lists m2 in
  L.map2 vertices m1_list m2_list




