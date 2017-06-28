(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998.
    Please see this book for definitions of terms, proofs, algorithms. *)

module L = Batteries.List
module A = Batteries.Array
module M = Owl.M

(** Given two lists of length n, return a list containing the 2^n lists
    containing each combination of elements from p and q at the same indices.
    e.g.
       vertices [1; 2; 3] [100; 200; 300]
    returns
       [[1; 2; 3];   [1; 2; 300];   [1; 200; 3];   [1; 200; 300];
        [100; 2; 3]; [100; 2; 300]; [100; 200; 3]; [100; 200; 300]] *)
let rec sequences p q =
  match p, q with
  | [], [] -> []
  | hp::[], hq::[] -> [[hp];[hq]]
  | hp::tp, hq::tq -> 
      let tailverts = vertices tp tq in
      L.concat [L.map (L.cons hp) tailverts; L.map (L.cons hq) tailverts]
  | _, _ -> raise (Failure "lists are not the same length")

(** Converts a 1xN matrix, i.e. row vector, into a list. *)
let vec_to_list m = A.to_list (M.to_array m)

(** Converts an MxN matrix into a list of M lists of length N. *)
let mat_to_lists m = A.to_list (A.map (A.to_list % M.to_array) (M.to_rows m))
  
(** Given the vector boundaries of an interval, lists its vertices assuming 
    it is tight. *)
let vertices p q =
  let (_, size) = M.shape p in
  let pl, ql = mat_to_lists p, mat_to_lists q in
  let vertices_at i pl ql =
    let pl', ql' = L.remove_at i pl, L.remove_at ql in
    let seqs = sequences pl' ql' in
    let sums = L.map (fun seq -> 1 - (L.fsum seq)) seqs in
    (* now take the ones that are nonnegative, and use them to select
     * the sequences that can be made stochastic, and then insert the
     * nonnegative value into the ith location in the sequence, and
     * ultimately make that into a vector. *)


(** Can also be used with vectors. *)
let mat_vertices m1 m2 =
  let m1_list, m2_list = mat_to_lists m1, mat_to_lists m2 in
  L.map2 vertices m1_list m2_list




