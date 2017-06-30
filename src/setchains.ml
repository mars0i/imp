(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998.
    Please see this book for definitions of terms, proofs, algorithms. *)

module B = Batteries
module L = Batteries.List
module A = Batteries.Array
module M = Owl.Mat
module U = Utils

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
      let tailverts = sequences tp tq in
      L.concat [L.map (L.cons hp) tailverts; L.map (L.cons hq) tailverts]
  | _, _ -> raise (Failure "lists are not the same length")

(** Converts a 1xN matrix, i.e. row vector, into a list. *)
let vec_to_list = B.(A.to_list % M.to_array)
let list_to_vec = B.(M.of_array % A.of_list)

(** Converts an MxN matrix into a list of M lists of length N. *)
let mat_to_lists m =
  let open Batteries in
  A.to_list (A.map (A.to_list % M.to_array) (M.to_rows m))

let slop = 0.0000001

(** List vertices in which the ith element is free *)
let vertices_at p q i =
  let min_at, max_at = (L.at p i) -. slop, (L.at q i) +. slop in
  let p', q' = L.remove_at i p, L.remove_at i q in
  let seqs = sequences p' q' in
  let sums = L.map (fun seq -> 1. -. (L.fsum seq)) seqs in
  let add_vertex sum seq acc =
    if sum >= min_at && sum <= max_at  (* assume min_at >= 0 *)
    then (U.insert_before i sum seq)::acc else acc
  in L.fold_right2 add_vertex sums seqs []

(* TODO Make me more efficient? *)
(** Given the vector boundaries of an interval, lists its vertices assuming 
    it is tight.  Might return two variants of the same vertex that
    differ only by float rounding errors.  If ~digits:digs is provided, 
    digs will specify number of decimal digits to round to.  If ~uniq:true
    is provided, duplicate vertices are combined.  Doing this usually
    makes sense only if ~digits is also provided. *)
let list_vertices ?(digits) ?uniq p q =
  let idxs = L.range 0 `To ((L.length p) - 1) in (* kludgey *)
  let verts = L.concat (L.map (vertices_at p q) idxs) in
  let verts' = match digits with
               | None -> verts
               | Some digs -> U.mapmap (U.roundto digs) verts
  in match uniq with
  | None | Some false -> verts'
  | Some true -> L.unique_cmp verts'

(** Convenience alias for list_vertices ~digits:3 ~uniq:true *)
let verts3 = list_vertices ~digits:3 ~uniq:true

(* Alternatives to Batteries.List.unique_cmp:
 * Batteries.List.unique (slower),
 * Core.List.dedup
 * Core.List.stable_dedup *)


(* FIXME broken *)
let mat_vertices p q = 
 L.map list_to_vec (list_vertices (vec_to_list p) (vec_to_list q))
