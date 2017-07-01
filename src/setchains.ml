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

let list_to_vec l =
  let a = A.of_list l in
  M.of_array a 1 (A.length a)

(** Converts an MxN matrix into a list of M lists of length N. *)
let mat_to_lists m =
  let open Batteries in
  A.to_list (A.map (A.to_list % M.to_array) (M.to_rows m))

(** Try to avoid inequality comparisons that are incorrect due to
    inherent floating point fudgeyness.  Used in vertices_at. *)
let slop = 0.0000001

(** List vertices in which the ith element is free *)
let vertices_at p q i =
  let min_at, max_at =  (L.at p i) -. slop,  (L.at q i) +. slop in
  let p', q' = L.remove_at i p, L.remove_at i q in
  let seqs = sequences p' q' in
  let sums = L.map (fun seq -> 1. -. (L.fsum seq)) seqs in
  let add_vertex sum seq acc =
    if sum >= min_at && sum <= max_at (* assume min >= 0, but note slop above. *)
    then (U.insert_before i sum seq)::acc else acc
  in L.fold_right2 add_vertex sums seqs []

(** Given the vector boundaries of an interval, lists its vertices assuming 
    it is tight.  Might return two variants of the same vertex that
    differ only by float rounding errors.  If ~digits:digs is provided, 
    digs will specify number of decimal digits to round to.  If ~uniq:true
    is provided, duplicate vertices are combined.  Doing this usually
    makes sense only if ~digits is also provided. *)
let list_vertices ?digits ?uniq p q =
  let idxs = L.range 0 `To ((L.length p) - 1) in
  let verts = L.concat (L.map (vertices_at p q) idxs) in
  let verts' = match digits with
               | None -> verts
               | Some digs -> U.mapmap (U.roundto digs) verts
  in match uniq with
  | None | Some false -> verts'
  | Some true -> L.unique_cmp verts'

(* Alternative uniq'ers: Batteries.List.unique_cmp, Batteries.List.unique, Core.List.dedup, Core.List.stable_dedup *)

(** Convenience alias for list_vertices ~digits:3 ~uniq:true *)
let verts3 = list_vertices ~digits:3 ~uniq:true

(* Return a list of vertices in the form of Owl vectors. See documentation
 * for list_vertices for additional information, including info on optional
 * args. *)
let vec_vertices ?digits ?uniq p q =
  L.map list_to_vec (list_vertices ?digits ?uniq p q)

let vec2vec_vertices ?digits ?uniq p q =
  vec_vertices ?digits ?uniq (vec_to_list p) (vec_to_list q)

(** Given a min and max matrices p and q for a tight interval, return a list of
    vertex matrices.  See documentation for list_vertices for additional info *)
let mat_vertices ?digits ?uniq p q =
  let p_rows, q_rows = A.to_list (M.to_rows p), A.to_list (M.to_rows q) in   (* lists of the row vectors from p and q *)
  let vec_verts = L.map2 (vec2vec_vertices ?digits ?uniq) p_rows q_rows in   (* A list of lists of vectors. Each list reps row vertices for one row *)
  let vec_vert_arrays = L.map A.of_list (L.n_cartesian_product vec_verts) in (* list of (ordered) arrays of vectors rep'ing rows of vertex matrices *)
  L.map M.of_rows vec_vert_arrays
