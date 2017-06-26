(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998 *)

module L = Batteries.List
module A = Batteries.Array

let vertices_wrt_i i p q =
  let size = A.length p in
  if size != A.length q then raise (Failure "vectors have different lengths")
  else if i >= size || i < 0 then raise (Failure "index negative or too large")
  else let pl, ql = A.to_list p, A.to_list q in
  let p_unfree, q_unfree = L.remove_at i pl, L.remove_at i ql in
  (* TODO here generate a list of lists with all possible combinations of
   * one element from either p or q at each unfree index *)
  (* Then for each such vector/list, compute 1 - (sum of values in it).
   * if this value is >= 0, then add to the output list the vector
   * containing this value inserted into the vector/list at the
   * missing position.  If the new value is < 0, then simply move on. *)


let vertices p q =
  (* concatenate the results of vertices_wrt_i for each i *)





(*
let vertices p q =
  let idxs = L.range 0 `To (L.length p) in
  let verts free =
    let unfree = L.remove idxs free in
    L.concat (L.map (fun i -> [L.at p i; L.at q i]) unfree

  L.concat 
  (L.init
*)


