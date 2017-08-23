(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998.
    Please see this book for definitions of terms, proofs, algorithms. *)

module B = Batteries
module L = Batteries.List
module A = Batteries.Array
module M = Owl.Mat
module U = Utils
module Pm = Probmat

(************************************************************)
(** Utility helper functions *)

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


(************************************************************)
(** Tight Interval Algorithm from p. 31: *)

(** Return a tightened version of the value at index idx in this_vec
    given other_vec.  *)
let tighten_one_coord relation idx this_vec other_vec =
  let this_elt = M.get this_vec 0 idx in
  let f _ i acc elt =  (* sum other_vec only at other idxes *)
    if i = idx then acc else acc +. elt in
  let other_sum = M.foldi f 0. other_vec in
  if relation (other_sum +. this_elt) 1. then this_elt
  else 1. -. other_sum

(** Returns a tight version of [this_vec], which could be either the lower
    or upper vector of the interval.  You must also pass [other_vec], i.e.
    whichever end vector of the interval is not this_vec. Relation should be
    (>=) if this_vec is the lower vector, and (<=) if it's the upper vector. *)
let tighten_vec relation this_vec other_vec =
  let _, n = M.shape this_vec in
  let tight_vec = M.empty 1 n in
  for i = 0 to (n - 1) do
    M.set tight_vec 0 i (tighten_one_coord relation i this_vec other_vec)
  done;
  tight_vec

(** Given a lower and upper vector, return a tight vector interval, i.e. a list
    containing an upper and a lower vector.  *)
let tighten_interval2 p q =
  let p' = tighten_vec (>=) p q in
  let q' = tighten_vec (<=) q p in
  [p'; q']

(** Given a vector interval, i.e. a list containing a lower and upper vector,
    return a tight vector interval, i.e. a list containing an upper and a lower 
    vector.  *)
let tighten_interval pq =
  match pq with
  | p::q::[] -> tighten_interval2 p q
  | _ -> raise (Failure "Not an interval.")

(* TODO add matrix interval tightener. *)
(* Need to write map_rows2 *)

  
(************************************************************)
(** Determine vertices: *)

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


(************************************************************)
(* Find vertices of a list of 2D stochastic vectors *)
(* TODO
let twoD_vertices vs = calculate min and max of first coord, or mins of each coord.  *)

(************************************************************)
(** Hi-Lo Method *)
(* in progress *)

(** Use separate compare functions for row and column vectors to avoid
    having a test for row vs. col inside the compare function. *)

(** Compare function for use by idx_sort for row vector *)
let row_vec_idx_cmp mat j j' =
  if M.get mat 0 j > M.get mat 0 j' then 1 
  else if M.get mat 0 j < M.get mat 0 j' then -1 
  else 0

(** Compare function for use by idx_sort for col vector *)
let col_vec_idx_cmp mat i i' =
  if M.get mat i 0 > M.get mat i' 0 then 1 
  else if M.get mat i 0 < M.get mat i' 0 then -1 
  else 0

(** Given a row or column vector, return a list of indexes in
    order of the numerical order of the values at those indexes.
    Automatically determines whether the argument is a row or a column 
    vector, raising an exception if neither. *)
let idx_sort v =
  let rows, cols = M.shape v in
  if rows > 1 && cols > 1 then raise (Failure "Argument is not a vector.");
  let size, idx_cmp = if rows = 1
                      then cols, row_vec_idx_cmp 
                      else rows, col_vec_idx_cmp in
  let idxs = A.of_list (L.range 0 `To (size - 1)) in
  A.fast_sort (idx_cmp v) idxs;
  idxs

(** Return the sum of all values in matrix except the one at i j. *)
let sum_except mat i j = 
  (M.sum mat) -. (M.get mat i j)

(* TODO generalize so that it can make qbar, too. *)
(** Given a column l vec and two rows vecs p and q, return a new pbar vector
    with high values from q where l is low and low values from p where l is
    high. *)
let recombine_p l p q =
  let size, _ = M.shape l in
  let idxs = idx_sort l in
  let pbar = M.clone p in  (* Note sum p should always ust be <= 0. *)
  let j = ref (-1) in
  let not_found = ref true in
  while !not_found do  (* forgive me for my sin--thought I'd try it this way *)
    incr j; if !j >= size then raise (Failure "bad vectors");
    M.print pbar;
    print_int !j;
    let i = idxs.(!j) in  (* we're not walking through pbar in its index order *)
    let sum_rest = sum_except pbar 0 i in
    let qi = M.get q 0 i in
    if  qi +. sum_rest <= 1. then M.set pbar 0 i qi (* FIXME logic wrong here. need to test also for >= 1. *)
    else (M.set pbar 0 i (1. -. sum_rest);
          not_found := true)
  done;
  pbar
(* FIXME to see the error, try these values:
let l = Owl.Mat.of_array [|0.3; 0.5; 0.4|] 3 1;;
let p = Owl.Mat.of_array [|0.1; 0.2; 0.4|] 1 3;;
let q = Owl.Mat.of_array [|0.4; 0.3; 0.6|] 1 3;;
  *)

(************************************************************)
(** Example 2.10 *)

let mp  = M.of_array [|0.35; 0.55; 0.25; 0.65|] 2 2;;
let mq  = M.of_array [|0.45; 0.65; 0.35; 0.75|] 2 2;;
let s0p = M.of_array [|0.4; 0.5|] 1 2;;
let s0q = M.of_array [|0.5; 0.6|] 1 2;;
let m = [mp; mq];;
let s0 = [s0p; s0q];;

let m_verts = mat_vertices ~uniq:true ~digits:3 mp mq;;
let s0_verts = mat_vertices ~uniq:true ~digits:3 s0p s0q;;

let w1 = Pm.cross_apply M.dot s0_verts m_verts;;  (* dot = ( *@ )  *)
let s1_verts = [(L.last w1); (L.first w1)];;  (* (0.29, 0.71), (0.4, 0.6); see Hartfiel *)

let w2 = Pm.cross_apply M.dot s1_verts m_verts;;


(* OBSOLETE:

(** Compare function for use by idx_sort (q.v.) *)
let vec_idx_cmp ?(rowvec=true) mat j j' =
  (* swap indexes if this is a column vector *)
  let i, i', j, j' = if rowvec
                     then 0, 0, j, j'
                     else j, j', 0, 0 in
  if M.get mat i j > M.get mat i' j' then 1 
  else if M.get mat i j < M.get mat i' j' then -1 
  else 0

(** Given a row or column vector, return a list of indexes in
    order of the numerical order of the values at those indexes.
    Automatically determines whether the argument is a row or a column 
    vector, raising an exception if neither. *)
let idx_sort v =
  let rows, cols = M.shape v in
  if rows > 1 && cols > 1 
  then raise (Failure "Matrix argument is not a vector.");
  let rowvec = (rows = 1) in 
  let size = if rowvec then cols else rows in
  let idxs = L.range 0 `To (size - 1) in
    L.fast_sort (vec_idx_cmp ~rowvec v) idxs

*)

(* OBSOLETE AND BROKEN 
(* For convenience, operate on lists.
 * Is there a more efficient algorithm?? Without summing repeatedly? *)
let recombine_columns p q =
  if L.fsum p > 1.0 then raise (Failure "p sums to > 1");
  let rec construct_vec p' q' =
    match p', q' with
    | [p1], [q1] -> 
        Printf.printf "%f, %f, %f\n" p1 q1 (L.fsum q');
        if (L.fsum q') >= 1. then p', q' (* FIXME not right *)
                    else raise (Failure "q is not >= 1")
    | p1::rest_p, q1::rest_q ->
        Printf.printf "%f, %f, %f\n" p1 q1 (q1 +. (L.fsum rest_p));
        if q1 +. (L.fsum rest_p) >= 1.  (* FIXME WAIT why would this be >= 1 if it's only part of the lists?? *)
        then p', q1::rest_p
        else
          let back_p, back_q = construct_vec rest_p rest_q in
          p1::back_p, p1::back_q
    | _::_, [] -> raise (Failure "Lists not same length")
    | [], _::_ -> raise (Failure "Lists not same length")
    | _, _ -> raise (Failure "Bad columns.")  (* can this happen? *)
  in construct_vec p q
*)
