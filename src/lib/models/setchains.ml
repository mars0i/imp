(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998.
    Please see this book for definitions of terms, proofs, algorithms. *)

module L = Batteries.List
module A = Batteries.Array
module M = Owl.Mat
module Pmap = Parmap

module U = Utils.Genl
module WF = Wrightfisher
module T = Tdists

(** One goal here is to create a "distlist", which is a LazyList of Lists 
    Owl row vector matrices representing probability distributions over 
    possible frequencies of alleles (or other organism types0 in a population. 
    The basic distlist creation functions are in tranmats.ml, but this module
    uses higher-level functions in Wrightfisher to create distlists. *)

(************************************************************)
(** Utility helper functions *)

let always _ = true

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

(************************************************************)

(** Reusable sanity check.  The arguments p and q should be Owl vectors,
    i.e. 1 x  n or n x 1 matrices. *)
let sanity_check_vec_interval p q =
  let p_shape = M.shape p in
  if p_shape <> (M.shape q) then raise (Failure "Vectors aren't same size");
  if (M.sum' p) > 1. then raise (Failure "Low vector sums to > 1");
  if (M.sum' q) < 1. then raise (Failure "High vector sums to < 1");
  if M.exists U.float_is_negative M.(q - p) then raise (Failure "High vector is not >= p vector everywhere");
  () (* redundant clarification *)

(************************************************************)
(** Tight Interval Algorithm from p. 31: *)

(** Return a tightened version of the value at index idx in this_vec
    given other_vec.  *)
let tighten_one_coord relation idx this_vec other_vec =
  let this_elt = M.get this_vec 0 idx in
  let other_sum = (M.sum' other_vec) -. (M.get other_vec 0 idx) in
  if relation (other_sum +. this_elt) 1. then this_elt
  else 1. -. other_sum

(* This can probably be sped up in the way that recombine was using 
   Evik Tak's suggestion: https://stackoverflow.com/a/46127060/1455243 *)
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
let tighten_vec_interval p q =
  sanity_check_vec_interval p q;
  let p' = tighten_vec (>=) p q in
  let q' = tighten_vec (<=) q p in
  (p', q')

(** Given a vector interval, i.e. a list containing a lower and upper vector,
    return a tight vector interval, i.e. a list containing an upper and a lower 
    vector.  *)
let tighten_vec_interval2 pq =
  let p, q = pq in
  tighten_vec_interval p q

(* This is very slow, taking time on the order of 2^(N/100) seconds on my MBA.
 * Almost all of the time is in tighten_vec; the splitting and concatenation
 * have negligible impact. *)
(** Matrix interval tightener *)
let tighten_mat_interval low high =
  (* sanity_check_vec_interval low high; *) (* needs to be different for nxn matrices *)
  let low_rows = M.to_rows low in
  let high_rows = M.to_rows high in
  let low_tight_vecs  = A.map2 (tighten_vec (>=)) low_rows high_rows in
  let high_tight_vecs = A.map2 (tighten_vec (<=)) high_rows low_rows in
  let low'  = M.concatenate low_tight_vecs in
  let high' = M.concatenate high_tight_vecs in
  (low', high')

  
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
  L.map U.list_to_vec (list_vertices ?digits ?uniq p q)

let vec2vec_vertices ?digits ?uniq p q =
  vec_vertices ?digits ?uniq (U.vec_to_list p) (U.vec_to_list q)

(** Given a min and max matrices p and q for a tight interval, return a list of
    vertex matrices.  See documentation for list_vertices for additional info *)
let mat_vertices ?digits ?uniq p q =
  let p_rows, q_rows = A.to_list (M.to_rows p), A.to_list (M.to_rows q) in   (* lists of the row vectors from p and q *)
  let _ = L.map2 sanity_check_vec_interval p_rows q_rows in
  let vec_verts = L.map2 (vec2vec_vertices ?digits ?uniq) p_rows q_rows in   (* A list of lists of vectors. Each list reps row vertices for one row *)
  let vec_vert_arrays = L.map A.of_list (L.n_cartesian_product vec_verts) in (* list of (ordered) arrays of vectors rep'ing rows of vertex matrices *)
  L.map M.of_rows vec_vert_arrays


(************************************************************)
(** Hi-Lo Method.  See Hartfiel section 2.4, pp.  46-54

This is an analog of matrix multiplication for intervals of stochastic
matrices.  The method calculates tight bounds for products of all
possible multiplications of matrices in two intervals, the original
interval, and one that's the result of a previous application of the
hi-lo method.  The core idea of this process is that the elements of a
new pair of lower and upper matrices must be individually calculated so
that each element estimates the minimum or maximum value for that
element from all of the possible products.

Utop usage example using [CredalsetIO.make_setchain_bounds_pdfs]:
{[
Module I = Models.CredalsetIO;;
Module S = Models.Setchains;;
Module T = Models.Tdists;;
Module W = Models.Wrightfisher;;
let fitns = W.([{w11=1.0; w12=0.5; w22=0.2}; {w11=1.0; w12=0.8; w22=0.1}]);;
let pmat, qmat = S.make_wf_interval 100 fitns;;
let bounds_mats =  S.lazy_bounds_mats_list pmat qmat;;
let intervals = S.lazy_prob_intervals_from_freq 50 bounds_mats;;
let tdistlists = T.add_gens intervals;;
I.make_setchain_bounds_pdfs ~rows:2 ~cols:3 "foo"
                            (T.sublist 1 12 tdistlists);;
]}
More variations are illustrated in setchainPDFs.ml. *)

(** Compare function for use by idx_sort for col vector *)
let col_vec_idx_cmp mat i i' =
  if M.get mat i 0 > M.get mat i' 0 then 1 
  else if M.get mat i 0 < M.get mat i' 0 then -1 
  else 0

(** Given column vector, return a list of indexes in
    order of the numerical order of the values at those indexes. *)
let idx_sort_colvec v =
  let size, _ = M.shape v in
  let idxs = L.range 0 `To (size - 1) in
  L.fast_sort (col_vec_idx_cmp v) idxs

(** "Recombination" functions (by analogy with genetic recombination) that
    take two vectors and create a new vector from parts of each of them,
    though in this case the order in which the elements are considered is
    not the linear order within vectors. 
    (Don't forget to tighten the arguments first.) 
    NOTE: Functions like this one that take a function argument such as 
    [relation] or [recomb] are mainly intended to be used for building 
    other functions.  Among other things, the functionals will usually 
    require arguments in a different order depending on which function 
    is passed.  *)

(* This version of recombine uses a suggestion by Evik Tak: https://stackoverflow.com/a/46127060/1455243 *)
(** This function is at the core of the hi-lo method.
    Given a relation (>=), a column l vec and two tight row vecs p and q s.t. 
    p<=q, return a stochastic row vec ("p bar") with high values from q where l
    is low and low values from p where l is high.  Or pass (<=), l, and tight
    row vecs s.t. p >= q to return a stoch row vec ("q bar") with low values 
    from p where l is low.  Note that the latter swaps the normal meanings of 
    p and q in Hartfiel, i.e. here the arguments should be (<=), l, q, p
    according to the normal senses of p and q. *)
let recombine relation p q p_sum idxs =
  (* This was a mistake because sometimes both p and q are nothing but a single 1.0 entry with zeros: *)
  (* if p = q then raise (Failure "recombine: Lower and upper tran matrices are identical"); *)
  let pbar = M.copy p in  (* p was created using M.row, so it's a view not a copy. *)
  let rec find_crossover idxs' psum =
    match idxs' with
    | i::idxs'' -> 
        let qi = M.get q 0 i in
        let sum_rest = psum -. (M.get pbar 0 i) in (* pbar begins <= 1 if p<=q, or >= 1 if p, q swapped *)
        let sum_rest_plus_qi = (sum_rest +. qi) in
        if relation sum_rest_plus_qi 1.
        then M.set pbar 0 i (1. -. sum_rest) (* return--last iter put it over/under *)
        else (M.set pbar 0 i qi;             (* still <= 1, or >=1; try next one *)
              find_crossover idxs'' sum_rest_plus_qi) 
    | [] -> raise (Failure "recombine: bad vectors") (* this should never happen *)
  in 
  find_crossover idxs p_sum;
  pbar

(** Given a [recomb] function, the original P and Q matrices [pmat] and [qmat],
    a previous tight bounds matrix [prev_bound_mat], pre-calculated row sums
    [pmat_row_sums] from [pmat] (more efficient to calculate once), and an 
    array of sorted lists of indexes [prev_mat_idx_lists], calculate the value
    at i j for the next tight bounds matrix, where i j is calculated from vector
    index [idx] and [width].  i.e. if we laid out a matrix one row after another
    in vector form, idx would be an index into it, and width is the row width
    of the original matrix.  (We pass [pmat_row_sums] and [idx_lists] even 
    though they could be calculated on demand from [pmat], [prev_bound_mat],
    to avoid repeatedly performing the same computations.) Used by [_hilo_mult].  
    SEE doc/nonoptimizedcode.ml for an older, perhaps clearer version.  *)
let calc_bound_val recomb pmat qmat prev_bound_mat pmat_row_sums prev_mat_idx_lists width idx =
  let i, j = U.flat_idx_to_rowcol width idx in
  let p_row_sum = M.get pmat_row_sums i 0 in
  let idxs = A.get prev_mat_idx_lists j in
  let p_row, q_row = M.row pmat i, M.row qmat i in (* row doesn't copy; it just provides a view *)
  let bar_row = recomb p_row q_row p_row_sum idxs in
  let prev_col = M.col prev_bound_mat j in (* col makes a copy *) (* TODO CHECK HAS THIS CHANGED IN OWL *)
  M.(get (bar_row *@ prev_col) 0 0)

(** Wrapper for calc_bound_val (which see), adding an additional, ignored 
    argument for Pmap.array_float_parmapi. *)
let calc_bound_val_for_parmapi recomb pmat qmat prev_bound_mat pmat_row_sums prev_mat_idx_lists width idx _ =
  calc_bound_val recomb pmat qmat prev_bound_mat pmat_row_sums prev_mat_idx_lists width idx

(** Given a [recomb] function, the original tight interval bounds pmat and
    qmat, and either the previous tight component lo or hi bound 
    [prev_bound_mat](as appropriate), and a vector [row_sums] of summed rows
    (so we don't have to do it again), return the next lo or hi tight component 
    bound.  This function normally uses [Parmap] to split the work between 
    additional cores.  If [~fork] is false with any value, it won't divide 
    the work between processes.
    NOTE:
      If recomb is [recombine (>=)], the arguments should be [pmat], [qmat], 
      and the previous lo matrix, along with the row sums.  
      If recomb is [recombine (<=)], the arguments should be (notice!) [qmat],
      [pmat], and the previous hi matrix, along with the row sums. *)
let _hilo_mult ?(fork=true) recomb pmat qmat prev_bound_mat row_sums = 
  let (rows, cols) = M.shape pmat in
  let len = rows * cols in
  let prev_mat_idx_lists = M.map_cols idx_sort_colvec prev_bound_mat in (* sorted list of indexes for each column *)
    if fork 
    then let bounds_array' = A.create_float len in
      M.of_array
        (Pmap.array_float_parmapi (* ~ncores:4 *)
          ~result:bounds_array' 
          (calc_bound_val_for_parmapi recomb pmat qmat prev_bound_mat row_sums prev_mat_idx_lists cols)
          bounds_array') (* this arg will be ignored *)
        rows cols
    else M.init rows cols (calc_bound_val recomb pmat qmat prev_bound_mat row_sums prev_mat_idx_lists cols)

(** Starting from the original [pmat] and [qmat] tight interval bounds, the 
    previous component tight lo bound, and a vector of p row sums [p_row_sums],
    make the next lo matrix.
    If [~fork] is present with any value, won't use Parmap to divide the 
    work between processes. *)
let _lo_mult ?(fork=true) pmat qmat prev_lo_mat p_row_sums =
  _hilo_mult ~fork (recombine (>=)) pmat qmat prev_lo_mat p_row_sums

(** Starting from the original [pmat] and [qmat] tight interval bounds, the 
    previous component tight hi bound, and a vector of q row sums [q_row_sums],
    make the next hi matrix.
    NOTE args are in same order as lo_mult.
    If [~fork] is present with any value, won't use Parmap to divide the 
    work between processes. *)
let _hi_mult ?(fork=true) pmat qmat prev_hi_mat q_row_sums =
  _hilo_mult ~fork (recombine (<=)) qmat pmat prev_hi_mat q_row_sums (* NOTE SWAPPED ARGS *)

(** Tip: The next few functions create a LazyList in which each element is
    constructed from the preceding one by a method that usually forks 
    multiple operating system processes.  If you abort the processing
    before it completes, you may end up with a partially constructed or
    otherwise somehow corrupted element in the LazyList.  Since a LazyList
    won't recalculate an element once it's been created the first time,
    your lazy list may become useless, and you'll have to regenerate it
    from scratch.  (That's my interpretation of something that happened once.) *)

(** Return pair of pairs: The first pair is the bounds matrices that were
    passed as the third argument [(lo,hi)], unchanged, and the next bounds
    matrix pair.  For use with [Batteries.LazyList.from_loop] *)
let next_bounds_mats ?(fork=true) pmat qmat p_row_sums q_row_sums (lo,hi) =
  let lo' = _lo_mult ~fork pmat qmat lo p_row_sums in
  let hi' = _hi_mult ~fork pmat qmat hi q_row_sums in
  (lo', hi')

(** lazy_bounds_mats [pmat] [qmat] returns a LazyList of bounds matrix pairs
    starting from the initial transition matrix interval defined [pmat] defined
    by [qmat] *)
let lazy_bounds_mats_list ?(fork=true) pmat qmat =
  let p_row_sums, q_row_sums = M.sum_cols pmat, M.sum_cols qmat in (* sum_cols means add col vecs, = new col vec w/ sum of ea row *)
  T.seq (pmat, qmat) (next_bounds_mats ~fork pmat qmat p_row_sums q_row_sums) always

(** Convenience version of lazy_bounds_mats_list that takes
    (pmat, qmat) as argument rather than pmat and qmat. *)
let lazy_bounds_mats_list_from_pair ?(fork=true) (pmat, qmat) =
  lazy_bounds_mats_list ~fork pmat qmat

(** Given a transition matrix interval [(lo_mat, hi_mat)] and a probability 
   interval that's represented by a single frequency [freq] to which all 
   probability is assigned for the single distribution in the interval, 
   lo-multiply hi-multiply the distribution times [lo_mat] and [hi_mat]
   respectively.  (When the probabilty interval is of this kind, this function
   should be more efficient than lo_mult and hi_mult, and even more efficient
   than the normal dot product, which does the same thing when the interval
   contains only one distribution.) *)
let freq_mult freq (lo_mat, hi_mat) =
  [M.row lo_mat freq; M.row hi_mat freq]

(** lazy_prob_intervals_from_freq [freq] [bounds_mats_list] expects an initial
    frequency for a single population and a LazyList of bounds matrix pairs,
    and returns a LazyList of probability intervals for each timestep. Like
    [lazy_prob_intervals] and [lazy_singleton_intervals] but more efficient
    if the probability interval consists of a single distribution that puts
    all probability on a single frequency.  Note that the first element will
    be the initial interval: a list containing two copies of a vector with
    1.0 at index freq and 0.0 everywhere else. *)
let lazy_prob_intervals_from_freq freq bounds_mats_list =
  let size, _ = M.shape (fst (T.hd bounds_mats_list)) in
  let init_dist = (WF.make_init_dist size freq) in
  T.cons [init_dist; init_dist] (T.map (freq_mult freq) bounds_mats_list)


(*****************************************************)
(** Functions for making matrix intervals *)

(** The next two functions make a matrix interval that will be big enough
    (but no bigger) to contain the probability distributions defined by applying a
    Wright-Fisher model to all of the the specified fitness triples.  Note
    that there is no reason that intervals must be defined in this way, even
    if you're thinking in terms of Wright-Fisher processes. *)

(** Make an interval from a popsize and list of fitness structures without
    verifying tightness. Raises an exception if the resulting lower and upper 
    matrices are identical.  (After tightening, such matrices will cause
    an error in [recombine].)
    *)
let make_wf_interval_no_tight_check popsize fitn_list =
  let tranmats = L.map (WF.make_tranmat popsize) fitn_list in
  let low = L.reduce M.min2 tranmats in
  let high = L.reduce M.max2 tranmats in
  (* Test for degenerate cases.  (Here rather than testing fitnesses, since 
     it's possible for the same probability relationships to be generated by 
     different fitness triples. *)
  if low = high then raise (Failure "make_wf_interval_no_tight_check: Lower and upper tran matrices are identical");
  low, high

(** Make an interval from a popsize and list of fitness structures, making
    sure that it is tight.  Intervals made this way should always be tight, 
    so this is just sanity check.  For large population sizes it might be 
    it might be significantly faster to call the version without the
    tightness test directly. *)
let make_wf_interval popsize fitn_list =
  let low, high = make_wf_interval_no_tight_check popsize fitn_list in
  let tight_low, tight_high = tighten_mat_interval low high in
  if (low, high) <> (tight_low, tight_high) (* This should not happen; such an interval should already be tight. *)
  then Printf.eprintf "\n[make_wf_interval] Note: had to tighten original Wright-Fisher-based interval.\n";
  tight_low, tight_high


(*
module T = Tdists

let lazy_tdists_from_freq ?(first_tick=0) freq bounds_mats_list =
  let f dists tdists_list = 
    let prev_gen = T.((T.hd tdists_list).gen) in
    T.{gen = prev_gen + 1; dists}
  in
  let lazy_intervals = lazy_prob_intervals_from_freq freq bounds_mats_list in
  let first_interval = T.hd lazy_intervals in
  let rest_intervals = T.tl lazy_intervals in
  T.lazy_fold_right f rest_intervals {t = 0; dists = first_interval}
*)

(*
regular list illustration:

type x = {a:int; b:int};;

fold_right (fun x ys -> let {a; b} = hd ys in {a=(a+1); b=x}::ys) z [{a=0;b=0}];;
- : x list =
[{a = 12; b = 0}; {a = 11; b = 1}; {a = 10; b = 2}; {a = 9; b = 3}; {a = 8; b = 4};
 {a = 7; b = 5}; {a = 6; b = 6}; {a = 5; b = 7}; {a = 4; b = 8}; {a = 3; b = 9};
 {a = 2; b = 10}; {a = 1; b = 11}; {a = 0; b = 0}]

*)


(** USAGE EXAMPLES: 
     module S = Models.Setchains;;
     module T = Models.Tranmats;;
     module W = Models.Wrightfisher;;
     let pmat, qmat = W.(S.make_wf_interval 100 [{w11=1.0; w12=0.5; w22=0.3}; {w11=0.2; w12=0.9; w22=1.0}]);;
     let bounds_mats = S.lazy_bounds_mats_list pmat qmat;;
     let distslist = S.lazy_prob_intervals_from_freq 50 bounds_mats;;
     let tdistslist = T.add_gens distslist;;
   or:
     let tdistslist = 
       T.add_gens (S.lazy_prob_intervals_from_freq 50
                   (S.lazy_bounds_mats_list_from_pair
                       W.(S.make_wf_interval 100 [{w11=1.0; w12=0.5; w22=0.3};
                                                  {w11=0.2; w12=0.9; w22=1.0}])));;
     (* however at present you need the one without ts, too *)

   Then to pass a sublist e.g. to make_setchain_bounds_pdfs use
     T.sublist 1 4 tdistlists
   or
     T.sublist 1 4 tdistlists
   e.g. like this:
     I.make_setchain_bounds_pdfs ~rows:2 ~cols:2 "yo" (T.sublist 1 4 tdistlists);;

*)
