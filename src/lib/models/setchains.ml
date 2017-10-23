(** Functions inspired by Hartfiel's _Markov Set-Chains_, Springer 1998.
    Please see this book for definitions of terms, proofs, algorithms. *)

module L = Batteries.List
module A = Batteries.Array
module LL = Batteries.LazyList
module M = Owl.Mat
module Pmap = Parmap

module U = Utils.Genl
module WF = Wrightfisher

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

(************************************************************)

(** Reusable sanity check.  The arguments p and q should be Owl vectors,
    i.e. 1 x  n or n x 1 matrices. *)
let sanity_check_vec_interval p q =
  let p_shape = M.shape p in
  if p_shape <> (M.shape q) then raise (Failure "Vectors aren't same size");
  if (M.sum p) > 1. then raise (Failure "Low vector sums to > 1");
  if (M.sum q) < 1. then raise (Failure "High vector sums to < 1");
  if M.exists U.float_is_negative M.(q - p) then raise (Failure "High vector is not >= p vector everywhere");
  () (* redundant clarification *)

(************************************************************)
(** Tight Interval Algorithm from p. 31: *)

(** Return a tightened version of the value at index idx in this_vec
    given other_vec.  *)
let tighten_one_coord relation idx this_vec other_vec =
  let this_elt = M.get this_vec 0 idx in
  let other_sum = (M.sum other_vec) -. (M.get other_vec 0 idx) in
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
(** Hi-Lo Method.  See Hartfiel section 2.4, pp.  46-54 *)

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

(* This version of recombine uses suggestion by Evik Tak: https://stackoverflow.com/a/46127060/1455243 *)
(** Given a relation (>=), a column l vec and two tight row vecs p and q s.t. 
    p<=q, return a stochastic row vec ("p bar") with high values from q where l
    is low and low values from p where l is high.  Or pass (<=), l, and tight
    row vecs s.t. p >= q to return a stoch row vec ("q bar") with low values 
    from p where l is low.  Note that the latter swaps the normal meanings of 
    p and q in Hartfiel, i.e. here the arguments should be (<=), l, q, p
    according to the normal senses of p and q. *)
let recombine relation p q p_sum idxs =
  let pbar = M.clone p in  (* p was created using M.row, so it's a view not a copy. *)
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
    | [] -> raise (Failure "bad vectors") (* this should never happen *)
  in 
  find_crossover idxs p_sum;
  pbar

(** Given tight row vecs p and q and sorted indexs, idxs, from a column vec,
    return stochastic vec lo with high values from q where l is low, and low 
    values from p where l is high.*)
let recombine_lo p q p_sum idxs = 
  recombine (>=) p q p_sum idxs

(** Given tight row vecs p and q and sorted indexs, idxs, from a column vec,
    return stochastic vec hi with high values from q where h is high, and low 
    values from p where h is low.*)
let recombine_hi p q p_sum idxs = 
  recombine (<=) q p p_sum idxs (* note swapped args *)

(** Calculate a pair of matrix indexes from an index into a vector and
    a row width for the matrix.  i.e. if we laid out a matrix, one row 
    after another in vector form, idx would be an index into it, and width
    woud be the row width of the original matrix.  The corresponding index
    pair is recovered by this function. *)
let flat_idx_to_rowcol width idx =
  let row = idx / width in
  let col = idx mod width in
  row, col

(** Given a [recomb] function ([recombine_lo] or [recombine_hi]), the original
    P and Q matrices [p_mat] and [q_mat], a previous
    tight bounds matrix [prev_bound_mat], and an array of sorted lists of 
    indexes [idx_lists], calculate the value at i j for the
    next tight bounds matrix, where i j is calculated from vector index [idx]
    and [width].  i.e. if we laid out a matrix one row after another in
    vector form, idx would be an index into it, and width is the row width
    of the original matrix.  (We pass [p_mat_row_sums] and [idx_lists] even 
    though they could be calculated on demand from [p_mat], [prev_bound_mat],
    to avoid repeatedly performing the same computations.) Used by [hilo_mult].  
    SEE doc/nonoptimizedcode.ml for an older, perhaps clearer version.  *)
let calc_bound_val recomb p_mat q_mat prev_bound_mat p_mat_row_sums prev_mat_idx_lists width idx =
  let i, j = flat_idx_to_rowcol width idx in
  let prev_col = M.col prev_bound_mat j in (* col makes a copy *)
  let p_row_sum = M.get p_mat_row_sums i 0 in
  let idxs = A.get prev_mat_idx_lists j in
  let p_row, q_row = M.row p_mat i, M.row q_mat i in (* row doesn't copy; it just provides a view *)
  let bar_row = recomb p_row q_row p_row_sum idxs in
  M.(get (bar_row *@ prev_col) 0 0)


(** Wrapper for calc_bound_val (which see), adding an additional, ignored 
    argument for Pmap.array_float_parmapi. *)
let calc_bound_val_for_parmapi recomb p_mat q_mat prev_bound_mat p_mat_row_sums prev_mat_idx_lists width idx _ =
  calc_bound_val recomb p_mat q_mat prev_bound_mat p_mat_row_sums prev_mat_idx_lists width idx

(** Given [recombine_lo] or [recombine_hi], the original tight interval bounds
    P and Q, and either the previous tight component lo or hi bound (as
    appropriate), return the next lo or hi tight component bound.
    This function normally uses [Parmap] to split the work between additional
    cores.  If [~fork] is present with any value, won't use Parmap to
    divide the work between processes.
    NOTE:
      If recomb is recombine_lo, the arguments should be P, Q, and the 
      previous lo matrix.  
      If recomb is recombine_hi, the arguments should be (notice!) Q, P,
      and the previous hi matrix. 
    SEE doc/nonoptimizedcode.ml for clearer versions of this function.  *)
let hilo_mult ?(fork=true) recomb p_mat q_mat prev_bound_mat = 
  let (rows, cols) = M.shape p_mat in
  let len = rows * cols in
  let p_mat_row_sums = M.sum_cols p_mat in (* sum_cols means add all of the column vectors together, which gives you a col vector containing a sum of each row *)
  let prev_mat_idx_lists = M.map_cols idx_sort_colvec prev_bound_mat in (* sorted list of indexes for each column *)
  let bounds_array =
    if fork 
    then let bounds_array' = A.create_float len in
         Pmap.array_float_parmapi (* ~ncores:4 *)
           ~result:bounds_array' 
           (calc_bound_val_for_parmapi recomb p_mat q_mat prev_bound_mat p_mat_row_sums prev_mat_idx_lists cols)
           bounds_array' (* this arg will be ignored *)
    else A.init len (calc_bound_val recomb p_mat q_mat prev_bound_mat p_mat_row_sums prev_mat_idx_lists cols)  (* TODO why make an array here? Just make a matrix. *)
  in M.of_array bounds_array rows cols


(** Starting from the original P and Q tight interval bounds and the previous
    component tight lo bound, make the netxt lo matrix.
    If [~fork] is present with any value, won't use Parmap to divide the 
    work between processes. *)
let lo_mult ?(fork=true) p_mat q_mat prev_lo_mat =
  hilo_mult ~fork recombine_lo p_mat q_mat prev_lo_mat

(** Starting from the original P and Q tight interval bounds and the previous
    component tight hi bound, make the netxt hi matrix.
    NOTE args are in same order as lo_mult.
    If [~fork] is present with any value, won't use Parmap to divide the 
    work between processes. *)
let hi_mult ?(fork=true) p_mat q_mat prev_hi_mat =
  hilo_mult ~fork recombine_hi p_mat q_mat prev_hi_mat (* note swapped args *)

(** Given [recombine_lo] or [recombine_hi], the original tight interval bounds
    P and Q, and either the previous tight component lo or hi bound (as 
    appropriate), return the kth lo or hi tight component bound.  Note that
    Hartfiel treats the initial state as # 1, not 0, but this function is
    0-based.  i.e. starting from the initial matrices, k=1 will produce the 
    first calculated bounds matrices after the initial state; Hartfiel calls 
    this (e.g.) L_2.
    NOTE args do not need to be swapped even if [reccombine_hi] is used; they
    will be swapped internall by recomb. *)
let rec make_kth_bounds_mat_from_prev recomb p_mat q_mat prev_bound_mat k =
  if k <= 0 then prev_bound_mat
  else let bound_mat = hilo_mult recomb p_mat q_mat prev_bound_mat in
  make_kth_bounds_mat_from_prev recomb p_mat q_mat bound_mat (k - 1)

(** TODO revise to use lo_mult *)
(** Starting from the original P and Q tight interval bounds and the previous
    component tight lo bound, make the kth lo matrix after the previous one. *)
let make_kth_lo_mat_from_prev p_mat q_mat prev_lo_mat k =
  make_kth_bounds_mat_from_prev recombine_lo p_mat q_mat prev_lo_mat k

(** TODO revise to use hi_mult *)
(** Starting from the original P and Q tight interval bounds and the previous
    component tight hi bound, make the kth hi matrix after the previous one.
    NOTE args are in same order as lo_mult. *)
let make_kth_hi_mat_from_prev p_mat q_mat prev_hi_mat k =
  make_kth_bounds_mat_from_prev recombine_hi p_mat q_mat prev_hi_mat k (* note swapped args *)

(** Convenience function to make both the kth lo and hi matrices from an 
    earlier pair of bound matrices and the original matrices. 
    [make_kth_bounds_mats_from_prev p_mat q_mat prev_lo prev_hi 0] returns
    [prev_lo] and [prev_hi].  
    [make_kth_bounds_mats_from_prev p_mat q_mat prev_lo prev_hi 1] returns
    the next bounds matrices after [prev_lo] and [prev_hi], and for k=30,
    you'll get the 30th pair after [prev_lo] and [prev_hi]  If [prev_lo]
    and [prev_hi] are the initial matrices, this will be what Hartfiel
    calls L_30 and H_30. *)
let make_kth_bounds_mats_from_prev p_mat q_mat prev_lo_mat prev_hi_mat k =
  (make_kth_lo_mat_from_prev p_mat q_mat prev_lo_mat k),
  (make_kth_hi_mat_from_prev p_mat q_mat prev_hi_mat k)

(** Convenience function to make both the kth lo and hi matrices.  
    [make_kth_bounds_mats p_mat q_mat 0] returns [p_mat] and [qmat].
    Note that this [make_kth_bounds_mats p_mat q_mat 1] returns the first bounds
    matrices after p_mat, q_mat, which Hartfiel calls L_2, H_2.  Note that this 
    function can *only* be used to create the kth matrices starting from the 
    initial matrices.  If you want to make kth matrices from earlier (k-n)th
    matrices, use [make_kth_bounds_mats_from_prev].*)
let make_kth_bounds_mats p_mat q_mat k = 
  make_kth_bounds_mats_from_prev p_mat q_mat p_mat q_mat k

(** TODO needs testing *)
(** Given a transition matrix interval [(lo_mat, hi_mat)] and a probability 
   interval containing a single distribution vector, lo-multiply and hi-multiply
   the distribution times [lo_mat] and [hi_mat] respectively *)
let prob_interval_mult p q (lo_mat, hi_mat) =
  [lo_mult p q lo_mat; hi_mult p q hi_mat]

(** TODO needs testing *)
(** Given a transition matrix interval [(lo_mat, hi_mat)] and a probability 
   interval, lo-multiply and hi-multiply the distribution times [lo_mat] and
   [hi_mat] respectively.  (When the probabilty interval is of this kind, this 
   function should be more efficient than lo_mult and hi_mult since it uses the
   normal matrix dot product rather than lo_mult and hi_mult, which are 
   equivalent to dot product when the interval contains only one element.  If
   in addition the element puts all probability on one value, freq_mult should
   be even more efficient.) *)
let singleton_interval_mult p (lo_mat, hi_mat) =
  M.([p *@ lo_mat; p *@ hi_mat])

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

(** Given a transition matrix interval [p_mat], [q_mat], return the kth 
    probability interval, assuming that the initial state was an interval 
    consisting of a single distribution that put all probability on one 
    frequency [init_freq].
    (To get the kth probability distribution interval from the kth
    lo and hi mats, you lo_mult the low bound of the distribution
    interval with the lo mat, and you hi_mult the high bound of the
    distribuition with the hi mat.  However, if the interval contains
    a single probability distribution, then you can just use normal
    dot product multiplication.  *And* if the single distribution
    puts all probability on one frequency, then the effect of the dot
    product of the vector with the matrix is simply to extract the matrix
    that's indexed by that frequency.  That's what this function does.) *)
let make_kth_dist_interval_from_freq init_freq p_mat q_mat k =
  freq_mult init_freq (make_kth_bounds_mats p_mat q_mat k)

(** Return pair of pairs: The first pair is the bounds matrices that were
    passed as the third argument [(lo,hi)], unchanged, and the next bounds
    matrix pair.  For use with [Batteries.LazyList.from_loop] *)
let next_bounds_mats_for_from_loop ?(fork=true) pmat qmat (lo,hi) =
  let lo', hi' = lo_mult ~fork pmat qmat lo, hi_mult ~fork pmat qmat hi in
  (lo,hi), (lo', hi')

(** Tip: The next few functions create a LazyList in which each element is
    constructed from the preceding one by a method that usually forks 
    multiple operating system processes.  If you abort the processing
    before it completes, you may end up with a partially constructed or
    otherwise somehow corrupted element in the LazyList.  Since a LazyList
    won't recalculate an element once it's been created the first time,
    your lazy list may become useless, and you'll have to regenerate it
    from scratch.  (That's my interpretation of something that happened once.) *)

(* TODO needs testing *)
(** lazy_bounds_mats [p_mat] [q_mat] returns a LazyList of bounds matrix pairs
    starting from the initial transition matrix interval defined [pmat] defined
    by [qmat] *)
let lazy_bounds_mats_list ?(fork=true) p_mat q_mat =
  LL.from_loop (p_mat, q_mat) (next_bounds_mats_for_from_loop ~fork p_mat q_mat)

(** Note functions in credalsetsPDF.ml expect to see a *list* not
   a pair for each tick. *)

(* TODO needs testing *)
(** lazy_prob_intervals [p] [q] [bounds_mats_list] expects two vectors defining
    a probability interval, and a LazyList of bounds matrix pairs, and returns
    a LazyList of probability intervals for each timestep.  Note that the first
    element will be the initial interval [p; q]. *)
let lazy_prob_intervals p q bounds_mats_list =
  LL.cons [p; q] (LL.map (prob_interval_mult p q) bounds_mats_list)

(** lazy_singleton_intervals [p] [bounds_mats_list] expects a vector 
    representing the sole probability distribution in an intial probability
    interval, and a LazyList of bounds matrix pairs.  It returns a LazyList of
    probability intervals for each timestep. Like [lazy_prob_intervals] but 
    should be more efficient if the probability interval consists of a single
    distribution.  If that distribution puts all probability on a single 
    frequency, then [lazy_prob_intervals_from_freq] should be more efficient. 
    Note that the first element will be the initial interval [p; p]. *)
let lazy_singleton_intervals p bounds_mats_list =
  LL.cons [p; p] (LL.map (singleton_interval_mult p) bounds_mats_list)

(* TODO needs testing *)
(** lazy_prob_intervals_from_freq [freq] [bounds_mats_list] expects an initial
    frequency for a single population and a LazyList of bounds matrix pairs,
    and returns a LazyList of probability intervals for each timestep. Like
    [lazy_prob_intervals] and [lazy_singleton_intervals] but more efficient
    if the probability interval consists of a single distribution that puts
    all probability on a single frequency.  Note that the first element will
    be the initial interval: a list containing two copies of a vector with
    1.0 at index freq and 0.0 everywhere else. *)
let lazy_prob_intervals_from_freq freq bounds_mats_list =
  let size, _ = M.shape (fst (LL.hd bounds_mats_list)) in
  let init_dist = (WF.make_init_dist size freq) in
  LL.cons [init_dist; init_dist] (LL.map (freq_mult freq) bounds_mats_list)


(***************************************)
(** Make example intervals *)

(** Make an interval from a popsize and list of fitness structures without
    verifying tightness. *)
let make_wf_interval_no_tight_check popsize fitn_list =
  let tranmats = L.map (WF.make_tranmat popsize) fitn_list in
  let low = L.reduce M.min2 tranmats in
  let high = L.reduce M.max2 tranmats in
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
  then Printf.eprintf "make_wf_interval: Wright-Fisher-based interval is not tight\n";
  tight_low, tight_high

(* example :
let p', q' = WF.(make_wf_interval 100 {w11=1.0; w12=0.3; w22=0.1} {w11=1.0; w12=0.9; w22=0.5});;
let p, q = tighten_mat_interval p' q';;
*)


(* NOT IN USE:
 * (** Use separate compare functions for row and column vectors to avoid
 *     having a test for row vs. col inside the compare function. *)
 * 
 * (** Compare function for use by idx_sort for row vector *)
 * let row_vec_idx_cmp mat j j' =
 *   if M.get mat 0 j > M.get mat 0 j' then 1 
 *   else if M.get mat 0 j < M.get mat 0 j' then -1 
 *   else 0
 * 
 * (** Given a row or column vector, return a list of indexes in
 *     order of the numerical order of the values at those indexes.
 *     Automatically determines whether the argument is a row or a column 
 *     vector, raising an exception if neither. *)
 * let idx_sort v =
 *   let rows, cols = M.shape v in
 *   let size, idx_cmp = if rows = 1
 *                       then cols, row_vec_idx_cmp 
 *                       else rows, col_vec_idx_cmp in
 *   let idxs = L.range 0 `To (size - 1) in
 *   L.fast_sort (idx_cmp v) idxs
 *) 
