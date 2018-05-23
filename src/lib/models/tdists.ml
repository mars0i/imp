module Mat = Owl.Mat
module L = Batteries.List
module LL = Batteries.LazyList
module G = Utils.Genl

let always_true _ = true

type tdists = {gen : int ; dists : Mat.mat list}

type tdistslist = tdists LL.t

(** tdists functions *) 

(** accessor, constructor functions: *)
let gen tds = tds.gen
let dists tds = tds.dists
let make gen dists = {gen ; dists}

(** tdistslist functions *) 

let hd = LL.hd
let tl = LL.tl
let last = LL.last
let next = LL.next
let is_empty = LL.is_empty
let at = LL.at
let fold_left = LL.fold_left
let map = LL.map
let seq = LL.seq
let cons = LL.cons
let fold_right = LL.lazy_fold_right
let take = LL.take
let take_while = LL.take_while
let drop = LL.drop
let to_list = LL.to_list
let rev = LL.rev
let from_loop = LL.from_loop
let nil = LL.nil



let lazy_range ?(step=1) start stop = 
  let open LL in
  if start = stop then make 1 start
  else let adjust_by, ineq =
    if stop > start then (+), (>) else (-), (<)
  in let rec aux curr stop' =
     if ineq curr stop' then nil
     else lazy (Cons (curr, aux (adjust_by curr step) stop'))
  in aux start stop

(** [lazy_ints ~every:n init_n] returns an infinite sequence of
    integers [~every] apart starting from [init_n].  [every]
    defaults to 1.  Giving it a negative value will produce a
    descending sequence. *)
let lazy_ints ?(skip=1) init_n =
  seq init_n (fun n -> n + skip) always_true

(** In [lazy_select accessor keys data], [keys] and [data] are lazy
    lists.  The function returns a lazy list of elements from [data]
    such that [accessor val] is equal to some element [keys].  Both 
    keys and the values by which elements of [data] are selected should
    be monotonically increasing numbers of the same kind, usually [int]s.
    For example, if [keys] and [data] are increasing sequences of integers,
    and [accessor] is the identity function, [lazy_select] returns the
    intersection of [keys] and [data]. *)
let lazy_select accessor keys data =
  let open LL in
  let rec sel ks ds =
    if is_empty ds || is_empty ks then Nil
    else let k, d = hd ks, hd ds in
         let d_key = accessor d in
         if k = d_key then Cons(d, (lzsel (tl ks) (tl ds)))
         else if k > d_key
	 then sel ks (tl ds)  (* let ds catch up *)
         else sel (tl ks) ds  (* let ks catch up *)
  and lzsel ks ds = lazy (sel ks ds)
  in lzsel keys data

(********************************************)
(** Iteration functions *)

(* Based on Batteries.LazyList.lazy_fold_right *)
(** [lazy_fold_right2 f l1 l2 init_val] folds function [f] over two lazy 
    lists [l1] and [l2], with initial value [init_val].  Note that for
    constructing lazy lists, one must use [Cons] and [nil] rather than
    [cons] and [nil], as with the [LazyList] eager fold functions:
    {[
      let natnos = seq 0 ((+) 1) (fun _ -> true);;
      let posints = seq 1 ((+) 1) (fun _ -> true);;
      let prods = G.lazy_fold_right2 (fun x y acc -> Cons(x*y, acc)) natnos posints TL.nil;;
      (to_list (take 10 prods));;
      - : int list = [0; 2; 6; 12; 20; 30; 42; 56; 72; 90]
    ]} *)
let lazy_fold_right2 f l1 l2 init_val =
  let rec aux rest1 rest2 =
    lazy begin
      match next rest1, next rest2 with
      | Cons (x1, t1), Cons (x2, t2) -> f x1 x2 (aux t1 t2)
      | Nil, Nil | Nil, _ | _, Nil -> Lazy.force init_val
    end
  in
aux l1 l2








(** Convenience abbreviation for [LazyList.(to_list (take n lazy_list))] *)
let take2list n lazy_list = (to_list (take n lazy_list))

(** Return a lazy list that's a sublist of the argument, from element start 
    (zero-based) to element finish, inclusive. *)
let sub_lazy_list start finish ll =
  take (finish - start + 1) (drop start ll)

(** Convenience function: Takes elements from start to finish, inclusive, 
    from a LazyList, and convert the result to a List. *)
let take_to_list start finish ll = 
  to_list (sub_lazy_list start finish ll)

let lazy_take_at_idxs ns ll =
  let f acc elt =
    let i, ns', ll' = acc in
    match ns' with
    | [] -> acc
    | n::nstl -> if i = n
                 then (i+1, nstl, cons elt ll')
                 else (i+1, ns', ll')
  in
  let _, _, result = fold_left f (0, ns, nil) ll in
  rev result





(** Create a power set of integers from a smaller power set.
    Given a sequence of sequences representing the power set of non-negative
    integers up through n-1, returns a pair consisting of the sequence of 
    sequences for up through n-1 and one up through n.  (This kind of pair 
    is what LazyList.from_loop wants.)  Assumes that the first element of 
    the first element of the power set sequence pset passed in is n-1. *)
let next_intsets pset =
  let n = 1 + L.hd (L.hd pset) in  (* Get next integer; previous one must be first in the first element. *)
  let addl_sets = (L.map (fun xs -> n :: xs) pset) in
  (pset, addl_sets @ pset)

(** Generate a list of subsequent integer power sets. 
    Return a lazy list of subsequent power sets of integers from 0 to n. 
    They can be retreived using e.g., to get the power set of integers
    up to 5: LazyList.at intsets 4 .  Note each power set is in the form
    of a regular list; only the top level list is lazy. *)
let make_intsets () = from_loop [[0]; []] next_intsets

(** A lazy list of integer power sets. *)
let algebra_sets = make_intsets ()

(*********** probabilities over algebras **********)


(** Set difference for ordered lists of integers.
    Given lists xs and ys that are both ordered in the same way (e.g. 
    monotonically ordered integers), return a list containing all
    elements of xs not in ys (with order preserved).  The elements
    in ys must be a (perhaps improper) subset of xs.
    This version by RichN at https://codereview.stackexchange.com/a/162407/61384 *)
let rec subtract_list xs ys =
  match xs, ys with 
  | [], _ -> []
  | _, [] -> xs
  | x::xs', y::ys' when x = y -> subtract_list xs' ys'
  | x::xs', _ -> x::(subtract_list xs' ys)

(** Set complement for ordered lists of integers.
    i.e. return the atoms representing the negation of the original set.
    Given a maximum element omega_max and a subset of a domain of atoms 
    represented by list of integers in decreasing order, return a list
    of the integers, between 0 and omega_max inclusive, that are not in
    the subset.  The returned list will also be in decreasing order.
    Note that omega_max is one less than the size of the domain. *)
let list_complement omega_max subset =
  let omega = L.range omega_max `Downto 0 in
  subtract_list omega subset

(** Given a vector of atom probabilities and a list of indexes representing
    atoms, returns sum of probabilities for the set containing those atoms. *)
let prob_sum probs atom_idxs =
  let add_prob sum idx = sum +. Mat.get probs 0 idx  (* Owl.Mat.get rather than .{i,j} to get type right *)
  in L.fold_left add_prob 0. atom_idxs

(** Return (1 - the sum of probs of extreme probs) for a set of atoms
    See (3) and (4) in Skulj. *)
let invert_prob_sum omega_max atom_extrema subset_idxs = 
  1. -. prob_sum atom_extrema (list_complement omega_max subset_idxs)

(*********** algebras of indexes representing atoms **********)

(** Given a probability vector, returns an alist of pairs for all
    elements in the algebra of sets based on the atoms represented by
    elements in the vector.  Each pair contains a list of indexes 
    representing atoms in a set followed by the probability of that set. *)
let algebra_probs probs = 
  let i = (snd (Mat.shape probs)) - 1 in
  let idx_sets = at algebra_sets i in
  let make_entry event = (event, prob_sum probs event) in
  L.map make_entry idx_sets 

(** Map prob_sum over each possible combination of atoms. *)
let simple_sums omega_max atom_extrema =
  L.map (prob_sum atom_extrema) (at algebra_sets omega_max)

(** Map invert_prob_sum over each possible combination of atoms. *)
let inverted_sums omega_max atom_extrema =
  L.map (invert_prob_sum omega_max atom_extrema) (at algebra_sets omega_max)

(** Calculate L values for all members of the algebra and return an
    (atoms, L-value) alist.  See (3) in Skulj. *)
let pri_f_field_lowers omega_max atom_mins atom_maxs =
  let mins = simple_sums omega_max atom_mins in
  let inverted_maxs = inverted_sums omega_max atom_maxs in
  let minmins = L.map2 max mins inverted_maxs in
  L.combine (at algebra_sets omega_max) minmins

(** Calculate U values for all members of the algebra and return an
    (atoms, U-value) alist.  See (4) in Skulj. *)
let pri_f_field_uppers omega_max atom_mins atom_maxs =
  let maxs = simple_sums omega_max atom_maxs in
  let inverted_mins = inverted_sums omega_max atom_mins in
  let maxmaxs = L.map2 min maxs inverted_mins in
  L.combine (at algebra_sets omega_max) maxmaxs

(*

(*********** Strings for printing **********)

(* Jane Street Core has a function that does this btw. *)
(** Given a function that will create strings from the elements of list l,
    return a string representing l. *)
let string_of_t_list string_of_t l = 
  "[" ^ String.concat "; " (L.map string_of_t l) ^ "]"

(** Given two functions that will create strings from the elements of the pair,
    return a string representing it. *)
let string_of_pair string_of_t1 string_of_t2 (t1, t2) =
  "(" ^ string_of_t1 t1 ^ ", " ^ string_of_t2 t2 ^ ")"

(** Return a string representing a list of integers. *)
let string_of_int_list = string_of_t_list string_of_int

(** Return a string representing a list of floats. *)
let string_of_float_list = string_of_t_list string_of_float

(** Convert one algebra entry, a record containing a list of indexes
   and a probability, into a string. *)
let string_of_alg_prob = string_of_pair string_of_int_list string_of_float

(** Convert a list of indexes, probability entries into a string. *)
let string_of_alg_probs = string_of_t_list string_of_alg_prob

(** Return a string representation of a pair of floats. *)
let string_of_float_pair = string_of_pair string_of_float string_of_float

(** Convert one algebra interval entry, a record containing a list of indexes
   and a list of a lower and an upper probability, into a string. *)
let string_of_alg_interval = 
  string_of_pair string_of_int_list string_of_float_pair

(** Convert a list of indexes, probability interval entries into a string. *)
let string_of_alg_intervals = string_of_t_list string_of_alg_interval

(*********** Convenience function for testing **********)

(** Convenience function for testing.  Should be passed the size of
    a set of atoms omega and the number of probability distributions
    over them that is wanted.  Returns
       ((ps, mins, maxs),
        (algs, min_alg, max_alg),
        (f_lowers, f_mins, f_inverted_maxs),
        (f_uppers, f_maxs, f_inverted_mins),
        f_intervals) *)
let generate_system omega_size num_dists =
  let omega_max = omega_size - 1 in
  (* make num_dists atomic dists *)
  let ps = L.init num_dists (fun _ -> unif_stoch_vec omega_size) in
  (* probabilities for algebras for each of the num_dists distributions *)
  let algs = L.map algebra_probs ps in (* alists mapping atom lists to probs *)
  (* min and max values of atomic probs across all num_dists distributions *)
  let mins = min_elts ps in
  let maxs = max_elts ps in
  (* min and max values of probs for each member of the algebra *)
  let min_alg = min_algebra_probs algs in
  let max_alg = max_algebra_probs algs in
  (* prob values for each member of the algebra computed using (3) in Skulj 
   * The first two are min'ed to produce the third. *)
  let f_mins = simple_sums omega_max mins in
  let f_inverted_maxs = inverted_sums omega_max maxs in
  let f_lowers = pri_f_field_lowers omega_max mins maxs in
  (* prob values for each member of the algebra computed using (4) in Skulj
   * The first two are max'ed to produce the third. *)
  let f_maxs = simple_sums omega_max maxs in
  let f_inverted_mins = inverted_sums omega_max mins in
  let f_uppers = pri_f_field_uppers omega_max mins maxs in
  (* interval enteries constructed from f_lowers, f_uppers *)
  let f_intervals = pri_f_field_intervals f_lowers f_uppers in
  (* return all of the above: *)
  ((ps, mins, maxs),
   (algs, min_alg, max_alg),
   (f_lowers, f_mins, f_inverted_maxs),
   (f_uppers, f_maxs, f_inverted_mins),
   f_intervals)



(* Transform lazy lists of dists to/from lazy lists of tdists: *)
let ints_from n = seq n ((+) 1) always_true

let add_gens ?(first_tick=0) dists_llist =
  map2 make (ints_from first_tick) dists_llist

let remove_gens tdists_llist = map dists tdists_llist

(** [sublist start_t finish_t tdists_llist] returns a lazy list that's
    a finite sublist of [tdists_llist], from the first element with 
    [gen] >= [start_gen] to the last element with [gen] <= [finish_gen].  
    Note that if the list is infinite and there are no elements satisfying
    both of these conditions, the function will try to run forever. *)
let sublist start_gen finish_gen tdists_llist =
  take_while (fun tds -> tds.gen <= finish_gen)
                (drop_while (fun tds -> tds.gen < start_gen)
		               tdists_llist)

(** In [select_by_gens generations tdists_llist], [generations] is a lazy
    list of integers in increasing order, and [tdists_llist] is a lazy
    list of tdists.  The function returns a lazy list contanining those 
    tdists whose generation numbers match the integers in [generations]. *)
let select_by_gens generations tdists_llist =
  G.lazy_select gen generations tdists_llist

*)
