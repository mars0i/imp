module Mat = Owl.Mat
module L = Batteries.List
module LL = Batteries.LazyList

let always_true _ = true

type t = {gen : int ; dists : Mat.mat list}
(** accessor, constructor functions: *)
let gen tds = tds.gen
let dists tds = tds.dists
let make_tdists gen dists = {gen ; dists}

(* Transform lazy lists of dists to/from lazy lists of tdists: *)
let ints_from n = LL.seq n ((+) 1) always_true

let add_gens ?(first_tick=0) dists_llist =
  LL.map2 make_tdists (ints_from first_tick) dists_llist

let remove_gens tdists_llist = LL.map dists tdists_llist

(** [sublist start_t finish_t tdists_llist] returns a lazy list that's
    a finite sublist of [tdists_llist], from the first element with 
    [gen] >= [start_gen] to the last element with [gen] <= [finish_gen].  
    Note that if the list is infinite and there are no elements satisfying
    both of these conditions, the function will try to run forever. *)
let sublist start_gen finish_gen tdists_llist =
  LL.take_while (fun tds -> tds.gen <= finish_gen)
                (LL.drop_while (fun tds -> tds.gen < start_gen)
		               tdists_llist)

(** In [select_by_gens generations tdists_llist], [generations] is a lazy
    list of integers in increasing order, and [tdists_llist] is a lazy
    list of tdists.  The function returns a lazy list contanining those 
    tdists whose generation numbers match the integers in [generations]. *)
let eager_select_by_gens generations tdists_llist =
  let open LL in
  let rec select gs tds =
    if is_empty tds || is_empty gs then nil
    else 
    let g, td = hd gs, hd tds in
    let tdg = td.gen in
    Printf.printf "%d %d\n" g tdg;
    if g = tdg then cons td (select (tl gs) (tl tds))
    else if g > tdg then select gs (tl tds) (* let tds catch up *)
    else select (tl gs) tds                 (* let gs catch up *)
  in select generations tdists_llist

(* THIS METHOD SEEMS TO WORK! *)
let lazy_select ?(accessor=(fun x -> x)) keys vals =
  let rec sel ks vs =
    if LL.is_empty vs || LL.is_empty ks then LL.Nil
    else let k, v = LL.hd ks, LL.hd vs in
         let v_key = accessor v in
         if k = v_key then LL.Cons(v, (lzsel (LL.tl ks) (LL.tl vs)))
         else if k > v_key then sel ks (LL.tl vs) (* let vs catch up *)
         else sel (LL.tl ks) vs                   (* let ks catch up *)
  and lzsel ks vs = lazy (sel ks vs)
  in lzsel keys vals

(* This doesn't work:
let select_by_gens generations tdists_llist =
  lazy_select ~accessor:gen generations tdists_llist
*)
