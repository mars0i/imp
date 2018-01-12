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

(** [tdists_sublist start_t finish_t tdists_llist] returns a lazy list that's
    a sublist of [tdists_llist], from the element with [gen]=[start_gen] to the 
    element with [gen]=[finish_gen], inclusive.  Note that if the list is infinite
    and there is no element with [gen]=[start_gen] or with [gen]=[finish_gen], the 
    function will run forever, or until the system is overloaded. *)
let sublist start_gen finish_gen tdists_llist =
  LL.take_while (fun tds -> tds.gen <= finish_gen)
                (LL.drop_while (fun tds -> tds.gen < start_gen)
		               tdists_llist)

let selection generations tdists_llist =
  let add_if_selected td gens_plus_tds =
	let gens, tds = gens_plus_tds in
	if td.gen = (L.hd generations)
        then (L.tl generations), (LL.cons td tds)
	else generations, tds
  in
  let empty_generations, new_tdists_llist = 
	LL.fold_right add_if_selected (generations, LL.nil) tdists_llist in
  new_tdists_llist

    
let select_gens generations tdists_llist =
  let open LL in
  let rec select gens tds =
    Printf.printf "%d %d\n" (hd gens) (hd tds);
    if is_empty tds then tds
    else let g, td = hd gens, hd tds in
         if g = td.gen
         then cons td (select (tl gens) (tl tds))
         else (select (tl gens) (tl tds))
  in select generations tdists_llist
