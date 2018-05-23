module Mat = Owl.Mat
module L = Batteries.List
module TL = Tdistslist
module G = Utils.Genl

let always_true _ = true

type t = {gen : int ; dists : Mat.mat list}
(** accessor, constructor functions: *)
let gen tds = tds.gen
let dists tds = tds.dists
let make gen dists = {gen ; dists}

(* Transform lazy lists of dists to/from lazy lists of tdists: *)
let ints_from n = TL.seq n ((+) 1) always_true

let add_gens ?(first_tick=0) dists_llist =
  TL.map2 make (ints_from first_tick) dists_llist

let remove_gens tdists_llist = TL.map dists tdists_llist

(** [sublist start_t finish_t tdists_llist] returns a lazy list that's
    a finite sublist of [tdists_llist], from the first element with 
    [gen] >= [start_gen] to the last element with [gen] <= [finish_gen].  
    Note that if the list is infinite and there are no elements satisfying
    both of these conditions, the function will try to run forever. *)
let sublist start_gen finish_gen tdists_llist =
  TL.take_while (fun tds -> tds.gen <= finish_gen)
                (TL.drop_while (fun tds -> tds.gen < start_gen)
		               tdists_llist)

(** In [select_by_gens generations tdists_llist], [generations] is a lazy
    list of integers in increasing order, and [tdists_llist] is a lazy
    list of tdists.  The function returns a lazy list contanining those 
    tdists whose generation numbers match the integers in [generations]. *)
let select_by_gens generations tdists_llist =
  G.lazy_select gen generations tdists_llist
