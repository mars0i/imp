
module M = Owl.Mat
module P = Pomap

#require "pomap";;
module P = Pomap_impl.Make(struct type el = int type ord = Unknown | Lower | Equal | Greater let compare x y = if x < 0 || y < 0 then Unknown else if x = y then Equal else if x < y then Lower else Greater end);;
let pm = P.empty;;
pm;;
P.add 2 3 pm;;
pm;;
P.fold (fun n acc -> (P.get_key n, P.get_el n)::acc) pm [] ;;
P.topo_fold (fun n acc -> (P.get_key n, P.get_el n)::acc) pm [] ;;
P.find 2 pm;;
let pm = P.add 2 3 pm;;
P.find 2 pm;;
open P;;
P.find 2 pm;;
P.get_el (P.find 2 pm);;
let (x,y) = P.find 2 pm;;
let (x,y) = P.find 2 pm in P.get_el y;;

module PM = Pomap_impl.Make(struct type el = M.mat type ord = Unknown | Lower | Equal | Greater let compare x y = if x = y then Equal else if x < y then Lower else if x > y then Greater else Unknown end);;

(*

In po_examples.ml and hasse.ml, the keys are the things that get
partially ordered.  These models just unit for *all* of the values.
Note that type el is the *key* type.

When you make an empty pomap, it has type 'a pomap, where 'a is the
value type.  When you add an element to it and return a new pomap
object, it will have the type of whatever kind of value you used.

*)
