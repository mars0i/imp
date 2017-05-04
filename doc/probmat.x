
(* NOT equation (3) in Skulj ?? *)
let flowers max_n min_alg max_alg =
  let flower (idxs, lower) = 
    let complement_upper = (L.assoc (list_complement max_n idxs) max_alg)
    in (idxs, (max lower (1. -. complement_upper)))
  in L.map flower min_alg

(* NOT equation (4) in Skulj ?? *)
let fuppers max_n min_alg max_alg =
  let fupper (idxs, upper) = 
    let complement_lower = (L.assoc (list_complement max_n idxs) min_alg)
    in (idxs, (min upper (1. -. complement_lower)))
  in L.map fupper max_alg


let rec subtract_list xs ys =
  match xs with 
  | [] -> []
  | x::more_xs -> match ys with
    | [] -> xs
    | y::more_ys -> 
      if x = y 
      then subtract_list more_xs more_ys
      else x::(subtract_list more_xs ys)
