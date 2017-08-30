(** Example 11.5 in Hermans & Skulj's "Stochastic Processes" in 
  * Introduction to Imprecise Probabilities. *)

module M = Owl.Mat
module P = Matutils.Probmat

let q1 = P.vec_from_list [0.5 ; 0.5 ; 0.0];;
let q2 = P.vec_from_list [0.0 ; 0.5 ; 0.5];;

let i = M.eye 3;;

let t  = P.mat_from_lists [[0. ; 1. ; 0.];
                           [0. ; 0. ; 1.];
                           [1. ; 0. ; 0.]]

let mix m1 m2 p = 
  let q = 1. -. p in
  M.( (p $* m1) + (q $* m2) );;

let cmix = mix q1 q2

let tmix = mix i t

let qit = M.((0.5 $* q1 *@ i) + (0.5 $* q2 *@ t));;

let ctmix a b = 
  M.((cmix a) *@ (tmix b))

(* Claim: There is no a and no b s.t. ctmix a b = qit *)
