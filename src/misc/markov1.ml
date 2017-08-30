(* experiments with pop freq (normal) Markov procs *)

(* Note that owl comes with a built-in plot module that can make 3D
 * wireframe plots and countour heatmaps as well as various 2D plots. *)

module LL = Batteries.LazyList;;
module P = Matutils.Probmat;;

let s = P.vec_from_int_list [0; 0; 0; 1; 0;]

let t = P.mat_from_lists [[1.0; 0.0; 0.0; 0.0; 0.0];
                          [0.5; 0.3; 0.2; 0.0; 0.0];
                          [0.0; 0.5; 0.3; 0.2; 0.0];
                          [0.0; 0.0; 0.5; 0.3; 0.2];
                          [0.0; 0.0; 0.0; 0.0; 1.0]]

let next v = M.(v *@ t)

let nextand v = (v, next v)

(* https://stackoverflow.com/a/15542110/1455243 *)
let rec applyn n func arg =
  if n <= 0 then arg
  else applyn (n-1) func (func arg)

let states = LL.from_loop s nextand

(* LL.to_list (LL.take 10 states) *)
