
module M = Owl.Mat
module G = Utils.Genl

module Command = Core.Command
module Spec = Core.Command.Spec;;

let n = 1000 in
let scale = 500. in
let int_iters = 10000 in
let float_iters = float_of_int int_iters in
let g = M.((uniform n n) /$ scale) in

Printf.printf "Matrix is %dx%d random uniform divided by %f\n" n n scale;

Printf.printf "divide-and-conquere int power:\n";
let _ = G.time2 G.int_mpow g int_iters in

Printf.printf "divide-and-conquere float power:\n";
let _ = G.time2 G.float_mpow g float_iters in ()

(*
Printf.printf "simple recursive int power:\n";
let _ = G.time2 G.dot_pow g int_iters in ()
*)
