(* based on https://mjambon.github.io/atdgen-doc/tutorial *)

open Yo_t

let yo1 = {gen = 1; dists = [0.2; 0.3; 0.17]}

let filename = "yo.dat"

let () =
  print_endline (Yo_j.string_of_yo yo1)

(* based on
 * https://mjambon.github.io/atdgen-doc/tutorial#inspecting-and-pretty-printing-json
 * *)

let () =
  Ag_util.Biniou.to_file Yo_b.write_yo filename yo1
  

