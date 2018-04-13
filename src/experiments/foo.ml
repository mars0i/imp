(* based on https://mjambon.github.io/atdgen-doc/tutorial *)

open Foo_t

let foo1 = {idx = 1; data = [0.2; 0.3; 0.17]}

let () = print_endline (Foo_j.string_of_foo foo1)

let filename = "foo.dat"

let () = Ag_util.Biniou.to_file Foo_j.write_foo filename foo1

let foo1_again = Ag_util.Json.from_file Foo_j.read_foo filename
  
let () = 
  Printf.printf "%s %B\n" (Foo_j.string_of_foo foo1_again) (foo1 = foo1_again)

