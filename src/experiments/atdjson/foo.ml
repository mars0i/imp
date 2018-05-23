(* based on https://mjambon.github.io/atdgen-doc/tutorial *)

open Foo_t

let myfoos = 
  {somefoos = [ {idx = 1; data = [0.2; 0.3; 0.17]};
                {idx = 2; data = [(-17.); 1.07; 2.**0.5; 27.; 15.55; (-1.)]};
                {idx = 3; data = [(-1.); 1.0; 2.**0.25; 27.8; 15.57; (-1.1)]}; ]
  }

let () = print_endline (Foo_j.string_of_foos myfoos)

let filename = "foo.dat"

let () = Ag_util.Json.to_file Foo_j.write_foos filename myfoos

let myfoos_again = Ag_util.Json.from_file Foo_j.read_foos filename
  
let () = 
  Printf.printf "%s, %B\n"
     (Foo_j.string_of_foos myfoos_again)
     (myfoos = myfoos_again)

