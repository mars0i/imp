(* based on https://mjambon.github.io/atdgen-doc/tutorial *)

open Foo_t

let myfoos = 
  {somefoos = [ {idx = 1; data = [0.2; 0.3; 0.17]};
                {idx = 2; data = [(-17.); 1.07; 2.**0.5; 27.; 15.55; (-1.)]};
                {idx = 3; data = [(-1.); 1.0; 2.**0.25; 27.8; 15.57; (-1.1)]}; ]
  }

(* let () = print_endline (Foo_b.string_of_foos myfoos) *)

let filename = "foo.dat"

let () = Ag_util.Biniou.to_file Foo_b.write_foos filename myfoos

let myfoos_again = Ag_util.Biniou.from_file Foo_b.read_foos filename
  
let () = 
  Printf.printf "Data read from file is same as original? %B\n"
     (myfoos = myfoos_again)

