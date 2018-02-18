
module U = Utils.Genl


let run_test f n =
  print_string "start gsl:\n";
  U.time0 (fun () -> for i = 1 to n do f 500 0.5 1000 done) ()
  print_string "done.\n"

let run_gsl = run_test Gsl.Randist.binomial_pdf
let run_owl = run_test Owl.Stats.binomial_pdf

let description = "time something"

let commandline = Command.basic ~summary:description
                                ~readme:(fun () -> description)
                                Spec.(empty +> anon ("N" %: int))
                                (fun () -> run_gsl; run_owl)

let wrap_run version build_info commandline = 
  Command.run ~version:version ~build_info:build_info commandline

wrap_run "0.0" "speed test of something commandline"

