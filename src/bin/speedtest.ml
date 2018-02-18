module U = Utils.Genl
module L = Batteries.List

let run_test name f m =
  let cpu_time, wall_time = Sys.time(), Unix.gettimeofday() in
  let result = L.init m (fun _ -> f 50000 ~p:0.5 ~n:100000) in
  Printf.printf "cpu: %fs, wall: %fs " (Sys.time() -. cpu_time) (Unix.gettimeofday() -. wall_time);
  Printf.printf "%s %d\n" name (L.length result);; (* length to convince compiler I care *)

run_test "gsl" Gsl.Randist.binomial_pdf 1000000;
run_test "owl" Owl.Stats.binomial_pdf   1000000;
