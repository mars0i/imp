module W = Wf.Wrightfisher
module U = Matutils.Utils
module S = Wf.Setchains
module Slow1 = Wf.Slow1setchains
module Slow3 = Wf.Slow3setchains

module Command = Core.Command
module Spec = Core.Command.Spec

let run_test n () =
  Printf.printf "N=%d\n\n" n;
  print_string "make initial interval: ";
  let p', q' = W.(U.time3 S.make_wf_interval n {w11=1.0; w12=0.3; w22=0.1}
                                               {w11=1.0; w12=0.9; w22=0.5}) in
  print_string "tighten it: ";
  let p, q = U.time2 S.tighten_mat_interval p' q' in

  print_string "\ncalculate lo2, hi2 using Parmap, using make_bounds_mat4: ";
  let parmap_start_time = Unix.gettimeofday() in (* yes silly to do this twice for same function call *)
  let lo2, hi2 = Matutils.Utils.time3 S.make_kth_bounds_mats p q 2 in
  let parmap_duration = Unix.gettimeofday() -. parmap_start_time in

  print_string "\ncalculate lo2, hi2 without Parmap using make_bounds_mat3: ";
  let noparmap_start_time = Unix.gettimeofday() in
  let lo2', hi2' = Matutils.Utils.time3 Slow3.make_kth_bounds_mats p q 2  in
  let noparmap_duration = Unix.gettimeofday() -. noparmap_start_time in

  Printf.printf "\nResults of both calculations are the same? %B\n" ((lo2, hi2) = (lo2', hi2'));
  Printf.printf "Ratio between make_bounds_mat4 and make_bounds_mat3 wall times ~= %f\n"
                (noparmap_duration /. parmap_duration);

  print_string "\nExtra check: calculate lo2, hi2 without Parmap using (earlier) make_bounds_mat1: ";
  let lo2'', hi2'' = Matutils.Utils.time3 Slow1.make_kth_bounds_mats p q 2  in
  Printf.printf "Results are the same? %B\n" ((lo2, hi2) = (lo2'', hi2''));

  print_string "\nOverall time:\n"



let description = "time setchain functions";;

let commandline = Command.basic ~summary:description
                                ~readme:(fun () -> description)
                                Spec.(empty +> anon ("N" %: int))
                                run_test;;

let wrap_run version build_info commandline = 
  Command.run ~version:version ~build_info:build_info commandline
in
U.time3 wrap_run "1.1" "setchaintest, (c) 2017 Marshall Abrams" commandline;

