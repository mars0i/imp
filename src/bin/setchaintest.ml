module W = Wf.Wrightfisher
module U = Matutils.Utils
module S = Wf.Setchains

let n = 300;;

Printf.printf "N=%d\n" n;;

let p', q' = W.(S.make_wf_interval n {w11=1.0; w12=0.3; w22=0.1} {w11=1.0; w12=0.9; w22=0.5});;

print_string "tighten: ";;
let p, q = U.time2 S.tighten_mat_interval p' q';;

let next_bounds_test mk_bounds p_mat q_mat = 
  mk_bounds S.recombine_lo p_mat q_mat p_mat,
  mk_bounds S.recombine_hi p_mat q_mat q_mat;;


print_string "current default lo, hi: ";;
let lo2, hi2 = Matutils.Utils.time3 S.make_kth_bounds_mats p q 2;;

print_string "mat1: ";;
let lo2', hi2' = Matutils.Utils.time3 next_bounds_test S.make_bounds_mat1 p q in
Printf.printf "lo: %B\n" (lo2 = lo2');
Printf.printf "hi: %B\n" (hi2 = hi2');;

print_string "mat2: ";;
let lo2', hi2' = Matutils.Utils.time3 next_bounds_test S.make_bounds_mat2 p q in
Printf.printf "lo: %B\n" (lo2 = lo2');
Printf.printf "hi: %B\n" (hi2 = hi2');;

print_string "mat3: ";;
let lo2', hi2' = Matutils.Utils.time3 next_bounds_test S.make_bounds_mat3 p q in
Printf.printf "lo: %B\n" (lo2 = lo2');
Printf.printf "hi: %B\n" (hi2 = hi2');;

print_string "mat4: ";;
let lo2', hi2' = Matutils.Utils.time3 next_bounds_test S.make_bounds_mat4 p q in
Printf.printf "lo: %B\n" (lo2 = lo2');
Printf.printf "hi: %B\n" (hi2 = hi2');;
