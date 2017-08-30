module L  = Batteries.List
module M = Owl.Mat
module P = Matutils.Probmat

let omega_size = 2
let num_dists = 2

let ((ps, mins, maxs),
     (algs, min_alg, max_alg),
     (f_lowers, f_mins, f_inverted_maxs),
     (f_uppers, f_maxs, f_inverted_mins),
     f_intervals) = Probmat.generate_system omega_size num_dists;;

print_string "\natomic dists:\n";;
L.iter (fun m -> M.print m) ps;;
print_string "\nmin dist:";;
M.print mins;;
print_string "\nmax dist:";;
M.print maxs;;
print_string "\nalgebra maxs:\n";;
Printf.printf "%s\n" (P.string_of_alg_probs max_alg);;
print_string "\nalgebra mins:\n";;
Printf.printf "%s\n" (P.string_of_alg_probs min_alg);;
print_string "\nalgebra uppers:\n";;
Printf.printf "%s\n" (P.string_of_alg_probs f_uppers);;
print_string "\nalgebra lowers:\n";;
Printf.printf "%s\n" (P.string_of_alg_probs f_lowers);;
print_string "\nalgebra intervals:\n";;
Printf.printf "%s\n" (P.string_of_alg_intervals f_intervals);;
print_string "\nYow!\n";;


