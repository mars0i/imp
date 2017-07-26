module Mat = Owl.Mat
module Lin = Owl.Linalg.Generic

let non_diag_elts m =
  let open Mat in
  (triu ~k:1 m) + (tril ~k:(-1) m)

let is_diag m =
  let a, b = Mat.shape m in
  if a <> b then raise (Failure "is_diag: Matrix is not square.");
  0. = Mat.sum (non_diag_elts m);;

(** Moore-Penrose inverse for diagonal matrices:
    https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Singular_value_decomposition_.28SVD.29 *)
let diag_pinv m =
  let a, b = Mat.shape m in
  if a <> b then raise (Failure "diag_pinv: Matrix is not square.");
  if not (is_diag m) then raise (Failure "diag_pinv: Matrix is not diagonal.");
  let m' = Mat.clone m in
  for i = 0 to a do
    let elt = Mat.get m' i i in
    if elt <> 0. then Mat.set m' i i (1. /. elt);
  done;
  (Mat.ctranspose m')

(** Moore-Penrose pseudoinverse:
    https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Singular_value_decomposition_.28SVD.29 *)
let pinv m =
  let u, s, vt = Lin.svd m in
  let open Mat in
  (ctranspose vt) *@ (diag_pinv s) *@ (ctranspose u)

