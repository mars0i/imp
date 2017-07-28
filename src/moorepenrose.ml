module Mat = Owl.Mat
module Lin = Owl.Linalg.Generic

(** Moore-Penrose inverse for diagonal matrices:
    https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Singular_value_decomposition_.28SVD.29 *)
let pinv_diag s rows cols =
  let _, len = Mat.shape s in
  let s' = Mat.(1. $/ s) in
  let diag_stub = Mat.diagm s' in
  (* pad diag_stub so that it has dims cols x rows (to transpose it): *)
  Mat.pad [[0; cols - len];[0; rows - len]] diag_stub


(** Moore-Penrose pseudoinverse:
    https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Singular_value_decomposition_.28SVD.29 *)
let pinv mat =
  let u, s, vt = Lin.svd ~thin:false mat in
  let open Mat in
  let rows, cols = shape mat in
  let v, ut = ctranspose vt, ctranspose u in
  let sigma_plus = pinv_diag s rows cols in
  (* let (x1, y1), (x2, y2), (x3, y3) = shape v, shape sigma_plus, shape ut in
  Printf.printf "(%d, %d); (%d, %d); (%d, %d)\n" x1 y1 x2 y2 x3 y3; *)
  v *@ sigma_plus *@ ut

