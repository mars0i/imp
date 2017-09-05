module W = Wf.Wrightfisher
module U = Matutils.Utils
module S = Wf.Setchains

let n = 300;;
let p', q' =
  W.(S.make_wf_interval n {w11=1.0; w12=0.3; w22=0.1} {w11=1.0; w12=0.9; w22=0.5});;
let p, q = U.time2 S.tighten_mat_interval p' q';;
let lo2, hi2 = Matutils.Utils.time3 S.make_kth_bounds_mats p q 2;;
