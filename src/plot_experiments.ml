
let f m sd x = Stats.Pdf.gaussian (x +. m) sd in
let h = Plot.create "yo.pdf" in
for i=1 to 50 do 
  let i' = float i in
  let m, sd = i' /. 1000., i' /. 1000. in
  Plot.plot_fun ~h (f m sd) (-0.2) 0.1;
done;
Plot.output h;;

let gauss mean sd x = Stats.Pdf.gaussian (x +. mean) sd;;
let xs = Mat.(((sequential 1 1000) -$ 499.) /$ 100.);;
let ys = Mat.(map (gauss 0. 1.) xs);;
let h = Plot.create "yo.pdf" in Plot.plot ~h xs ys ; Plot.output h;;


open Owl;;
open Owl.Plot;;
open Batteries;;
let gauss mean sd x = Stats.Pdf.gaussian (x +. mean) sd;;
let num_plots = 30;;
let third_num_plots = num_plots / 3;;
let idxs = List.range 1 `To num_plots;;
let xs, ys = Mat.meshgrid (-9.) 4. (-9.) 4. 1000 num_plots;;
(* let means = List.map (fun x -> 
                       let x' = float x in
		       if x <= third_num_plots then x' /. 10.
		       else x' /. 8.)
		     idxs;; *)
let means = List.map (fun x -> let x' = float x in x' /. 10.) idxs;;
let sds = List.map (fun x -> (float x) /. 5.) idxs;;
let distlist = List.map2 (fun mean sd -> Mat.(map (gauss mean sd) (row xs 1))) means sds;;
let zs = List.reduce Mat.concat_vertical (List.rev distlist);;
List.map Mat.shape [xs; ys; zs];; (* all three shapes s/b same *)
let h = Plot.create ~m:2 ~n:1 "foo.pdf" in
  Plot.set_background_color h 255 255 255;
  Plot.subplot h 0 0;
  Plot.set_foreground_color h 0 0 0;
  Plot.plot ~h xs zs;
  Plot.subplot h 1 0;
  Plot.set_foreground_color h 0 0 0;
  Plot.mesh ~h ~spec:[NoMagColor; ZLine X; Azimuth 15.; Altitude 30.] xs ys zs;
  Plot.output h;;

