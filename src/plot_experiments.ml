
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


(*
 * For large sets of 2D plots, there are too few distinguishable colors
 * to differentiate curves by color, but giving each curve its own plot
 * produces too many plots for convenient viewing.  It's sometimes useful
 * to see your 2D curves both in the normal way and also laid out in three 
 * dimensions, with each 2D plot having its own plane.  Whether this is
 * useful will depend on the curves, of course, but it can be especially 
 * useful when you want to compare different families of curves.
 * 
 * In the example, notice that overlapping curves make 2D plots confusing in
 * some regions.  The 3D plots clarify the patterns of curves.  There are 
 * nevertheless relationships that can be seen more precisely in the 2D plots,
 * so it's useful to display both the 2D and 3D plots.
 *)
open Owl;;
open Batteries;; (* for List.range, List.reduce *)
(* Generate contrived data.  Your data may provide a better illustration. *)
let gauss mean sd x = Stats.Pdf.gaussian (x +. mean) sd;;
let num_points = 1000;;
let num_plots = 30;;
let idxs = List.range 1 `To num_plots;;
let xs, ys = Mat.meshgrid (-9.) 4. (-9.) 4. num_points num_plots;;
let third_num_plots = num_plots / 3;;
let half_num_plots = num_plots / 2;;
let means = List.map (fun a -> 
                       let a' = float a in
             		       if a <= third_num_plots then a' /. 5.
             		       else (a' /. 5.) -. 2.)
                     idxs;;
let sds =   List.map (fun a -> 
                       let a' = float a in
             		       if a <= half_num_plots then a' /. 5.
             		       else (a' /. 10.))
                     idxs;; 
let distlist = List.map2 (fun mean sd -> Mat.(map (gauss mean sd) (row xs 1)))
                         means sds;;
let zs = List.reduce Mat.concat_vertical (List.rev distlist);;
let shifts = Mat.(repeat ~axis:1 ((sequential num_plots 1) /$ 5.)  num_points);;
let zs' = Mat.(zs + shifts);;

(* Plot families of 2D curves in both two (left) and three (right) dimensions. *)
let open Owl.Plot in
  (* configuration specifications for the mesh and plot functions. 
   * "ZLine X" is what allows the 2D plots to appear as separate curves;
   * otherwise surfaces would seem to connect them. The RGB spec overrides
   * the foreground color, which otherwise would be used for the plot lines. *)
  let mesh_spec = [NoMagColor; RGB (0,0,0); ZLine X; Azimuth 15.; Altitude 30.] in
  let plot_spec = [RGB (0,0,0)] in
  (* create page and add plots: *)
  let h = create ~m:2 ~n:2 "foo.png" in
  set_background_color h 255 255 255;
  (* first family, using zs: *)
  subplot h 0 0;
  set_foreground_color h 140 140 140; (* gray grid, tick marks *)
  plot ~h ~spec:plot_spec xs zs; (* note zs--not ys--gives the y-axis values *)
  subplot h 0 1;
  set_foreground_color h 140 140 140;
  mesh ~h ~spec:mesh_spec xs ys zs; (* here zs gives z-axis values. *)
  (* second family, using zs': *)
  subplot h 1 0;
  set_foreground_color h 140 140 140;
  plot ~h ~spec:plot_spec xs zs';
  subplot h 1 1;
  set_foreground_color h 140 140 140;
  mesh ~h ~spec:mesh_spec xs ys zs';
  output h;;

(* List.map Mat.shape [xs; ys; zs; zs'];; *) (* all shapes s/b same *)
