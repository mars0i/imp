(*
 * For large sets of 2D plots, there are too few distinguishable colors
 * to differentiate curves by color, but giving each curve its own plot
 * produces too many plots for convenient viewing.  It's sometimes useful
 * to see look at such curves both overlaid on one 2D plot, and also laid 
 * out as if in three dimensions, with each 2D plot seeming to be in its 
 * its own plane.  Whether this is useful depends on the curves, of course.
 * It can be especially useful when you want to compare related families of
 * curves.
 * 
 * In the example below, notice that overlapping curves make 2D plots confusing
 * in some regions.  The 3D plots clarify the patterns in each family. There 
 * are nevertheless relationships that can be seen more precisely in the 2D 
 * plots, so it's useful to display both plots for each family of curves.
 *
 * mars0i (Marshall Abrams)
 *)

let open Owl in
let open Batteries in (* for List.range, List.reduce *)

(* Generate contrived data.  Your data may provide a better illustration! *)
let gauss mean sd x = Stats.Pdf.gaussian (x +. mean) sd in
let num_points = 1000 in
let num_plots = 30 in
let idxs = List.range 1 `To num_plots in
let xs, ys = Mat.meshgrid (-9.) 4. (-9.) 4. num_points num_plots in
let third_num_plots = num_plots / 3 in
let half_num_plots = num_plots / 2 in
let means = List.map (fun a -> let a' = float a in
             	               if a <= third_num_plots then a' /. 5.
                       	       else (a' /. 5.) -. 2.)
                     idxs
in
let sds = List.map (fun a -> let a' = float a in
                             if a <= half_num_plots then a' /. 5.
                             else (a' /. 10.))
                   idxs
in
let distlist = List.map2 (fun mean sd -> Mat.(map (gauss mean sd) (row xs 1)))
                         means sds
in
let zs = List.reduce Mat.concat_vertical (List.rev distlist) in
let shifts = Mat.(repeat ~axis:1 ((sequential num_plots 1) /$ 5.) num_points) in
let zs' = Mat.(zs + shifts) in

(* Plot two families of curves in both two (left) and three (right) versions. *)
(* Note configuration specifications for the mesh and plot functions:
 * "ZLine X" is what allows the 2D plots to appear as separate curves;
 * otherwise surfaces would seem to connect them. The RGB spec overrides
 * the foreground color, which otherwise would be used for the plot lines. *)
let open Owl.Plot in
  (* mesh and plot configuration specifications: *)
  let mesh_spec = [NoMagColor; RGB (0,0,0); ZLine X; 
                   Azimuth 15.; Altitude 30.] in
  let plot_spec = [RGB (0,0,0)] in
  let font_size = 2.75 in
  (* create page and add plots: *)
  let h = create ~m:2 ~n:2 "2D3Dexample.png" in
  set_background_color h 255 255 255;
  (* first family of curves, using zs: *)
  subplot h 0 0;
  set_font_size h font_size;
  set_foreground_color h 120 120 120; (* gray grid, tick marks *)
  plot ~h ~spec:plot_spec xs zs; (* note zs--not ys--gives the y-axis values *)
  subplot h 0 1;
  set_font_size h font_size;
  set_foreground_color h 120 120 120;
  mesh ~h ~spec:mesh_spec xs ys zs; (* here zs gives z-axis values. *)
  (* second family of curves, using zs': *)
  subplot h 1 0;
  set_font_size h font_size;
  set_foreground_color h 120 120 120;
  plot ~h ~spec:plot_spec xs zs';
  subplot h 1 1;
  set_font_size h font_size;
  set_foreground_color h 120 120 120;
  mesh ~h ~spec:mesh_spec xs ys zs';
  output h
