run("8-bit");

 run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 stack");

run("Set Measurements...", "area center stack redirect=None decimal=3");

run("Shape Filter", "area=150-1500 area_convex_hull=0-Infinity perimeter=0-Infinity perimeter_convex_hull=0-Infinity feret_diameter=0-Infinity min._feret_diameter=0-Infinity long_side_min._bounding_rect.=0-Infinity short_side_min._bounding_rect.=0-Infinity aspect_ratio=1-Infinity area_to_perimeter_ratio=0-Infinity circularity=0-Infinity elongation=0-1 convexity=0-1 solidity=0-1 num._of_holes=0-1 thinnes_ratio=0-1 contour_temperatur=0-1 fractal_box_dimension=0-2 option->box-sizes=2,3,4,6,8,12,16,32,64 add_to_manager draw_holes black_background fill_results_table stack");
