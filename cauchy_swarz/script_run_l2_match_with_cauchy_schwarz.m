debug=0;
threshold=0.01;
up_to_wind_per=20;
load('boxes_smpl1_upto5.mat'); 
exact=0;
run_l2_match_with_cauchy_schwarz(points,filters_boxes,debug,threshold,exact);