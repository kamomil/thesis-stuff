debug=0;
threshold=10000;
up_to_wind_per=20;
load('../mat_files/boxes_sample1_upto5per_surf_im.mat')     
run_l2_match_with_cauchy_schwarz(points,filters_boxes,debug,threshold);