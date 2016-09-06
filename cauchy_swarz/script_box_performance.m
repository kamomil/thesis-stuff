debug=0;
threshold=0.01;
up_to_wind_per=20;
%load('../mat_files/boxes_sample1_upto5per_surf_im.mat') 
load('boxes_smpl1_upto5.mat')
exact=0;
boxes_performance(points,filters_boxes,debug,threshold);