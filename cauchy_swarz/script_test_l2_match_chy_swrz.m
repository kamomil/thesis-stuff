debug=1;
threshold=0.01;
up_to_wind_per=20;
load('boxes_smpl1_upto5.mat');   
% small=load('boxes1_aspect_3_to_Inf.mat');   
% large=load('boxes1_aspect_2_to_Inf.mat');   
% legend_small='aspect 3 to inf';
% legend_large='aspect 2 to inf';

test_l2_match_with_cauchy_schwarz(points,filters_boxes,debug,threshold);