
debug=1;
threshold=0.85;
up_to_wind_per=20;

load('boxes_sapl1_upto1.mat');   
% small=load('boxes1_aspect_3_to_Inf.mat');   
% large=load('boxes1_aspect_2_to_Inf.mat');   
% legend_small='aspect 3 to inf';
% legend_large='aspect 2 to inf';

test_ncc_match_with_cauchy_schwarz(points,filters_boxes,debug,threshold);