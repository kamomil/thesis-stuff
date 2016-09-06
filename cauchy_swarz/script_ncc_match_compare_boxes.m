
debug=1;
threshold=1-0.01;
up_to_wind_per=20;


small=load('boxes2_aspect_1_to_Inf.mat');
large=load('boxes1_aspect_2_to_Inf.mat');   
legend_small='aspect -inf to inf';
legend_large='aspect 2 to inf';

ncc_match_comper_boxes(points,large.filters_boxes,legend_large,small.filters_boxes,legend_small,debug,threshold);