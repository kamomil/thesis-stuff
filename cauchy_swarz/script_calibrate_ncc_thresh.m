

thresholds=0.4:0.02:0.99;%[0.99,0.97,0.95,0.93,0.9,0.87,0.85,0.83,0.8,0.77,0.75,0.73,0.7,0.67];
noises_sigmas=[0.05,0.1,0.2,0.4];
iter=6;


%load('boxes1_aspect_2_to_Inf.mat');   
load('boxes_sapl1_upto1.mat');   

rand_points =[30,13,37,2,22,6,51,33,18,8,44,19,16,36,11,41,40,38,46,32,31,43,24,25,15,17,35,26,1,28];



  


%rand_points=[42 21 49 37 29 13 52 10 51 30 40 17 3 44 26 12 5 27 41 32];
rand_pats =  [9 13  6 20 19  8 16 14 4 5 7 12 3 11 2];
      
filename='thresh_stat_evning.mat';
calibrate_ncc_thresh2(filename,points,filters_boxes,thresholds, noises_sigmas.^2,iter,rand_pats,rand_points)

