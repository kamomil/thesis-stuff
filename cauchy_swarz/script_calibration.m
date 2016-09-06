clear all
load('../generate_points_and_learn_filters/boxes_sapl1_upto1.mat')
thresholds=[0,0.0001,0.001,0.01,0.1];
noises=[0,0.01,0.05,0.1];
noises=[0,0.1];
iminds=[5,10,15,20,25,30,35,40,45,50];
iminds=[5,10,15];
cauchy_schwarz_calibration(points,iminds,filters_boxes,0,thresholds,noises );