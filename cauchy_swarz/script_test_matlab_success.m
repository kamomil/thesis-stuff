




%load('boxes1_aspect_2_to_Inf.mat');   
load('boxes_sapl1_upto1.mat');   
load('frac_stat_mon.mat')

file='matlab_success.mat';

res=test_matlab_success(file,points,noises_vars(1:3));