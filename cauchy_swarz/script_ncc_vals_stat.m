
           
noises_sigmas=[0.05,0.1,0.2];


%load('boxes1_aspect_2_to_Inf.mat');
load('boxes_sapl1_upto1.mat');


rand_points=[42 21 49 37 29 13 52 10 51 30 40 17 3 44 26 12 5 27 41 32];

rand_points=[43 22 50 38 30 14 53 11 52 31 41 18 4 45 27 13 6 28 1 33];

rand_pats = [9 3 4 7 15 1 11 10 14 2];

ncc_value_stat('ncc_vals_stat.mat',points, noises_sigmas.^2,rand_pats,rand_points)


