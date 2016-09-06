
%noises_sigmas=[0,0.05,0.1,0.2,0.3];
%this thresholds are due to
% thresholds=  [0.99,0.99,0.99,0.99,0.99,0.99 ;
%               0.98,0.98,0.98,0.98,0.98,0.98 ;
%               0.92 ,0.82 ,0.8 ,0.9 ,0.9 ,0.86 ;
%               0.68 ,0.7,0.66,0.74,0.72,0.76 ;
%               0.5,0.66,0.62,0.66,0.66,0.72];

thresholds=  [0.99,0.99,0.99,0.99,0.99,0.99 ;
              0.98,0.98,0.98,0.98,0.98,0.98 ;
              0.8 ,0.8 ,0.8 ,0.8 ,0.8 ,0.8 ;
              0.68 ,0.68,0.68,0.68,0.68,0.68 ;
              0.6,0.6,0.6,0.6,0.6,0.7];
          
%sew to thresh_stat_evning.mat
thresholds = [ 0.9400    0.9800    0.9800    0.9800    0.9800    0.9800 ;
               0.7400    0.8600    0.9200    0.9600    0.9000    0.8600 ;
               0.5800    0.7000    0.760    0.720    0.760    0.740 ;
               0.4600    0.5800    0.6000    0.6000    0.660    0.640 ];
           
           
                  
 thresholds = [ 0.9200    0.9200    0.9200    0.920    0.920   0.920 ;
                0.7200    0.7200    0.7200    0.720    0.720   0.720 ;
                0.5600    0.5600    0.5600    0.560    0.560   0.560 ;
                0.4200    0.4200    0.4200    0.420    0.420   0.420 ];
           
noises_sigmas=[0.05,0.1,0.2];

iters=10;
%load('boxes1_aspect_2_to_Inf.mat');
load('boxes_sapl1_upto1.mat');


rand_points=[42 21 49 37 29 13 52 10 51 30 40 17 3 44 26 12 5 27 41 32];

rand_points=[43 22 50 38 30 14 53 11 52 31 41 18 4 45 27 13 6 28 1 33];

rand_pats = [9 3 4 7 15 1 11 10 14 2];

filename='frac_stat_zevel2.mat';
calibrate_ncc_frac_for_direct(filename,points,filters_boxes, noises_sigmas.^2,iters,rand_pats,rand_points,thresholds)

