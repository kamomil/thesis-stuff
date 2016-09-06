l1=load('frac_stat_sun.mat')


l2=load('frac_stat_tmp3.mat')


for i=1:2%length(l2.measures)
for j=1:10
%l1.measures(i,j).run_time-l2.measures(i,j).run_time
l1.measures(i,j).fracs-l2.measures(i,j).fracs
end
end