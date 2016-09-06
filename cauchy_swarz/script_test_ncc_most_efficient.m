



clear all
%load('boxes1_aspect_2_to_Inf.mat');   
load('boxes_sapl1_upto1.mat');   
%load('frac_stat_mon.mat')

%sew to thresh_stat_evning.mat
thresholds = [ 0.9400    0.9800    0.9800    0.9800    0.9800    0.9800 ;
               0.7400    0.8600    0.9200    0.9600    0.9000    0.8600 ;
               0.5800    0.7000    0.760    0.720    0.760    0.740 ;
               0.4600    0.5800    0.6000    0.6000    0.660    0.640 ];
           
thresholds = [ 0.9200    0.9200    0.9200    0.920    0.920   0.920 ;
               0.7200    0.7200    0.7200    0.720    0.720   0.720 ;
               0.5600    0.5600    0.5600    0.560    0.560   0.560 ;
               0.4200    0.4200    0.420    0.420    0.420   0.420 ];

thresholds = [  0.900    0.900    0.9200    0.920    0.920   0.920 ;
                0.600    0.600    0.600    0.60    0.60   0.60 ;
                0.400    0.500    0.500    0.50    0.50   0.50];
noises_vars=[0.05,0.1,0.2].^2;
means =[

    0.8944    0.9333    0.9545    0.9470    0.9341    0.9216;
    0.7363    0.8130    0.8585    0.8327    0.8198    0.7971;
    0.4982    0.6011    0.6624    0.6162    0.6155    0.5863];

stds = [

    0.0775    0.0598    0.0402    0.0266    0.0692    0.0890;
    0.1537    0.1220    0.0887    0.0726    0.1283    0.1506;
    0.1557    0.1635    0.1334    0.1142    0.1581    0.1742];
    
    
% thresholds =[
% 
%     0.8169    0.8735    0.9142    0.9203    0.8650    0.8327;
%     0.5825    0.6909    0.7698    0.7600    0.6915    0.6465;
%     0.3425    0.4376    0.5290    0.5019    0.4574    0.4121];

% thresholds =[
% 
%     0.7394    0.8136    0.8740    0.8937    0.7958    0.7437;
%     0.4288    0.5689    0.6811    0.6874    0.5631    0.4958;
%     0.1867    0.2742    0.3956    0.3877    0.2993    0.2379];

%max_iter=[10,20,40,250,300,400];
%max_iter=[40,50,200,400,600,800];
max_iter=[100,100,400,600,600,800];

frac_mus_for_sigmarange1_30_30 = [

    0.0022    0.0017    0.0001    0.0001    0.0010    0.0006;
    0.0202    0.0158    0.0007    0.0003    0.0041    0.0032;
    0.1569    0.1681    0.0086    0.0048    0.0188    0.0246];


frac_stds_for_sigmarange1_30_30 = [

    0.0045    0.0041    0.0001    0.0001    0.0047    0.0012;
    0.0164    0.0214    0.0025    0.0010    0.0146    0.0065;
    0.1253    0.1732    0.0124    0.0101    0.0356    0.0384];

    
my_frac = [frac_stds_for_sigmarange1_30_30(:,1:2)/3, frac_stds_for_sigmarange1_30_30(:,3:end)/5];

frac= frac_mus_for_sigmarange1_30_30-my_frac;


%max_iter=[30,30,100,300,400,500];

%load('runtime3.mat')
str_iters=num2str(max_iter);
str_iters(str_iters==' ')='_';
stdrange=0.5;
file=['runtime_mean_minus_',num2str(stdrange),'stds_iters',str_iters,'_my_frac.mat'];
thresholds=means-stdrange*stds;
res=debug_for_opencv('bla',points,filters_boxes,thresholds,noises_vars(1:3),max_iter,frac);
%res=test_ncc_most_efficient('bla',points,filters_boxes,thresholds,noises_vars(1:3),max_iter,frac);

% load(file);
% for i=1:length(points)
%     for k=1:length(points(i).pats)
%         %[i,k]
%         %res(i,k)
%          if isempty(res(i,k).matlab_pass) 
%              res(i,k).matlab_pass = [ 0 0 0 ]';
%          end
%     end
% end
% save(file,'res','noises_vars','points','filters_boxes','thresholds','max_iter');         

% m=load('matlab_success.mat');
% load(file)
% for i=1:53
%     for k=1:20
%         res(i,k).matlab_pass=m.res(i,k).fft_pass;
%     end
% end
% save(file,'res','noises_vars','points','filters_boxes','thresholds','max_iter');