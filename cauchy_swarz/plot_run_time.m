clear all
load('frac_stat_zevel.mat')
load('boxes_sapl1_upto1.mat')

[mr,mc]=size(measures);
sizes = [16    32    48    64    96   128];
sizes_nm=zeros(1,6);
fracs=struct;
[l1,l2]=size(measures(1,1).fracs);
for i=1:6
    fracs(i).f=double(zeros(l1,l2));
    fracs(i).r_per_pixel=zeros(l1,l2);
    fracs(i).dconv_r_per_pixel=zeros(l1,1);
    
    fracs(i).fracs_from_before=zeros(l1,l2);
end


for pointsidx=1:mr
    for patidx=1:mc
        
        szidx=find(sizes==points(rand_points(pointsidx)).pats(rand_pats(patidx)).sz(1));
        I=imread(points(rand_points(pointsidx)).im_name);
        total_ind_num=(size(I,1)-sizes(szidx)+1)*(size(I,2)-sizes(szidx)+1);
        sizes_nm(szidx)=sizes_nm(szidx)+1;
        fracs(szidx).f=fracs(szidx).f+double(measures(pointsidx,patidx).fracs)/double(total_ind_num);
        ind_num=measures(pointsidx,patidx).fracs;
        ind_num=[ones(l1,1)*total_ind_num,ind_num(:,1:end-1)];
        if(szidx==3 && ind_num(2,end)==0)
             pointsidx
             patidx
         end
        ind_num(ind_num==0)=1;
        fracs(szidx).fracs_from_before=fracs(szidx).fracs_from_before+double(measures(pointsidx,patidx).fracs)./double(ind_num);
        
        fracs(szidx).r_per_pixel=fracs(szidx).r_per_pixel+double(measures(pointsidx,patidx).run_time)./double(ind_num);
        
        ind_num_left_to_dconv=measures(pointsidx,patidx).fracs(:,end);
        
        ind_num_left_to_dconv(ind_num_left_to_dconv==0)=1;
        
        %measures(pointsidx,patidx).dconv_run_time
        %ind_num_left_to_dconv
        fracs(szidx).dconv_r_per_pixel=fracs(szidx).dconv_r_per_pixel+double(measures(pointsidx,patidx).dconv_run_time')./double(ind_num_left_to_dconv);
        
    end
end

legend_text=cell(9,1);
to_plot=zeros(0,l2);
next=1;
figure(1)

plot(1:l2,fracs(1).f(1,:)/sizes_nm(1),'b-o')
hold on
plot(1:l2,fracs(1).f(2,:)/sizes_nm(1),'g-o')
hold on
plot(1:l2,fracs(1).f(3,:)/sizes_nm(1),'r-o')
hold on
plot(1:l2,fracs(4).f(1,:)/sizes_nm(1),'b-*')
hold on
plot(1:l2,fracs(4).f(2,:)/sizes_nm(1),'g-*')
hold on
plot(1:l2,fracs(4).f(3,:)/sizes_nm(1),'r-*')
hold on
plot(1:l2,fracs(6).f(1,:)/sizes_nm(1),'b')
hold on
plot(1:l2,fracs(6).f(2,:)/sizes_nm(1),'g')
hold on
plot(1:l2,fracs(6).f(3,:)/sizes_nm(1),'r')
title('fraction of indexes left vs iteration')
legend_text{1}=['size 16, var ',num2str(noises_vars(1))];
legend_text{2}=['size 16, var ',num2str(noises_vars(2))];
legend_text{3}=['size 16, var ',num2str(noises_vars(3))];
legend_text{4}=['size 64, var ',num2str(noises_vars(1))];
legend_text{5}=['size 64, var ',num2str(noises_vars(2))];
legend_text{6}=['size 64, var ',num2str(noises_vars(3))];
legend_text{7}=['size 128, var ',num2str(noises_vars(1))];
legend_text{8}=['size 128, var ',num2str(noises_vars(2))];
legend_text{9}=['size 128, var ',num2str(noises_vars(3))];


legend(legend_text)

hold off

% for i=[1,4,6]
% fracs(i).f=fracs(i).f/sizes_nm(i); 
% to_plot=[to_plot ; fracs(i).f];
% legened_text{next}=['var ',num2str(noises_vars(1)),' sizes ',num2str(sizes(i))];
% next=next+1;
% legened_text{next}=['var ',num2str(noises_vars(2)),' sizes ',num2str(sizes(i))];
% next=next+1;
% legened_text{next}=['var ',num2str(noises_vars(3)),' sizes ',num2str(sizes(i))];
% next=next+1;
% end
% 
% 
% plot(1:l2,fracs(i).f)
% plot(1:l2,to_plot)
% title('fraction of indexes left vs iteration')
% legend(legened_text)
%legend('16','32','48','64','96','128');
%legend(num2str(noises_vars(1)), num2str(noises_vars(2)) , num2str(noises_vars(3)) , num2str(noises_vars(4)) )

%hold off

as=zeros(4,0);

figure(2)
for i=1:6
fracs(i).r_per_pixel=fracs(i).r_per_pixel/sizes_nm(i);    
%figure(i)
hold on
as=[as,fracs(i).r_per_pixel(:,2:end)];
plot(1:l2,fracs(i).r_per_pixel/1000)
title('run time per pixel')
axis([1 l2 0 1])
%legend('16','32','48','64','96','128');
%legend(num2str(noises_vars(1)), num2str(noises_vars(2)) , num2str(noises_vars(3)) , num2str(noises_vars(4)) )
end
hold off
a0s=as(:,1);
as=as(:,2:end);
var(as(:));
a=mean(as(:));
a0=mean(a0s(:));


figure(3)
alphs=zeros(l1,l2);
for i=1:6
fracs(i).fracs_from_before=fracs(i).fracs_from_before/sizes_nm(i);    
alphs=alphs+fracs(i).fracs_from_before;
%figure(i)
hold on
title('alpha');
plot(1:l2,fracs(i).fracs_from_before);

%legend('16','32','48','64','96','128');
%legend(num2str(noises_vars(1)), num2str(noises_vars(2)) , num2str(noises_vars(3)) , num2str(noises_vars(4)) )
end

hold off
alphs=alphs/6;
figure(31)
plot(1:l2,alphs)
title('fraction of indexes from the indexe iteration before')
legend(['var ',num2str(noises_vars(1))],['var ',num2str(noises_vars(2))],['var ',num2str(noises_vars(3))])
alpha0=alphs(:,1);
alpha=mean(alphs(:,2:5),2);

dconvs_run_time_per_pixel=zeros(6,l1);
for i=1:6
fracs(i).dconv_r_per_pixel=fracs(i).dconv_r_per_pixel/sizes_nm(i);
dconvs_run_time_per_pixel(i,:)=fracs(i).dconv_r_per_pixel'/(sizes(i)^2);
end
figure(6)
plot(1:l1,dconvs_run_time_per_pixel);
title('direct donvolution runtime per pixel vs noise var');
legend('16','32','48','64','96','128');
%legend(num2str(noises_vars(1)), num2str(noises_vars(2)) , num2str(noises_vars(3)) , num2str(noises_vars(4)) )
b=mean(dconvs_run_time_per_pixel(:));






