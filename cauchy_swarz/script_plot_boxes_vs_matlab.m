stdrange=1;
%max_iter=[10,20,40,250,300,400];
%max_iter=[30,30,100,300,400,500];
max_iter=[40,50,200,400,600,800];
%max_iter=[100,100,400,600,600,800];
str_iters=num2str(max_iter);
str_iters(str_iters==' ')='_';
%file=['runtime_mean_minus_',num2str(stdrange),'stds_iters',str_iters,'_frac_mu_minus_1sigma.mat']
%file=['runtime_mean_minus_',num2str(stdrange),'stds_iters',str_iters,'.mat']
%file='runtime_mean_minus_0.5stds_iters100__100__400__600__600__800_frac_mu_minus_0.2sigma.mat';
file='runtime_mean_minus_0.5stds_iters100__100__400__600__600__800_my_frac.mat';
%file='runtime_mean_minus_2stds.mat';
load(file)
%load('runtime2.mat')
%load('runtime_mean_minus_2stds.mat')
%filters_boxes  max_iter       noises_vars    points         res            thresholds
%
% res =
%
% 53x20 struct array with fields:
%
%     pass
%     boxes_time
%     matlab_time
res
sizes=[16,32,48,64,96,128];
sizes_nm=zeros(length(noises_vars),length(sizes));
sigmas=sqrt(noises_vars);
passed=zeros(length(sigmas),length(sizes));
boxes_time=zeros(length(sigmas),length(sizes));
matlab_time=zeros(length(sigmas),length(sizes));
frac=zeros(length(sigmas),length(sizes));
frac2=zeros(length(sigmas),length(sizes));

maxed_passed=zeros(length(sigmas),length(sizes));
maxed_boxes_time=zeros(length(sigmas),length(sizes));
maxed_matlab_time=zeros(length(sigmas),length(sizes));

maxed_boxes_time_img_idx=zeros(length(sigmas),length(sizes));
maxed_boxes_time_pat_idx=zeros(length(sigmas),length(sizes));

% maxed_matlab_time_img_idx=zeros(length(sigmas),length(sizes));
% maxed_matlab_time_pat_idx=zeros(length(sigmas),length(sizes));

maxed_matlab_time_patidx=zeros(length(sigmas),length(sizes));
maxed_matlab_time_imgidx=zeros(length(sigmas),length(sizes));
        


for i=1:length(points)
    for k=1:length(points(i).pats)
        [i,k]
        %res(i,k)
%          if res(i,k).matlab_pass(3)==0
%              continue;
%          end
        
         
         
         matpass=(res(i,k).matlab_pass==1);
        szidx=find(sizes==points(i).pats(k).sz(1));
        sizes_nm(matpass,szidx)=sizes_nm(matpass,szidx)+1;
       %if sum(sum(res(i,k).pass(:,:)==0))>0
%             [i,k]
%             I=mat2gray(imread(points(i).im_name));
%             tl=points(i).pats(k).top_left;
%             sz=points(i).pats(k).sz(1);
%             figure(10001)
%             imshow(I(tl(1):tl(1)+sz-1,tl(2):tl(2)+sz-1))
            
            %         if k>1
            %             points(i).pats(k)=points(i).pats(k-1);
            %             filters_boxes(i,k).boxes=filters_boxes(i,k-1).boxes;
            %             res(i,k).pass=res(i,k-1).pass;
            %             res(i,k).boxes_time=res(i,k-1).boxes_time;
            %             res(i,k).fracs=res(i,k-1).fracs;
            %             res(i,k).matlab_time=res(i,k-1).matlab_time;
            %         else
            %             points(i).pats(k)=points(i).pats(k+1);
            %             filters_boxes(i,k).boxes=filters_boxes(i,k+1).boxes;
            %             filters_boxes(i,k).boxes=filters_boxes(i,k+1).boxes;
            %             res(i,k).pass=res(i,k+1).pass;
            %             res(i,k).boxes_time=res(i,k+1).boxes_time;
            %             res(i,k).fracs=res(i,k+1).fracs;
            %             res(i,k).matlab_time=res(i,k+1).matlab_time;
            %         end
            %
            %         figure(10002)
            %         imshow(I)
            %         figure(10003)
            %         for m=1:20
            %             tl=points(i).pats(m).top_left;
            %             sz=points(i).pats(m).sz(1);
            %             subplot(5,4,m)
            %             imshow([ [I(tl(1):tl(1)+sz-1,tl(2):tl(2)+sz-1), zeros(sz,128-sz)] ; zeros(128-sz,128)])
            %
            %             title(num2str(sz))
            %         end
        %    input('')
        %end
        passed(matpass,szidx)=passed(matpass,szidx)+res(i,k).pass(matpass);
        boxes_time(matpass,szidx)=boxes_time(matpass,szidx)+res(i,k).boxes_time(matpass);
        matlab_time(matpass,szidx)=matlab_time(matpass,szidx)+res(i,k).matlab_time(matpass);
        
        frac(matpass,szidx)=frac(matpass,szidx)+res(i,k).fracs(matpass);
        frac2(matpass,szidx)=frac2(matpass,szidx)+res(i,k).fracs(matpass).^2;
        
        
        
        
        
        max([maxed_boxes_time(matpass,szidx) , res(i,k).boxes_time(matpass)],[],2);
        max([maxed_boxes_time(matpass,szidx) , res(i,k).boxes_time(matpass)],[],2)
        %[maxed_boxes_time(matpass,szidx),idx]=max([maxed_boxes_time(matpass,szidx) , res(i,k).boxes_time(matpass)],[],2);
        
        [m,idx]=max([maxed_boxes_time(:,szidx) , res(i,k).boxes_time],[],2);
        maxed_boxes_time(matpass,szidx)=m(matpass);
        
        maxed_boxes_time(matpass,szidx)
        
        %maxed_boxes_time_img_idx(idx==2,szidx)=i;
        %maxed_boxes_time_pat_idx(idx==2,szidx)=k;text(5.75,sin(2.5),txstr,'HorizontalAlignment','right')
        
        %[maxed_matlab_time(:,szidx),idx]=max([maxed_matlab_time(:,szidx) , res(i,k).matlab_time],2);
        size(idx)
        size(matpass)
        maxed_matlab_time_patidx(matpass & idx==2,szidx)=k;
        maxed_matlab_time_imgidx(matpass & idx==2,szidx)=i;
        
        %maxed_matlab_time_img_idx(idx==2,szidx)=i;
        %maxed_matlab_time_pat_idx(idx==2,szidx)=k;
        
        %maxed_matlab_time(:,szidx)=maxed_matlab_time(:,szidx)+res(i,k).matlab_time;
        
        
    end
end
sizes_nm
passed

maxed_matlab_time
maxed_boxes_time
maxed_matlab_time_patidx
maxed_matlab_time_imgidx
        
% for i=1:length(sizes)
%     figure(3*round(stdrange)+max_iter(5)*600)
%     imshow(get_pattern(points,maxed_boxes_time_img_idx(1,i),k));
%     title(num2str())
%     input('');
% end




c={'b','g','r','c','m','k'};
lines={'-o','--+'};
figure(3*round(100*stdrange)+max_iter(1)*1000)
legend_text=cell(length(sizes),1);
for i=1:length(sizes)
    legend_text{i}=['size ',num2str(sizes(i)),' iters',num2str(max_iter(i))];
    plot(sigmas,passed(:,i)'./sizes_nm(:,i)',c{i});
    hold on
    
end
legend(legend_text);
axis([sigmas(1) sigmas(end) 0.9 1])
title(['passed fraction, sigma range: ',num2str(stdrange)],'FontSize',20);
xlabel('Standart deviation','FontSize',20);
ylabel('Fraction','FontSize',20);
%uicontrol(3*round(100*stdrange)+max_iter(1)*1000,'style','text','pos',[0.1 0.1 0.9 0.1],'string','whatever')

%text(0.15,0.92,['sigma range ',num2str(stdrange),' max iters ', num2str(max_iter)],'HorizontalAlignment','right')

hold off


figure(7*round(100*stdrange)+max_iter(1))
legend_text=cell(2*length(sizes),1);

for i=1:length(sizes)
    legend_text{i}=['boxes size ',num2str(sizes(i)),' iters',num2str(max_iter(i))];
    plot(sigmas,boxes_time(:,i)'./sizes_nm(:,i)',[c{i},lines{1}]);
    hold on
    
end
for i=1:length(sizes)
    legend_text{6+i}=['matlab size ',num2str(sizes(i))];
    plot(sigmas,matlab_time(:,i)'./sizes_nm(:,i)',[c{i},lines{2}]);
    hold on
end

legend(legend_text);
axis([sigmas(1) sigmas(end) 0 4])
%axis([noises_vars(1) noises_vars(end) 0.5 1])
title(['boxes and FFT run time , sigma range: ',num2str(stdrange)],'FontSize',20);
xlabel('Standart deviation','FontSize',20);
ylabel('Run Time in seconds','FontSize',20);
hold off

%Y=(boxes_time./(ones(length(sigmas),1)*sizes_nm))';
Y=(boxes_time./sizes_nm)';

%Y=[Y , (mean(matlab_time./(ones(length(sigmas),1)*sizes_nm)))'];
Y=[Y , (mean(matlab_time./sizes_nm))'];

figure(56425+3*round(100*stdrange)  +max_iter(1)*1000)
bar(sizes,Y,'grouped');
ylim([0 3])

labels=cell(1,length(sizes));
for i=1:length(sizes)
    labels{i}=['sz: ',num2str(sizes(i)),' max iters:',num2str(max_iter(i))];
end
legend('noise 0.05','noise 0.1','noise 0.2','FFT')
title(['boxes and FFT runtime , sigma range: ',num2str(stdrange)],'FontSize',20)
xlabel('pattern size', 'FontSize',20)
ylabel('Run time (sec)', 'FontSize',20)
set(gca,'XTickLabel',labels);


figure(58825+3*round(100*stdrange)  +max_iter(1)*1000)

%min_frac=zeros(length(sigmas),length(sizes));
frac_mus=zeros(length(sigmas),length(sizes));
frac_stds=zeros(length(sigmas),length(sizes));

for i=1:length(sizes)
    legend_text{i}=['size ',num2str(sizes(i))];
    %plot(sigmas,frac(:,i)'./sizes_nm(:,i)',c{i});
    mu=frac(:,i)'./sizes_nm(:,i)';
    mu2=frac2(:,i)'./sizes_nm(:,i)';
    errorbar(sigmas,mu,sqrt(mu2-(mu.^2)),c{i});
    frac_mus(:,i)=mu';
    frac_stds(:,i)=(sqrt(mu2-(mu.^2)))';
    hold on
end
frac_mus
frac_stds
legend(legend_text);
axis([sigmas(1) sigmas(end) 0 0.3])
%axis([noises_vars(1) noises_vars(end) 0.5 1])
title(['fraction of indexes left in for direct convolution, sigma range:  ',num2str(stdrange)],'FontSize',20);
xlabel('Standart deviation','FontSize',20);
ylabel('fraction','FontSize',20);
hold off

