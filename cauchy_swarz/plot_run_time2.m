clear all
range1=1:5:50;
range2=1:20:300;

load('frac_stat2.mat')
load('boxes_sapl1_upto1.mat')


[mr,mc]=size(measures);
sizes = [16    32    48    64    96   128];
sizes_nm=struct;
fracs=struct;
l1=length(noises_vars);
l2s=length(range1);
l2l=length(range2);

for i=1:3
    fracs(i).f=double(zeros(l1,l2s));
    fracs(i).run_time=zeros(l1,l2s);
    fracs(i).pass=zeros(l1,l2s);
    sizes_nm(i).nm=zeros(1,l2s);
end
for i=4:6
    fracs(i).f=double(zeros(l1,l2l));
    fracs(i).run_time=zeros(l1,l2l);
    fracs(i).pass=zeros(l1,l2l);
    sizes_nm(i).nm=zeros(1,l2l);
end


for pointsidx=1:mr
    for patidx=1:mc
        
        boxes=filters_boxes(rand_points(pointsidx),rand_pats(patidx)).boxes;
        box_arr=boxes(:,1:4);
        l=length(boxes(:,5));
        
        szidx=find(sizes==points(rand_points(pointsidx)).pats(rand_pats(patidx)).sz(1));
        range=range1;
        if szidx>3
            %[l,szidx]
            range=range2;
        end
        I=imread(points(rand_points(pointsidx)).im_name);
        total_ind_num=(size(I,1)-sizes(szidx)+1)*(size(I,2)-sizes(szidx)+1);
        sizes_nm(szidx).nm=sizes_nm(szidx).nm+(l>=range);
        
        fracs(szidx).f=fracs(szidx).f+ ( (double(measures(pointsidx,patidx).fracs)/double(total_ind_num)) .* (ones(3,1) * (l>=range)));
%         ind_num=measures(pointsidx,patidx).fracs;
%         ind_num=[ones(l1,1)*total_ind_num,ind_num(:,1:end-1)];
% %         if(szidx==3 && ind_num(2,end)==0)
% %              pointsidx
% %              patidx
% %          end
%         ind_num(ind_num==0)=1;
%         %fracs(szidx).fracs_from_before=fracs(szidx).fracs_from_before+double(measures(pointsidx,patidx).fracs)./double(ind_num);
%         
        fracs(szidx).run_time=fracs(szidx).run_time+ (double(measures(pointsidx,patidx).run_time) .* (ones(3,1) * (l>=range) ));
        %size(measures(pointsidx,patidx).pass)
        %size(range)
        fracs(szidx).pass=fracs(szidx).pass + ( measures(pointsidx,patidx).pass & (ones(3,1) * (l>=range) ) );
        
%         if(length(find((measures(pointsidx,patidx).pass==0)))>0)
%             [pointsidx,patidx,sizes(szidx)]
%             measures(pointsidx,patidx).pass
%             tl=points(rand_points(pointsidx)).pats(rand_pats(patidx)).top_left;
%             br=tl+sizes(szidx)-1;
%             
%             pattern=I(tl(1):br(1),tl(2):br(2),:);
%             imshow(pattern)
%             %input('')
%         end
     
    end
end

c={'b','g','r','c','k','m'};
sp={'-o','-*','-x','-.','-+','-s'};
legend_text=cell(3*length(noises_vars),1);
figure(1)
next=1;
for i=1:3
    
    for nidx=1:length(noises_vars)
        %sizes_nm(i).nm
        fracs(i).f(nidx,:)=fracs(i).f(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range1,fracs(i).f(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend_text=cell(3*length(noises_vars),1);
figure(2)
next=1;
for i=4:6    
    for nidx=1:length(noises_vars)
        fracs(i).f(nidx,:)=fracs(i).f(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range2,fracs(i).f(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend_text=cell(3*length(noises_vars),1);
figure(3)
next=1;
for i=1:3
    for nidx=1:length(noises_vars)
        fracs(i).run_time(nidx,:)=fracs(i).run_time(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range1,fracs(i).run_time(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend_text=cell(3*length(noises_vars),1);
figure(4)
next=1;
for i=4:6

    for nidx=1:length(noises_vars)
        fracs(i).run_time(nidx,:)=fracs(i).run_time(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range2,fracs(i).run_time(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


legend_text=cell(3*length(noises_vars),1);
figure(5)
next=1;
for i=1:3

    for nidx=1:length(noises_vars)
        fracs(i).pass(nidx,:)
        sizes_nm(i).nm
        '  '
        fracs(i).pass(nidx,:)=fracs(i).pass(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range1,fracs(i).pass(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend_text=cell(3*length(noises_vars),1);
figure(6)
next=1;
for i=4:6

    for nidx=1:length(noises_vars)
        fracs(i).pass(nidx,:)
        sizes_nm(i).nm
        '  '
        fracs(i).pass(nidx,:)=fracs(i).pass(nidx,:)./sizes_nm(i).nm;    
        hold on
        plot(range2,fracs(i).pass(nidx,:),[c{nidx},sp{i}]);
        legend_text{next}=['size ',num2str(sizes(i)),' var ',num2str(noises_vars(nidx))];
        next=next+1;
    end
end
hold off
legend(legend_text);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






