function [  ] = test_l2_match_with_cauchy_schwarz(points,filters_boxes, debug,threshold)
%d=dir('640x960_images');

%filters_boxes: [30x10 struct]
%points: [1x30 struct]
%filters_indices: [30x10 double]


npoints=length(points);

count=0;
ac_t1=0;
ac_t2=0;
ac_t3=0;

for i=1:npoints
    
    
    I=rgb2gray(imread(points(i).im_name));
    
    I(1,1:10)
    I=mat2gray(I);
%     if debug
%         figure(996), imshow(I), title('the image')
%     end
    
    
    
    
    for k=1:length(points(i).pats)
        
        
        sz=points(i).pats(k).sz;
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        
        % disp(['-----------------------computing for ',num2str(i),'th image: ',points(i).im_name,'----------']);
        % disp(['-----------------------',num2str(k),'th filter ----------']);
        
        %disp(' ');
        
        %[pattern,tl] = get_filter_from_surf(I,surfs(k), [12 12]);
        
        [n_p,m_p]=size(pattern);
        %numel(I(:))
%         if debug
%             figure(995), imshow(pattern), title('the pattern')
%         end
%         
        
        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        
        
           
        
        
        
        %t0=tic;
        max_iter=20;
        %[Lm1,ind_left1,ith1,jth1,Lr1,Lt1]=l2_match_with_cauchy_schwarz_bound1(I,pattern,n_p,m_p,box_arr(2:end,:),w_arr(2:end),threshold,tl);
        [Lm1,ind_left1,ith1,jth1,Lr1,Lt1,Lrf,Lrs,Ltf,Lts]=l2_match_with_cauchy_schwarz_bound1(I,pattern,n_p,m_p,box_arr(2:end,:),w_arr(2:end,:),threshold,tl,max_iter);
        
        percent_from_before=[ind_left1,ind_left1(end)]./[numel(Lm1),ind_left1];
        %[sorted,sortedidx]=sort(percent_from_before);
        localidx=local_minimas(percent_from_before);
        localidx(1:min([10,length(localidx)]))
       % figure(20)
        %plot(1:length(percent_from_before),percent_from_before);
        %plot(1:length(ind_left),ind_left/numel(L));
        
        figure(35)
        plot(1:length(Lt1),Lr1);
        hold on
        plot(1:length(Lt1),Lt1,'--gs','LineWidth',2);
        hold on
        plot(1:length(Lt1),threshold*n_p*m_p*ones(1,length(Lt1)),'--','LineWidth',2);
        l1='';
        title(['Lower bound random indexes, true index (dashed line = threshold: ',l1]);
        hold off
        
        figure(40)
        
        plot(1:length(Lrf),Lrf,'b');
        hold on
        plot(1:length(Lrf),Lrs,'r');
        hold on
        plot(1:length(Lrf),Ltf,'b--');
        hold on
        plot(1:length(Lrf),Lts,'r--');
        hold on
        plot(1:length(Lrf),-2*(Ltf+Lts),'g--');
        hold on
        
        plot(1:length(Lrf),-2*(Lrf+Lrs),'g');
        hold on
        
        hold on
        plot(1:length(Lt1),threshold*n_p*m_p*ones(1,length(Lt1)),'--','LineWidth',2);
%         plot(1:length(Lt1),[0,cumsum(bc2(1:length(Lt1)-1))],'--rs','LineWidth',2);
%         hold on
        
        
        legend('rand first','rand sec','true first','true sec','-sum(true first,true sec)','-sum(rand first,rand sec)');
        hold off
        
%         %%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%
%          boxes2=filters_boxes2(i,k).boxes;
%         box_arr2=boxes2(:,1:4);
%         w_arr2=boxes2(:,5);
%         [Lm2,ind_left2,ith2,jth2,Lr2,Lt2]=l2_match_with_cauchy_schwarz_bound1(I,pattern,n_p,m_p,box_arr2(2:end,:),w_arr2(2:end),threshold,tl);
%         percent_from_before=[ind_left2,ind_left2(end)]./[numel(Lm2),ind_left2];
%         figure(36)
%         plot(1:length(Lt2),Lr2);
%         hold on
%         plot(1:length(Lt2),Lt2,'--gs','LineWidth',2);
%         hold on
%         plot(1:length(Lt2),threshold*n_p*m_p*ones(1,length(Lt2)),'--','LineWidth',2);
%         title(['Lower bound random indexes, true index (dashed line = threshold: ',l2]);
%         hold off
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%          figure(20)
%         %plot(1:length(percent_from_before),percent_from_before);
%         plot(1:length(ind_left1),[ind_left1 ; ind_left2]/numel(Lm1));
%         axis([1 length(ind_left1) 0 1])
%        legend(l1,l2)
       
       figure(20)
        %plot(1:length(percent_from_before),percent_from_before);
        plot(1:length(ind_left1),ind_left1/numel(Lm1));
        axis([1 length(ind_left1) 0 1])
       
        
%        figure(21)
%        plot(1:length(w_arr),[w_arr' ; w_arr2']);
%        title('coefficients')
%        legend(l1,l2)
%        
       figure(21)
       plot(1:length(w_arr),w_arr');
       title('coefficients')
  
       
%         bnum=min([7,length(localidx)]);
         %rf1=reconstruct_filter(box_arr1(localidx(1:bnum),:),w_arr1(localidx(1:bnum),:),n_p,m_p,0);
         rf1=reconstruct_filter(box_arr,w_arr,n_p,m_p,0);
        % rf2=reconstruct_filter(box_arr2,w_arr2,n_p,m_p,0);
%         figure(97)
%         imshow(mat2gray(rf))
        if debug
            figure(995), imshow([pattern,rf1]), title('the pattern')
        end
         figure(30)
         imagesc([Lm1])
        input('')
        %close all
        %toc_t0=toc(t0)
        
        %input('')
        
        %t01=tic;
        %[b,w]=sort_boxes_by_size(box_arr,w_arr);
        %Lm=l2_match_with_cauchy_schwarz_bound1(single(I),filter,win_sz(k,1),win_sz(k,2),box_arr,single(w_arr),single(threshold) ,single(reconst_filter),debug,exact);
        %[left_ind2,Lm]=l2_match_with_cauchy_schwarz_bound1(I,pattern,win_sz(k,1),win_sz(k,2),b,w,threshold ,reconst_filter,debug,exact);
        
        %toc_t01=toc(t0)
        
        %plot(1:length(left_ind1),[left_ind1/numel(L) ; left_ind2]);
        %legend('orig','sorted');
        %input('')
        %[b,w]=sort_boxes_by_size(boxes,w_arr)
        
        
        %'done Lm'
        
        
        %         t1=tic;
        %
        %         %Lc=cauchy_swarz_l2_match(I,int32(box_arr),w_arr,threshold,reconst_filter);
        %         I(1:10,1)
        %         L=mex_cauchy_l2(single(I),int32(box_arr),single(w_arr),single(threshold),single(reconst_filter),[1]);
        %         L(1:10,1:10)
        %         L = double(L);
        %
        %         %  Lc=cauchy_swarz_l2_match(I,int32(box_arr),single(w_arr),single(threshold),single(reconst_filter),[]);
        %         toc_t1=toc(t1);
        %         ac_t1=ac_t1+toc_t1;
        %
        %         t2=tic;
        %         %L2=mex_cauchy_l2(single(I),int32(box_arr),single(w_arr),single(threshold),single(reconst_filter),[]);
        %         %disp(' EEEEEEEEEEEE ')
        %         %L2(1:50,1:50)
        %         %L2 = double(L2);
        %         toc_t2=toc(t2);
        %         ac_t2=ac_t2+toc_t2;
        %
        %          t3=tic;
        %         %%%%%%%%%%%%%% MATLAB'S CONV2 %%%%%%%%%%%%%%%
        %          pat_norm2 = sum(pattern(:).^2);
        %          wins_norm2 = box_corr_c(I.^2,int32([1,n_p,1,m_p]),1,n_p,m_p);
        %          corr_surf = conv2(I,rot90(pattern,2),'valid');
        %          Lconv2= pat_norm2 + wins_norm2 - 2*corr_surf;
        % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         toc_t3=toc(t3);
        %         ac_t3=ac_t3+toc_t3;
        %
        %         [toc_t1 toc_t2 toc_t3]
        %
        %         sum(sum(abs(Lm-L)))
        %         count=count+1;
        %          %imagesc([L,L2,Lm])
        % %
        % %         [sum(sum(Lm-Lc)) sum(sum(Lm-Lopt)) sum(sum(Lc-Lopt))]
        %          disp('idx ptr fft');
        %
        %          %[ac_t1 ac_t2 ac_t3]/count
        % %         toc_t3
        % %
        % %
        % %
        %          mm=min(min(Lm));
        %          ml=min(min(L));
        %          ml2=min(min(L2));
        %          mmat=min(min(Lconv2));
        %
        %          %mc=min(min(Lc));
        %          %mo=min(min(Lopt));
        % %         m2=min(min(Lconv2));
        % %
        %          [mi mj]=find(Lm==mm,1);
        %          [li lj]=find(L==ml,1);
        %          [l2i l2j]=find(L2==ml2,1);
        %          [mati matj]=find(Lconv2==mmat,1);
        % %
        %          [ tl(1) mi li l2i mati ; tl(2) mj lj l2j matj]
        %          %[ tl(1) mi ci oi; tl(2) mj cj oj]
        %          [mm ml ml2 mmat]
        
        %         x=try_mex_with_opencv(single(I),single(filter));
        %         imshow(x)
        %          input('');
        
        %cum_t=cum_t+toc_t1;
        %disp(['avarage time to cauchy swarz: ',num2str(cum_t/k)]);
        
        %         if norm([fi,fj]-tl)<=5
        %             figure(1),hold on, plot(arr), title(['FOUND! threshold=',num2str(threshold),' distance=',num2str(m)])
        %         else
        %             figure(1),hold on, plot(arr), title(['NOT FOUND! threshold=',num2str(threshold),' distance=',num2str(m)])
        %         end
    end
    input('')
end
%x
end

function [filter,tl] = get_filter_from_surf(I,s, std_sz)

win_sz=[ round(s.Scale*std_sz(1)) , round(s.Scale*std_sz(2)) ];
center=[s.Location(:,2) s.Location(:,1)];

tl=round(center)-floor(win_sz/2);
br=tl+win_sz-1;

filter=I(tl(1):br(1),tl(2):br(2) );
end

function [b,w]=sort_boxes_by_size(boxes,w_arr)

s=zeros(length(w_arr),1);
for i=1:length(w_arr)
    s(i)=(boxes(i,2)-boxes(i,1)+1)*(boxes(i,4)-boxes(i,3)+1);
end

[s,ix]=sort(s,'descend');

b=boxes(ix,:);

w=w_arr(ix);
end


