function [  ] = test_ncc_match_with_cauchy_schwarz(points,filters_boxes, debug,threshold)
%d=dir('640x960_images');

%filters_boxes: [30x10 struct]
%points: [1x30 struct]
%filters_indices: [30x10 double]


npoints=length(points);

sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
for i=1:npoints
    
     if sum(sd==i)>0
        continue;
    end
    
    I=rgb2gray(imread(points(i).im_name));
    
    I=mat2gray(I);
%     figure
%     imshow(I<0.05)
    %     if debug
    %         figure(996), imshow(I), title('the image')
    %     end
    
    sigma = 0.2;    
    I=imnoise(I,'gaussian',0,sigma^2);
    
    
    for k=1:length(points(i).pats)
        
        [i,k]
        
        sz=points(i).pats(k).sz;
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        
        
        [n_p,m_p]=size(pattern);
        
        
        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        box_arr=box_arr(2:end,:);
        w_arr=w_arr(2:end);
        
        normp=norm(pattern,'fro');
        max_iter=15;
        
        
%         tione = sub2ind(size(I)-size(pattern)+1, tl(1),tl(2));
%         if tione>5
%             given_indexes=((tione-4):(tione-4))';
%         else
%             given_indexes=((tione+4):(tione+4))';
%         end
        
        gidx=[1,1]; 
    
%         if(i==1 && k==5)
%             gidx=[15,475];
%         end
%         
%         if(i==1 && k==6)
%             gidx=[61,709];
%         end
%         if(i==3 && k==2)
%             gidx=[425,353];
%         end
%         if(i==3 && k==5)
%             gidx=[404,400];
%         end
        given_indexes=(sub2ind(size(I)-size(pattern)+1, gidx(1),gidx(2)))';
        frac_for_direct=0.001;
         residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
         
         J=(rand(1))*I;
         J(J>1)=1;
         figure(120)
         imshow([I ; J])
                                           frac_to_direct=0.0001;                                                          
       [U,i,j,ind_left]=ncc_match_cauchy_schwarz_most_efficient(   J,n_p,m_p,box_arr,w_arr,threshold,residuals,max_iter,frac_to_direct);
        [Lm2,ind_left2,ith2,jth2,Lr2,Lg2,Lt2,Urf2,Urs2,Utf2,Uts2]=ncc_match_with_cauchy_schwarz_bound1(I,max_iter,n_p,m_p,box_arr,w_arr,threshold,residuals,tl,given_indexes);
        
%         [U,ith2,jth2,vals,ind_num] = ncc_match_cauchy_with_mex(I,n_p,m_p,box_arr,w_arr,threshold,residuals,max_iter,pattern,frac_for_direct);
%         if(ismember(tl,[ith2,jth2],'rows'))
%             'GOOD'
%         end

plot(ind_left/numel(U))
            
        %[Lm,ind_left,ith,jth]=l2_match_with_cauchy_schwarz_bound1(0.5*ones(size(I)),pattern,n_p,m_p,box_arr,w_arr,threshold);
        %percent_from_before=[ind_left1,ind_left1(end)]./[numel(Lm1),ind_left1];
        %[sorted,sortedidx]=sort(percent_from_before);
        %length(Lt1)
        %length(Lr1)
        
        bc2=get_Bc2(box_arr,w_arr);
        %disp('done 1')
        figure(35)
        plot(1:length(Lt1),Lr1);
        hold on
        plot(1:length(Lt1),Lt1,'--gs','LineWidth',2);
        hold on
        plot(1:length(Lt1),Lt2,'--rs','LineWidth',2);
        hold on
        
%         plot(1:length(Lt1),Lg1,'--rs','LineWidth',2);
%         hold on
        plot(1:length(Lt1),threshold*ones(1,length(Lt1)),'--','LineWidth',2);
        title(['Lower bound random indexes, true index (dashed line = threshold: ']);
        hold off
        
        %disp('done 1')
        figure(40)
        
        plot(1:length(Urf),Urf,'b');
        hold on
        plot(1:length(Urf),Urs,'r');
        hold on
        plot(1:length(Urf),Utf,'b--');
        hold on
        plot(1:length(Urf),Uts,'r--');
        hold on
        plot(1:length(Urf),Utf-Uts,'g--');
        hold on
        plot(1:length(Lt1),[0,cumsum(bc2(1:length(Lt1)-1))],'--rs','LineWidth',2);
        hold on
        
        
        legend('rand first','rand sec','true first','true sec','diff(true first,true sec)','sum_i |Bi|ci^2');
        hold off
        
        
        figure(20)
        %plot(1:length(percent_from_before),percent_from_before);
        %plot(1:length(ind_left1),ind_left1/numel(Lm1));
        plot(1:max_iter,ind_left1(1:max_iter)/numel(Lm1));
        %axis([1 length(ind_left1) 0 1])
        axis([1 max_iter 0 1])
        %legend(l1)
        
%         figure(21)
%         plot(1:length(w_arr(2:10)),w_arr(2:10)' );
%         title('coefficients')
        %legend(l1)
        
%         figure(345)
%         imshow(plot_boxes(box_arr(1:10,:),w_arr(1:10),n_p,m_p));
%         %         bnum=min([7,length(localidx)]);
%         %rf1=reconstruct_filter(box_arr1(localidx(1:bnum),:),w_arr1(localidx(1:bnum),:),n_p,m_p,0);
%         rf1=reconstruct_filter(box_arr,w_arr,n_p,m_p,0);
%         %rf2=reconstruct_filter(box_arr2,w_arr2,n_p,m_p,0);
%         %         figure(97)
%         %         imshow(mat2gray(rf))
%         if debug
%             %figure(995), imshow([pattern,rf1,I(gidx(1):(gidx(1)+sz(1)-1),gidx(2):(gidx(2)+sz(2)-1))]), title('the pattern')
%             
%             figure(995), imshow([pattern,rf1]),title('the pattern')
%         end
       input('')
    end
    input('')
end
%x
end

function [r] = get_residual_pat_norms(pat,box_arr,w_arr)

r=zeros(1,length(w_arr)+1);
for i=1:length(w_arr)
    r(i)=norm(pat,'fro');
    
    f_r=box_arr(i,1);
    t_r=box_arr(i,2);
    f_c=box_arr(i,3);
    t_c=box_arr(i,4);
    
    pat(f_r:t_r,f_c:t_c)=pat(f_r:t_r,f_c:t_c)-w_arr(i);
end
r(end)=norm(pat,'fro');

end

function [r] = get_Bc2(box_arr,w_arr)

r=zeros(1,length(w_arr));
for i=1:length(w_arr)
    
    
    f_r=box_arr(i,1);
    t_r=box_arr(i,2);
    f_c=box_arr(i,3);
    t_c=box_arr(i,4);
    
    r(i)=(t_r-f_r+1)*(t_c-f_c+1)*w_arr(i)*w_arr(i);
 
end


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


