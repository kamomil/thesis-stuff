function [  ] = boxes_performance(points,filters_boxes,debug,threshold)
%d=dir('640x960_images');

%filters_boxes: [30x10 struct]
%points: [1x30 struct]
%filters_indices: [30x10 double]


npoints=length(points);

for i=2:npoints
    
    pats=points(i).pats;
    I=mat2gray(rgb2gray(imread(points(i).im_name)));
    
    if debug
        figure(996), imshow(I), title('the image')
    end
    
    for k=1:length(pats) 
    
        sz=points(i).pats(k).sz;
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        pattern=I(tl(1):br(1),tl(2):br(2));
        [n_p,m_p]=size(pattern);
      
        if debug
            figure(995), imshow(pattern), title('the pattern')
        end
        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        
        w_arr=w_arr(1:min([50,length(w_arr)]));
        box_arr=box_arr(1:length(w_arr),:);
        
        reconst_filter=reconstruct_filter(box_arr,w_arr,n_p,m_p,0);
        reconst_filter_o=reconst_filter;
        
        box_arr_learned_order=zeros(size(box_arr));
        w_arr_learned_order=zeros(size(w_arr));
        
        [n_p,m_p]=size(pattern);
        
        [L1,arr1,i_indexes,j_indexes]=l2_match_with_cauchy_schwarz_bound1(I,pattern,n_p,m_p,box_arr,w_arr,threshold);
         figure
        plot(0:length(w_arr)-1,arr1/numel(L1) ); 
        %input('')
        for bidx=1:length(w_arr)
           % [length(w_arr),bidx]
           bidx
            arr=zeros(1,length(w_arr));
            
            % run on points 8,..,20
            for imgidx=8:9%12
                imgidx
                Inow=mat2gray(rgb2gray(imread(points(imgidx).im_name)));
                arr=arr+boxes_performance_on_image(Inow,pattern,n_p,m_p,box_arr,w_arr,threshold ,reconst_filter,1,0);
            end
            %arr
            arr=arr/2;
            %figure,
            %stem(arr)
            
            %title(['box idx=',num2str(bidx),' length(arr)=',num2str(length(arr))]);
            %input('')
            %close all
            [M,maxidx]=max(arr);
            disp([num2str(bidx),' ',num2str(maxidx)]);
            
            box_arr_learned_order(bidx,:)=box_arr(maxidx,:);
            w_arr_learned_order(bidx)=w_arr(maxidx);
            
             f=zeros(n_p,m_p);
            f(box_arr(maxidx,1):box_arr(maxidx,2),box_arr(maxidx,3):box_arr(maxidx,4))=w_arr(maxidx);
            reconst_filter=reconst_filter-f;
            
            box_arr(maxidx,:)=zeros(1,4);
            w_arr(maxidx)=0;
            
            
        end
        
        for bb=find(w_arr_learned_order==0,1,'first'):length(w_arr)
            [iw,jw,w]=find(w_arr,1,'first');
            w_arr_learned_order(bb)=w;
            box_arr_learned_order(bb,:)=box_arr(iw,:);
            w_arr(iw)=0;
            box_arr(iw,:)=0;
        end
            
                
        %[w_arr_learned_order,box_arr_learned_order]
        
         [L2,arr2,i_indexes,j_indexes]=l2_match_with_cauchy_schwarz_bound1(I,pattern,n_p,m_p,box_arr_learned_order,w_arr_learned_order,threshold);
        arr2
         figure
        plot(0:length(w_arr)-1,arr2/numel(L1) ); 
        figure
        imagesc([L1,L2]);
        
         m1=min(min(L1));
         m2=min(min(2));
          
         [mi, mj]=find(L1==m1,1);
         [li, lj]=find(L2==m2,1);
          
         [ tl(1), mi, li ; tl(2), mj, lj ]
         
         
        input('');
        close all;  
            
        %         end
    end
    
end
%x
end
