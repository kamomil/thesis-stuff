function [  ] = run_l2_match_with_cauchy_schwarz(points,filters_boxes,debug,threshold,exact)
%d=dir('640x960_images');

%filters_boxes: [30x10 struct]
%points: [1x30 struct]
%filters_indices: [30x10 double]


npoints=length(points);

for i=1%:npoints
    
    I=mat2gray(rgb2gray(imread(points(i).im_name)));
    
    if debug
        figure(996), imshow(I), title('the image')
    end
    cum_t=0;
    for k=1%:length(surfs)
    
        
        disp(' ');
        disp(['-----------------------computing for ',num2str(i),'th image: ',points(i).im_name,'----------']);
        disp(['-----------------------',num2str(k),'th filter ----------']);
        
        
        
        %[filter,tl] = get_filter_from_surf(I,surfs(k), [12 12]);
        
        sz=points(i).pats(k).sz(1);
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
    
        pattern=I(tl(1):br(1),tl(2):br(2));

        if debug
            figure(995), imshow(pattern), title('the pattern')
        end
        
        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        
        %[box_arr,w_arr]=sort_boxes_by_w(box_arr,w_arr);
        reconst_filter=reconstruct_filter(box_arr,w_arr,sz(1),sz(2),0);
        
        
        t1=tic;
        
        
        box_arr_o=box_arr;
        
        w_arr_o=w_arr;
        
        for imgidx=1:npoints
            
            box_arr=box_arr_o;
        
            w_arr=w_arr_o;
        
            figure,
            Inow=mat2gray(rgb2gray(imread(points(imgidx).im_name)));
            %[L,arr1]=l2_match_with_cauchy_schwarz_bound1(I,filter,win_sz(k,1),win_sz(k,2),box_arr,w_arr,threshold ,reconst_filter,debug,1);
            [arr,L]=l2_match_with_cauchy_schwarz_bound1(Inow,pattern,sz(k,1),sz(k,2),box_arr,w_arr,threshold ,reconst_filter,debug,0);
            diff=conv(arr,[-1,1],'valid');
            [~,ix]=sort(diff,'descend');
            ix=[ix,max(ix)+1];
            
            box_arr=box_arr(ix,:);
    
            w_arr=w_arr(ix);
            
            [arr2,L]=l2_match_with_cauchy_schwarz_bound1(Inow,pattern,sz(k,1),sz(k,2),box_arr,w_arr,threshold ,reconst_filter,debug,0);
            size(arr)
            box_sizes=(box_arr(:,2)-box_arr(:,1)+1).*(box_arr(:,4)-box_arr(:,3)+1);
            size(box_sizes)
            
        
            
            
            
            
            %stem(0:length(w_arr), [0,box_sizes'/box_sizes(1)]);
        %set(gca,'XTick',0:1:length(w_arr))
        %hold on
        %stem(1:length(w_arr)-1, w_arr(2:end)','Color','r');
            plot(0:length(w_arr),[[arr/numel(L),0] ; [arr2/numel(L),0]]); 
            %stem(0:length(w_arr),[arr/numel(L),0],'Color','r');
            %set(gca,'XTick',0:1:length(w_arr))
            %hold on
            %stem(0:length(w_arr), [0,w_arr'],'Color','g');
        %plot(1:length(w_arr)+1,[arr/numel(L),0],'Color','r');
        
        
        %reconst_filter=reconstruct_filter(box_arr,w_arr,win_sz(k,1),win_sz(k,2),1);
        input('');
        end
        close all
        %plot(1:length(arr1),[arr1;arr2])
        %legend('exact','app');
        %input('')
        toc_t1=toc(t1);
        cum_t=cum_t+toc_t1;
        disp(['avarage time to cauchy swarz: ',num2str(cum_t/k)]);
        
        % compute template squared norm
        pat_norm2 = sum(pattern(:).^2);
        
        
        
        
        %D = pat_norm2 + wins_norm2 - 2*L;
        
        m=min(min(L));
        [fi,fj]=find(L==m);
        
        if norm([fi,fj]-tl)<=3
            %figure(1),hold on, plot(arr), title(['FOUND! threshold=',num2str(threshold),' distance=',num2str(m)])
        else
            %figure(1),hold on, plot(arr), title(['NOT FOUND! threshold=',num2str(threshold),' distance=',num2str(m)])
        end
    end
end
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

function [b,w]=sort_boxes_by_w(boxes,w_arr)
    
    
    [w,ix]=sort(w_arr);
    
    b=boxes(ix,:);
    
    
end