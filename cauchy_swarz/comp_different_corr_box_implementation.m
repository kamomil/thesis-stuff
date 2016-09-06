function [  ] = comp_different_corr_box_implementation(noises_vars)

boxes_measures  = struct;
matlab_measures = struct;
%d=dir('640x960_images');

%filters_boxes: [30x10 struct]
%points: [1x30 struct]
%filters_indices: [30x10 double]
load('boxes_up_to_5per_for_filters_std40_uptosize100.mat');

npoints=length(points);
nfilters=size(filters_boxes,2);
for i=1:npoints
    
    surfs=points(i).surfs;
    I=mat2gray(rgb2gray(imread(points(i).im_name)));
    %figure(996), imshow(I), title('the image')
    win_sz=round(surfs.Scale*[40 40]);
    %centers=[surfs.Location(:,2) surfs.Location(:,1)];
    
    for k=1:nfilters
        
        idx=filters_indices(i,k);
        disp(['-----------------------computing for ',num2str(i),'th image: ',points(i).im_name,'----------']);
        disp(['-----------------------',num2str(k),'th filter is at index: ',num2str(idx),'----------']);
        
        disp(' ');
        
        [filter] = get_filter_from_surf(I,surfs(idx), [40 40]);
        
        %figure(995), imshow(filter), title('the pattern')
        
        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        
        [n_p,m_p]=size(filter);
        
        int_img=integral_image(I);
        
        b1=0;
        b2=0;
        b3=0;
        for compidx=1:1000
            
            ind=rand(size(I)-size(filter)+1);
            
            ind=ind>0.5;
            t1=tic;
            C1=box_corr_given_indices_slow(int_img,box_arr,w_arr,n_p,m_p,ind);
            b1=b1+toc(t1);
            
            t2=tic;
            C2=box_corr_given_indices2(int_img,box_arr,w_arr,n_p,m_p,ind);
            b2=b2+toc(t2);
            
            t3=tic;
            C3=box_corr_given_int_img(int_img,box_arr,w_arr,n_p,m_p);
            b3=b3+toc(t3);
            
            norm(C1-C2,'fro') %correct!
            
            disp(['avarage time to slow: ',num2str(b1/compidx)]);
            disp(['avarage time to 2: ',num2str(b2/compidx)]);
            disp(['avarage time to full box_corr: ',num2str(b3/compidx)]);
           
        end
    end
end
end

function [filter] = get_filter_from_surf(I,s, std_sz)
 
    win_sz=[ round(s.Scale*std_sz(1)) , round(s.Scale*std_sz(2)) ];
    center=[s.Location(:,2) s.Location(:,1)];
    
    tl=round(center)-floor(win_sz/2);
    br=tl+win_sz-1;
    
    filter=I(tl(1):br(1),tl(2):br(2) );
end


