function [] = add_residual_norms_to_boxes_mat_file(boxesfile)



 % filters_boxes: [53x20 struct]
 % points: [1x53 struct]
 % sample_method: 1
    
load(boxesfile)

npoints=length(points);
for i=1:npoints
    
    i
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);
    
    
    
    for k=1:length(points(i).pats)
        

        sz=points(i).pats(k).sz;
        disp([num2str(i),' ',num2str(k),' '])
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        

        boxes=filters_boxes(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);

        filters_boxes(i,k).boxes=[boxes, residuals];
        
    end
    save(boxesfile,'filters_boxes','points','sample_method');

end

end

function [r] = get_residual_pat_norms(pat,box_arr,w_arr)

r=zeros(length(w_arr),1);
for i=1:length(w_arr)
    f_r=box_arr(i,1);
    t_r=box_arr(i,2);
    f_c=box_arr(i,3);
    t_c=box_arr(i,4);
    pat(f_r:t_r,f_c:t_c)=pat(f_r:t_r,f_c:t_c)-w_arr(i);
    r(i)=norm(pat,'fro');
end
end

