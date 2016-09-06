function [ res ] = calibrate_ncc_thresh2(filename,points,filters_boxes1, thresholds, noises_vars,iter,rand_pats,rand_points)

thresholds=sort(thresholds);
res=struct;
sizes=[16,32,48,64,96,128];


for i=1:length(sizes)
    res(i).sz=sizes(i);
    res(i).nm=zeros(length(noises_vars),length(thresholds));
    res(i).nm_total=0;
end

%rand_points=randperm(length(points),20)
%rand_pats=randperm(20,10)
for i=1:length(rand_points)
    
    I=rgb2gray(imread(points(rand_points(i)).im_name));
    I=mat2gray(I);
    
    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(rand_pats)
        
        sz=points(rand_points(i)).pats(rand_pats(k)).sz;
        res(sizes==sz(1)).nm_total=res(sizes==sz(1)).nm_total+1;
        
        disp(['im ',num2str(rand_points(i)),' k ',num2str(rand_pats(k))]);
        %[i,k,sz(1)]
        tl=points(rand_points(i)).pats(rand_pats(k)).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        
        [n_p,m_p]=size(pattern);
        boxes=filters_boxes1(rand_points(i),rand_pats(k)).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        box_arr=box_arr(2:end,:);
        w_arr=w_arr(2:end);
        
        
        frac_for_direct=0;
        
        szidx=find(sizes==sz(1));
        
        max_iter=iter;
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
        %         t01=tic;
        
        for noiseidx=1:length(noises_vars)
            J=imnoise(I,'gaussian',0,noises_vars(noiseidx));
            for thidx=1:length(thresholds)
                %        [Lm,ith,jth,ind_left]=ncc_match_cauchy_schwarz_most_efficient(J,n_p,m_p,box_arr,w_arr,thresholds(thidx),residuals,max_iter,frac_for_direct);
                [U,ith,jth,vals] = ncc_match_cauchy_with_mex_no_direct_conv(J,n_p,m_p,box_arr,w_arr,thresholds(thidx),residuals,max_iter,pattern,frac_for_direct);
                
                if(ismember(tl,[ith,jth],'rows'))
                    %disp('mem')
                    res(szidx).nm(noiseidx,thidx)=res(szidx).nm(noiseidx,thidx)+1;
                else
                    %disp('not mem')
                    break;
                end
                
            end
        end
        
        
        save(filename,'res','thresholds','noises_vars','rand_pats','rand_points','iter');
        
    end
    %input('')
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

function [b,w]=sort_boxes_by_size(boxes,w_arr)

s=zeros(length(w_arr),1);
for i=1:length(w_arr)
    s(i)=(boxes(i,2)-boxes(i,1)+1)*(boxes(i,4)-boxes(i,3)+1);
end

[s,ix]=sort(s,'descend');

b=boxes(ix,:);

w=w_arr(ix);
end


