function [ measures ] = calibrate_ncc_frac_for_direct2(filename,points,filters_boxes, noises_vars,rand_pats,rand_points,thresholds,range1,range2)

measures=struct;

sizes=[16,32,48,64,96,128];
% for i=1:length(rand_points)
%     for k=1:length(rand_pats)
%         measures(i,k).fracs=zeros(length(noises_vars),iters_num);
%         measures(i,k).run_time=zeros(length(noises_vars),iters_num);
%         measures(i,k).dconv_run_time=zeros(length(noises_vars),1);
%         measures(i,k).pass=zeros(4,1);
%     end
% end

%rand_points=randperm(length(points),20)
%rand_pats=randperm(20,10)
for i=1:length(rand_points)
    
    I=rgb2gray(imread(points(rand_points(i)).im_name));
    I=mat2gray(I);
    
    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(rand_pats)
        
        
        sz=points(rand_points(i)).pats(rand_pats(k)).sz;
        
        
        disp(['im ',num2str(rand_points(i)),' k ',num2str(rand_pats(k))]);
        %[i,k,sz(1)]
        tl=points(rand_points(i)).pats(rand_pats(k)).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        %imshow(pattern)
        [n_p,m_p]=size(pattern);
        boxes=filters_boxes(rand_points(i),rand_pats(k)).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        boxnum=length(w_arr);
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        box_arr=box_arr(2:end,:);
        w_arr=w_arr(2:end);
        
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
        
        
        for nidx=1:length(noises_vars)
            J=imread(['noise_var',num2str(noises_vars(nidx)),'/',points(rand_points(i)).im_name]);
            J=mat2gray(J);
            
            range=range2;
            if sz(1)<64
                range=range1;
            end
            
            for iter_idx=1:length(range)
                if boxnum>=range(iter_idx)
                    t1=tic;
                    [U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,thresholds(nidx,find(sizes==sz(1))),residuals,range(iter_idx),pattern,1);
                    toc_t1=toc(t1);
                    t2=tic;
                    [U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,thresholds(nidx,find(sizes==sz(1))),residuals,range(iter_idx),pattern,1);
                    toc_t2=toc(t2);
                    toc_t=min([toc_t1,toc_t2]);
                else
                    ith=0;
                    jth=0;
                    toc_t=0;
                end
                %toc_t
                %length(ith)
                %iter_idx
                measures(i,k).fracs(nidx,iter_idx)=length(ith);
                measures(i,k).run_time(nidx,iter_idx)=toc_t;
                measures(i,k).iters_num(nidx,iter_idx)=range(iter_idx);
               
                if(ismember(tl,[ith,jth],'rows'))
                    %if(ismember(sub2ind(size(I)-size(pattern)+1,tl(1),tl(2)),lith))
                    measures(i,k).pass(nidx,iter_idx)=1;
                    %disp('mem');
                else
                    measures(i,k).pass(nidx,iter_idx)=0;
                end
                %input('')
            end
            
        end
        
        
        save(filename,'measures','noises_vars','rand_pats','rand_points','thresholds','range1','range2');
        
    end
    %input('')
end
%x
end

function t=get_thresh1(sz,sgma)

k0=[0.7340    0.5194    0.3769    0.2923    0.1999    0.1513];
[szn,sznidx]=nearest_nighbour(sz,sizes);

k=k0(sznidx)*(sgma/0.05);
t=k*sqrt(1+(szn*sgma)^2);

end

function t=get_thresh2(sz,sgma)


[szn,sznidx]=nearest_nighbour(sz,sizes);

k=k0(sznidx)*(sgma/0.05);
t=k*sqrt(1+(szn*sgma)^2);

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

function t  = best_iter_num_binary_search(J,n_p,m_p,box_arr,w_arr,threshold,residuals,max_iter,pattern)

iters=1:maxt_iter;



t1=tic
[U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,threshold,residuals,max_iter,pattern,1);
upper_bound_time=toc(t1)
upper_bound_iter=max_iter;

t2=tic
[U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,threshold,residuals,1,pattern,1);
lower_bound_time=toc(t2)
lower_bound_iter=1;





[U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,threshold,residuals,floor(upper_bound_iter/2),pattern,1);

cur_iter=floor(max_iter/2);





end