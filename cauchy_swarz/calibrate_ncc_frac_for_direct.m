function [ measures ] = calibrate_ncc_frac_for_direct(filename,points,filters_boxes, noises_vars,iters_num,rand_pats,rand_points,thresholds)

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
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        box_arr=box_arr(2:end,:);
        w_arr=w_arr(2:end);
        
        
        frac_for_direct=1;
        
        max_iter=iters_num;
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
        
        
        for nidx=1:length(noises_vars)
            J=imread(['noise_var',num2str(noises_vars(nidx)),'/',points(rand_points(i)).im_name]);
            J=mat2gray(J);
            %        [Lm,ith,jth,ind_left]=ncc_match_cauchy_schwarz_most_efficient(J,n_p,m_p,box_arr,w_arr,thresholds(thidx),residuals,max_iter,frac_for_direct);
            %[U,ith,jth,vals,ind_fracs,run_times,dconv_run_time] = ncc_match_cauchy_with_mex(J,n_p,m_p,box_arr,w_arr,thresholds(noiseidx,find(sizes==sz(1))),residuals,max_iter,pattern,frac_for_direct);
            
               
            %if(linear==1)
                [U,ith,jth,vals,ind_fracs,run_times,dconv_run_time] = linear_idx_ncc_match_cauchy_with_mex(J,n_p,m_p,box_arr,w_arr,thresholds(nidx,find(sizes==sz(1))),residuals,max_iter,pattern,frac_for_direct);
                %ind_fracs
                %else
                %[U,ith,jth,vals,ind_fracs,run_times,dconv_run_time] = ncc_match_cauchy_with_mex(J,n_p,m_p,box_arr,w_arr,thresholds(nidx,find(sizes==sz(1))),residuals,max_iter,pattern,frac_for_direct);
            %run_times2
            %end1
            %[ dconv_run_time1 ; dconv_run_time2]
            %[size(ith1) , size(ith2)]
            %if(sum(size(ith1) ~= size(ith2)))
            % [size(ith1) size(ith2)]
            %    error('BUG!');
            %end
            %if(sum(ith1 ~= ith2)>0)
            %    error('BUG!');
            %end
            %dconv_run_time
            measures(i,k).fracs(nidx,:)=double(ind_fracs);
            measures(i,k).run_time(nidx,:)=run_times;
            measures(i,k).dconv_run_time(nidx,1)=dconv_run_time;
            if(ismember(tl,[ith,jth],'rows'))
                %if(ismember(sub2ind(size(I)-size(pattern)+1,tl(1),tl(2)),lith))
                measures(i,k).pass(nidx)=1;
                %disp('mem');
            else
                disp('ERR');
                [sz(1) nidx thresholds(nidx,find(sizes==sz(1)))]
            end
            %input('')
            
            
        end
        
        
        save(filename,'measures','iters_num','noises_vars','rand_pats','rand_points','thresholds');
        
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

