function [ res ] = calibrate_ncc_thresh(filename,points,filters_boxes1, thresholds, noises_vars,iters)


npoints=length(points);

sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
next=0;
res=struct;
sizes=[16,32,48,64,96,128];

for k=1:length(iters)
    for i=1:length(sizes)
        res(i,k).sz=sizes(i);
        res(i,k).ind_av_per=zeros(length(noises_vars),length(thresholds));
        res(i,k).nm=zeros(length(noises_vars),length(thresholds));
    end
end

for i=1:npoints
    
    if sum(sd==i)>0
        continue;
    end
   
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);

    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(points(i).pats)
        
        next=next+1;         
        
        sz=points(i).pats(k).sz;
        disp(['im ',num2str(i),' k ',num2str(k)]);
        %[i,k,sz(1)]
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        
        [n_p,m_p]=size(pattern);          
        boxes=filters_boxes1(i,k).boxes;
        box_arr=boxes(:,1:4);
        w_arr=boxes(:,5);
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        box_arr=box_arr(2:end,:);
        w_arr=w_arr(2:end);
         
        
        frac_for_direct=0;
        
        szidx=find(sizes==sz(1));
        
       
%         t01=tic;
        for thidx=1:length(thresholds)
            for noiseidx=1:length(noises_vars)
                for itersidx=1:length(iters)
                    max_iter=iters(itersidx)+1;
                    J=imnoise(I,'gaussian',0,noises_vars(noiseidx));
                    [Lm,ith,jth,ind_left]=ncc_match_cauchy_schwarz_most_efficient(J,n_p,m_p,box_arr,w_arr,thresholds(thidx),get_residual_pat_norms(pattern,box_arr,w_arr),max_iter,frac_for_direct);       
                    if(ismember(tl,[ith,jth],'rows'))
                        res(szidx,itersidx).nm(noiseidx,thidx)=res(szidx,itersidx).nm(noiseidx,thidx)+1;   
                        res(szidx,itersidx).ind_av_per(noiseidx,thidx)=res(szidx,itersidx).ind_av_per(noiseidx,thidx)+ind_left(iters(itersidx))/numel(Lm);
                      
                    end
                end
            end
        end

     
       save(filename,'res','thresholds','noises_vars','iters');
      
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


