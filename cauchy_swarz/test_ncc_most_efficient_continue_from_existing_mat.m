function [ res ] = test_ncc_most_efficient_continue_from_existing_mat(filename,points,filters_boxes,thresholds,noises_vars,max_iter,res)


npoints=length(points);

%sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
next=0;

sizes=[16,32,48,64,96,128];
% res=struct;
% for i=1:length(points)
%     for j=1:20
%         res(i,j).pass=zeros(length(noises_vars),1);
%         res(i,j).boxes_time=zeros(length(noises_vars),1);
%         res(i,j).matlab_time=zeros(length(noises_vars),1);
%     end
% end

i=1;
found=0;
while(i<length(points)+1)
    k=1;
    while(k<21)
        if (sum(res(i,k).matlab_time)==0)
            found=1;
            break;
        end
        k=k+1;
    end
    if found==1
        break;
    end
    i=i+1;
end

if i==(length(points)+1)
    return
end
i 
ifrom=i;


for i=ifrom:npoints
    
    %   if sum(sd==i)>0
    %       continue;
    %   end
    
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);
    
    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(points(i).pats)
        
        next=next+1; 
        sz=points(i).pats(k).sz;
%         if sz(1) < 128
%             continue
%         end
        disp([num2str(i),' ',num2str(k),' '])
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
         
%         t4=tic;
%         ncc_map2 = normxcorr2(pattern, I);
%         toc_t4=toc(t4);
%         ncc_map2=ncc_map2(sz(1):end-sz(1)+1,sz(2):end-sz(2)+1);        
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
   
        for nidx=1:length(noises_vars)
        % for nidx=length(noises_vars)
            J=imread(['noise_var',num2str(noises_vars(nidx)),'/',points(i).im_name]);
            J=mat2gray(J);
            
            t3=tic;
            %         %%%%%%% MATLABS NCC %%%%%%%%%%%%%%%%%%
            wins_sums= box_corr(J,[1,n_p,1,m_p],1,n_p,m_p);
            wins_mu=wins_sums/(n_p*m_p);
            wins_norm2 = box_corr(J.^2,[1,n_p,1,m_p],1,n_p,m_p);
            norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
            cor=conv2(J,rot90(pattern,2),'valid');
            ncc_map1=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
            %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc_t3=toc(t3);
             res(i,k).matlab_time(nidx)=toc_t3;
        
            thresh=thresholds(nidx,find(sizes==sz(1)));
            iters=max_iter(sizes==sz(1));
            
            t1=tic;
            [U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,thresh,residuals,iters,pattern,1);
            toc_t1=toc(t1);
            
            t2=tic;
            [U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(J,n_p,m_p,box_arr,w_arr,thresh,residuals,iters,pattern,1);
            toc_t2=toc(t2);
            
            res(i,k).boxes_time(nidx)=min([toc_t1,toc_t2]);
                    
            if(ismember(tl,[ith,jth],'rows'))
                %disp('GOOD')
                res(i,k).pass(nidx)=1;
            else
                %disp('BADDD')
            end
            %disp([num2str(next),' my time: ',num2str(toc_t01),' MATLAB time: ',num2str(toc_t3),' MATLAB time2: ',num2str(toc_t4)]);
            
            %[min([toc_t1,toc_t2]) sz(1) nidx iters thresh max(max(vals)); toc_t3 sz(1) nidx frac_for_direct thresh  max(max(ncc_map1)); toc_t4 sz(1) nidx frac_for_direct thresh max(max(ncc_map2)) ]
            %[v,idxv]=max(vals);
            %[mati1,matj1]=find(ncc_map1==max(max(ncc_map1)),1,'first');
            %[mati2,matj2]=find(ncc_map2==max(max(ncc_map2)),1,'first');
            %[ith(idxv),jth(idxv); mati1 , matj1 ; mati2 , matj2 ; tl(1),tl(2)]
            %input('')
            
        end
      
         save(filename,'res','noises_vars','points','filters_boxes','thresholds','max_iter');
       
        
    end
    %input('')
end
save(['backup_',filename],'res','noises_vars','points','filters_boxes','thresholds','max_iter');
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

