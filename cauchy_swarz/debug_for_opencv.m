function [ res ] = debug_for_opencv(filename,points,filters_boxes,thresholds,noises_vars,max_iter,frac_for_direct)

npoints=length(points);

for i=1:npoints
    
    %   if sum(sd==i)>0
    %       continue;
    %   end
    
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);
    
    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(points(i).pats)
        
        disp(num2str([i,k]));
        
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
        
        residuals=get_residual_pat_norms(pattern,box_arr,w_arr);
        
        %residuals(1:10)
        
        %sub2ind(size(I),tl(1),tl(2))
        
        thresh=0.9;
        iters=100;
        frac=0.01;
        
        %t1=tic;
        [U,ith,jth,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(I,n_p,m_p,box_arr,w_arr,thresh,residuals,iters,pattern,frac);
        %toc_t1=toc(t1);
        
        C=normxcorr2(pattern,I);
        C=C(sz(1):end-sz(1)+1,sz(2):end-sz(2)+1);
        %         %%%%%%% MATLABS NCC %%%%%%%%%%%%%%%%%%
        int_I=integral_image(I);
        int_I2=integral_image(I.^2);
        
        wins_sums= box_corr(int_I,[1,n_p,1,m_p],1,n_p,m_p);
        wins_mu=wins_sums/(n_p*m_p);
        wins_norm2 = box_corr(int_I2,[1,n_p,1,m_p],1,n_p,m_p);
        norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
        
        cor=conv2(I,rot90(pattern,2),'valid');
        C2=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
        C2(norms_of_wins_minus_mu<0.01)=0;
        %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [mf1,imf1]=max(C(:));
        [mf2,imf2]=max(C2(:));
        [m,idx]=max(vals);
        orig=sub2ind(size(C),tl(1),tl(2));
        disp(['orig ',num2str(sub2ind(size(C),tl(1),tl(2)))]);
        %disp(['fft1 ',num2str(imf1)]);
        %disp(['fft2 ',num2str(imf2)]);
        disp(['boxs ',num2str(sub2ind(size(U),ith(idx),jth(idx)))]);
        figure(1)
        imagesc(U)
        figure(2)
        imagesc(C)
        if(ismember(tl,[ith,jth],'rows'))
           disp('is member'); 
        end
        if(imf1 ~= orig)
            disp('mf1 differ')
            
        end
        if(imf2 ~= orig)
            disp('mf2 differ')
        end
        
        
        
        
        
        
    end
end
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

