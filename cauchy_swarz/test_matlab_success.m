function [ res ] = test_matlab_success(filename,points,noises_vars)


npoints=length(points);

%sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
next=0;

sizes=[16,32,48,64,96,128];
res=struct;
for i=1:length(points)
    for j=1:20
       
        res(i,j).fft_pass=zeros(length(noises_vars),1);
    end
end
    
for i=1:npoints
    
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
       
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        
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
            ncc_map1(norms_of_wins_minus_mu<0.001)=0;
            %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc_t3=toc(t3);
            
             fft_max = max(ncc_map1(:));
             [fft_max_i, fft_max_j]=find(ncc_map1==fft_max);
             if(ismember(tl,[fft_max_i, fft_max_j] ,'rows') )
                 res(i,k).fft_pass(nidx)=1;
                % tl
                % [fft_max_i, fft_max_j]
                 %disp('mem')
             else
                 disp(['not mem ',num2str(nidx)])
             end
        end
         save(filename,'res','noises_vars','points');
    end
    %input('')
end
save(['backup_',filename],'res','noises_vars','points');
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

