function [  ] = compare_different_fft_implementations(points)


npoints=length(points);

%sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
next=0;

sizes=[16,32,48,64,96,128];
% res=struct;
% for i=1:length(points)
%     for j=1:20
%         res(i,j).pass=zeros(length(noises_vars),1);
%         res(i,j).boxes_time=zeros(length(noises_vars),1);
%         res(i,j).fracs=zeros(length(noises_vars),1);
%         res(i,j).matlab_time=zeros(length(noises_vars),1);
%     end
% end
    
for i=1:npoints
     
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);
    
    I=imnoise(I,'gaussian',0,0.05);
    
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
        pattern=pattern-mean(pattern(:));
            
            t1=tic;
            %         %%%%%%% MATLABS NCC 1 %%%%%%%%%%%%%%%%%%
            wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
            wins_mu=wins_sums/(n_p*m_p);
            wins_norm2 = box_corr(I.^2,[1,n_p,1,m_p],1,n_p,m_p);
            norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
            cor=conv2(I,rot90(pattern,2),'valid');
            ncc_map1=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
            %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc_t1=toc(t1);
             %res(i,k).matlab_time(nidx)=toc_t1;
        
            
            
            
            t2=tic;
            %         %%%%%%% MATLABS NCC %%%%%%%%%%%%%%%%%%
            %wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
            wins_sums=conv2(I,ones(n_p,m_p),'valid');
            wins_mu=wins_sums/(n_p*m_p);
            
            
            %wins_norm2 = box_corr(I.^2,[1,n_p,1,m_p],1,n_p,m_p);
            wins_norm2= conv2(I.^2,ones(n_p,m_p),'valid');
            norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
            cor=conv2(I,rot90(pattern,2),'valid');
            ncc_map2=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
            %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc_t2=toc(t2);
            
            t3=tic;
            ncc_map3 = normxcorr2(pattern, I);
            ncc_map3=ncc_map3(sz(1):end-sz(1)+1,sz(2):end-sz(2)+1);
            toc_t3=toc(t3);
            
            ncc_map1(1:10,1:10)
            ncc_map3(1:10,1:10)
            
            [toc_t1 toc_t2 toc_t3]
            
            imagesc([ncc_map1, ncc_map2, ncc_map3])
            title(num2str([toc_t1 toc_t2 toc_t3]))
            
            input('')
            
            
            
            
            
            %disp([num2str(next),' my time: ',num2str(toc_t01),' MATLAB time: ',num2str(toc_t3),' MATLAB time2: ',num2str(toc_t4)]);
            
            %[min([toc_t1,toc_t2]) sz(1) nidx iters thresh max(max(vals)); toc_t3 sz(1) nidx frac_for_direct thresh  max(max(ncc_map1)); toc_t4 sz(1) nidx frac_for_direct thresh max(max(ncc_map2)) ]
            %[v,idxv]=max(vals);
            %[mati1,matj1]=find(ncc_map1==max(max(ncc_map1)),1,'first');
            %[mati2,matj2]=find(ncc_map2==max(max(ncc_map2)),1,'first');
            %[ith(idxv),jth(idxv); mati1 , matj1 ; mati2 , matj2 ; tl(1),tl(2)]
            %input('')
            
        
      
         %save(filename,'res','noises_vars','points','filters_boxes','thresholds','max_iter');
       
        
    end
    %input('')
end
%save(['backup_',filename],'res','noises_vars','points','filters_boxes','thresholds','max_iter');
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

