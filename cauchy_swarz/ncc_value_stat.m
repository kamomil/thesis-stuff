function [ measures ] = ncc_value_stat(filename,points, noises_vars,rand_pats,rand_points)

measures=struct;
sizes=[16,32,48,64,96,128];
 for i=1:length(noises_vars)
     for k=1:length(sizes)
         measures(i,k).vals=zeros(1,0);
     end
 end

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
        
        pattern=pattern-sum(pattern(:))/(numel(pattern));
        
        for nidx=1:length(noises_vars)
            J=imread(['noise_var',num2str(noises_vars(nidx)),'/',points(rand_points(i)).im_name]);
            J=mat2gray(J);
            
            %         %%%%%%% MATLABS NCC %%%%%%%%%%%%%%%%%%
            wins_sums= box_corr(J,[1,n_p,1,m_p],1,n_p,m_p);
            wins_mu=wins_sums/(n_p*m_p);
            wins_norm2 = box_corr(J.^2,[1,n_p,1,m_p],1,n_p,m_p);
            norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
            cor=conv2(J,rot90(pattern,2),'valid');
            ncc_map=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
            %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            measures(nidx,sizes==sz(1)).vals  = [measures(nidx,sizes==sz(1)).vals  ncc_map(tl(1),tl(2))];
            
        end
        
        
        save(filename,'measures','noises_vars','rand_pats','rand_points');
        
    end
    %input('')
end
%x
end
