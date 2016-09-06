function [ res ] = nnc_of_pats_with_noises_themself(points,noises_vars,rand_pats,rand_points)

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
        
        [n_p,m_p]=size(pattern);
        
        normelized_pattern=(pattern-sum(pattern(:))/(numel(pattern)));
        normelized_pattern=normelized_pattern/norm(normelized_pattern,'fro');
        
        
        for noiseidx=1:length(noises_vars)
            pnoised=imnoise(pattern,'gaussian',0,noises_vars(noiseidx));
            ncc=sum(sum(pnoised.*normelized_pattern));
            ncc=ncc/norm(pnoised,'fro');
            [noises_vars(noiseidx) sz(1) ncc]
            input('');
        end        
    end
    %input('')
end
%x
end

