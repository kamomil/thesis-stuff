function [ inds ] = sample_patterns(points,size_range,max_per_img,num)

inds=zeros(num,2);

indsidx=1;

for i=1:length(points)
        disp(['-----------------------computing for ',num2str(i),'th image: ',points(i).im_name,'----------']); 
    from_this_img=0;
    surfs=points(i).surfs;
    I=mat2gray(rgb2gray(imread(points(i).im_name)));
    k=1;
    while(indsidx<=num && k<=length(surfs) && from_this_img<max_per_img)
    
        [filter,tl,br] = get_filter_from_surf(I,surfs(k), [12 12]);    
        [n_p,m_p]=size(filter);
        if(n_p>=size_range(1) && n_p<=size_range(2))
            inds(indsidx,:)=[i,k];
            indsidx=indsidx+1;
            from_this_img=from_this_img+1;
        end
        k=k+1;
    end

end

inds=inds(1:indsidx-1,:);
end

function [filter,tl,br] = get_filter_from_surf(I,s, std_sz)

win_sz=[ round(s.Scale*std_sz(1)) , round(s.Scale*std_sz(2)) ];
center=[s.Location(:,2) s.Location(:,1)];

tl=round(center)-floor(win_sz/2);
br=tl+win_sz-1;

filter=I(tl(1):br(1),tl(2):br(2) );
end



