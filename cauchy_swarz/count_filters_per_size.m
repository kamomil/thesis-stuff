

function [ counts ] = count_filters_per_size( points,fromsz)

 counts=zeros(size(fromsz));
for i=1:length(points)

    surfs=points(i).surfs;
    
    
    
    win_sz= round(surfs.Scale*[12 12]);
  
    for k=1:length(surfs)

        if(win_sz(k,1)<fromsz(1))
            counts(1)=counts(1)+1;
            continue;
        end
        
        l=1;
        while(l<=length(fromsz) && win_sz(k,1)>=fromsz(l))
            l=l+1;
        end
        counts(l-1)=counts(l-1)+1;
    end

end
end

