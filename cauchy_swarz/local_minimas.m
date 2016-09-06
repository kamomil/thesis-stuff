function [ minimas ] = local_minimas( v )



minimas=zeros(1,length(v));
next=1;
for j=1:length(v)
    if(is_local_min(j,v))
        minimas(next)=j;
        next=next+1;
    end
end

minimas=minimas(1:(next-1));

end

function b=is_local_min(i,v)

if(i==1)
    if(v(i+1)>v(i))
        b=1;
    else
        b=0;
    end
    return;
end

if(i==length(v))
    if(v(i-1)>v(i))
        b=1;
    else
        b=0;
    end
    return;
end

if((v(i-1)<v(i) || v(i-1)==v(i)))
    b=0;
    return;
end

j=i+1;
while(j<length(v)+1 && v(j)==v(i))
   j=j+1; 
end

if(j==length(v)+1 || v(j)>v(i))
    b=1;
else
    b=0;
end

end