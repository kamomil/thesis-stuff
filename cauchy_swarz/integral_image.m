function [ iI ] = integral_image(I)

I = cumsum(cumsum(I,2),1);
iI = [zeros(1,size(I,2)+1) ; [zeros(size(I,1),1) , I]];

end

