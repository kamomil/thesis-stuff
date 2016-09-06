%w_arr the array of coefficients for the boxes
%box_arr of size [k,4] where k is the number boxes, each box represented by
%[start row, end row , start column , end column]
%
%integral_image - should be surrended by zeros
%
%it calculate the correlation on the given indices
%ind_mat - a logical matrix, the size (img)-size(pattern)+1 (result of
%valid conv). with 1 where the correlation should be calculated and zeros
%elsewhere
%
%
%
function [C] = box_corr_given_indices(integral_img,box_arr,w_arr,n_p,m_p, ind_mat)

[n,m] = size(integral_img);
n=n-2;
m=m-2;
[i1,i2]=size(ind_mat);
if (i1 ~= n-n_p+1 || i2 ~= m-m_p+1)
    error('the indices matrix should be size(img)-size(pattern)+1 (result of valid conv)');
end

%create the ind_mat to be the size of the integral image
tmp=zeros(n,m);
tmp(1:n-n_p+1,1:m-m_p+1)=ind_mat;
ind_mat=tmp;
ind_mat=[zeros(1,size(ind_mat,2)+2); [zeros(size(ind_mat,1),1) ind_mat zeros(size(ind_mat,1),1)]; zeros(1,size(ind_mat,2)+2)];
ind_mat=logical(ind_mat);

% initialize result matrix
C = zeros(n-n_p+1,m-m_p+1);

% cumulate box responses
k = size(box_arr,1);

%ind_mat
for i = 1:k
    a = box_arr(i,1);
    b = box_arr(i,2);
    c = box_arr(i,3);
    d = box_arr(i,4);
    
    z=zeros(n+2,m+2);
    z(ind_mat)=w_arr(i) * ( integral_img(circshift(ind_mat,[b-1,d-1])) ...
        - integral_img(circshift(ind_mat,[b-1,c-2])) ...
        - integral_img(circshift(ind_mat,[a-2,d-1])) ...
        + integral_img(circshift(ind_mat,[a-2,c-2])) );
    
    C = C +z(2:n-n_p+2,2:m-m_p+2);
end


