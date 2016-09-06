%w_arr - the array of coefficients for the boxes
%box_arr - of size [k,4] where k is the number boxes, each box represented by
%startRow,endRow,startCol,endCol relative to the filter.
%[n_p , m_p] filter size
%top_lefts - Nx2 matrix of the top-left coordinates of the places we want to
%calculate the convolution
%sizes - Nx2 matrix of the sizes of the places we want to
%calculate the convolution

function [C] = box_corr_given_indices2(integral_img,box_arr,w_arr,n_p,m_p,ind_mat)


[n,m] = size(integral_img);
n=n-2;
m=m-2;
[i1,i2]=size(ind_mat);
if (i1 ~= n-n_p+1 || i2 ~= m-m_p+1)
    error('the indices matrix should be size(img)-size(pattern)+1 (result of valid conv)');
end

%create the ind_mat to be the size of the integral image
[tpi, tpj] = find(ind_mat);

C = zeros(n-n_p+1,m-m_p+1);



for pidx=1:length(tpi)
    
    x_start = tpi(pidx);
    y_start = tpj(pidx);
    
    
    arr_a = box_arr(:,1);
    arr_b = box_arr(:,2);
    arr_c = box_arr(:,3);
    arr_d = box_arr(:,4);
    
    % cumulate box responses
    k = size(box_arr,1); % == numel(w_arr)
    for i = 1:k
        a = arr_a(i);
        b = arr_b(i);
        c = arr_c(i);
        d = arr_d(i);
        
        
        C(tpi(pidx),tpj(pidx)) = C(tpi(pidx),tpj(pidx)) ...
            + w_arr(i) * ( integral_img(x_start+b,y_start+d) ...
            - integral_img(x_start+b,y_start+c-1) ...%bottom left
            - integral_img(x_start+a-1,y_start+d) ...
            + integral_img(x_start+a-1,y_start+c-1) );
        
        
    end  
end

end





