

function [U,i,j,vals] = linear_idx_ncc_match_cauchy_with_mex_clean(img,n_p,m_p,box_arr,w_arr,threshold,residual_pat_norms,max_iter,pat,frac_for_direct)
%ind_left=0;
iter_num=min([max_iter,length(w_arr)]);
box_arr=box_arr(1:iter_num,:);
w_arr=w_arr(1:iter_num);

%box_arr=int32(box_arr);
%%int_img = integral_image(img);
I=integral_image(img);
%I(1:10,1:10)
I2=integral_image(img.^2);
%I2(1:10,1:10)
%I2(end)
wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
wins_mu=wins_sums/(n_p*m_p);
wins_norm2 = box_corr(I2,[1,n_p,1,m_p],1,n_p,m_p);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));




% wins_sums= box_corr(int_I,[1,n_p,1,m_p],1,n_p,m_p);
%         wins_mu=wins_sums/(n_p*m_p);
%         wins_norm2 = box_corr(int_I2,[1,n_p,1,m_p],1,n_p,m_p);
%         norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
        
%norms_of_wins_minus_mu(1:5,1:5)
%'max norms'

%max(max(norms_of_wins_minus_mu))
%[n_p,m_p]
%norms_of_wins_minus_mu(end)
%figure(1)
%imshow(wins_mu);
%figure(2)
%m2g = mat2gray(norms_of_wins_minus_mu);
%m2g(1:5,1:5)
%imshow(m2g);

%size(norms_of_wins_minus_mu)
%wins_mu(1:5,1:5)
%wins_mu(end)
%size(norms_of_wins_minus_mu)

%'max I'
%max(max(I))

%'max I2'
%max(max(I2))
input('')
box_arr=box_arr';

%disp(['in linear_idx_ncc_match_with_mex_clean  frac_for_direct=',num2str(frac_for_direct),' sz',num2str(n_p)])
%size(I)
[U,lidx,vals,ind_num]=mex_cauchy_ncc_linear_idx_clean(single(img),int32(box_arr),single(w_arr),single(threshold),...%first 4
    single(pat),single(norms_of_wins_minus_mu),single(wins_mu),single(residual_pat_norms),single(frac_for_direct),single(I));%five

lidx=lidx(1:ind_num)+1;
[i,j] = ind2sub(size(U),lidx);
vals=vals(1:ind_num);
end


