

function [U,i,j,vals,ind_frac,run_time,dconv_run_time] = linear_idx_ncc_match_cauchy_with_mex(img,n_p,m_p,box_arr,w_arr,threshold,residual_pat_norms,max_iter,pat,frac_for_direct)
%ind_left=0;
iter_num=min([max_iter,length(w_arr)]);
box_arr=box_arr(1:iter_num,:);
w_arr=w_arr(1:iter_num);

box_arr=int32(box_arr);
%%int_img = integral_image(img);
I=integral_image(img);
I2=integral_image(img.^2);
wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
wins_mu=wins_sums/(n_p*m_p);
wins_norm2 = box_corr(I2,[1,n_p,1,m_p],1,n_p,m_p);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
box_arr=box_arr';

[U,lidx,vals,ind_frac,run_time,dconv_run_time,ind_num]=mex_cauchy_ncc_linear_idx(single(img),int32(box_arr),single(w_arr),single(threshold),single(pat),single(norms_of_wins_minus_mu),single(wins_mu),single(residual_pat_norms),single(frac_for_direct),int32(max_iter),I);

lidx=lidx(1:ind_num)+1;
[i,j] = ind2sub(size(U),lidx);
vals=vals(1:ind_num);
end


