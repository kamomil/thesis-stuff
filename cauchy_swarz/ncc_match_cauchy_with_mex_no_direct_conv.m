
function [U,i,j,vals,ind_frac,run_time] = ncc_match_cauchy_with_mex_no_direct_conv(img,n_p,m_p,box_arr,w_arr,threshold,residual_pat_norms,max_iter,pat,frac_for_direct)
%ind_left=0;
iter_num=min([max_iter,length(w_arr)]);
box_arr=box_arr(1:iter_num,:);
w_arr=w_arr(1:iter_num);

box_arr=int32(box_arr);
%%int_img = integral_image(img);
wins_sums= box_corr(img,[1,n_p,1,m_p],1,n_p,m_p);
wins_mu=wins_sums/(n_p*m_p);
wins_norm2 = box_corr(img.^2,[1,n_p,1,m_p],1,n_p,m_p);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
box_arr=box_arr';
 %img(end-1:end,end-1:end)
% size(img)

% box_arr(1:5)
% size(box_arr)
% threshold=0.99;
% w_arr(1:5)
% pat(1:3)
% frac_for_direct
% residual_pat_norms(1:3)
%[U,i,j,vals,ind_num]=mex_cauchy_ncc(img,box_arr,w_arr,threshold,pat,norms_of_wins_minus_mu,wins_mu,residual_pat_norms,frac_for_direct);
[U,i,j,vals,ind_frac,run_time,ind_num]=mex_cauchy_ncc_no_direct_conv(single(img),int32(box_arr),single(w_arr),single(threshold),single(pat),single(norms_of_wins_minus_mu),single(wins_mu),single(residual_pat_norms),single(frac_for_direct),int32(max_iter));
i=i(1:ind_num)+1;
j=j(1:ind_num)+1;
vals=vals(1:ind_num);
end


