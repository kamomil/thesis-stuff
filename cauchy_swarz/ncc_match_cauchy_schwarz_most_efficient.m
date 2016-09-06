
function [U,i,j,ind_left] = ncc_match_cauchy_schwarz_most_efficient(img,n_p,m_p,box_arr,w_arr,threshold,residual_pat_norms,max_iter,frac_to_direct )
%ind_left=0;
iter_num=min([max_iter,length(w_arr)]);
box_arr=box_arr(1:iter_num,:);
w_arr=w_arr(1:iter_num);

box_arr=int32(box_arr);
%int_img = integral_image(img);
I=integral_image(img);
I2=integral_image(img.^2);
wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
wins_mu=wins_sums/(n_p*m_p);
wins_norm2 = box_corr(I2,[1,n_p,1,m_p],1,n_p,m_p);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));

C=box_corr_c(I,box_arr(1,:),w_arr(1),n_p,m_p);

f_r=box_arr(1,1);
t_r=box_arr(1,2);
f_c=box_arr(1,3);
t_c=box_arr(1,4);
support_sz=(t_r-f_r+1)*(t_c-f_c+1);

C=C-double(support_sz)*w_arr(1)*wins_mu;
t=2;

U=C./norms_of_wins_minus_mu;
ind_left=zeros(1,iter_num);
ind_left(1)=numel(U);
threshold*residual_pat_norms(1)-residual_pat_norms(t)
ind=U>(threshold*residual_pat_norms(1)-residual_pat_norms(t));
[i,j]=find(ind);
cur_frac = length(i)/numel(U);
%disp(['t=',num2str(t),' ind num=',num2str(length(i)),' cur_frac=',num2str(cur_frac)]);
while (t ~= iter_num+1 && cur_frac>frac_to_direct)

    ind_left(t)=length(i);
    C=C+box_corr_given_indices_c(I,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));
    f_r=box_arr(t,1);
    t_r=box_arr(t,2);
    f_c=box_arr(t,3);
    t_c=box_arr(t,4);
    support_sz=(t_r-f_r+1)*(t_c-f_c+1);
    
    C=C-double(support_sz)*w_arr(t)*wins_mu;
    t=t+1;
    threshold*residual_pat_norms(1)-residual_pat_norms(t)
    U(ind)=C(ind)./norms_of_wins_minus_mu(ind);%+residual_pat_norms(t)/residual_pat_norms(1);
    ind=U>(threshold*residual_pat_norms(1)-residual_pat_norms(t));
    [i,j]=find(ind);
    
    cur_frac = length(i)/numel(U);
 %   disp(['t=',num2str(t),' ind num=',num2str(length(i)),' cur_frac=',num2str(cur_frac)]);
end
