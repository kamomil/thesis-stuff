
function [L,ind_left,i,j,L_for_rand_indexes,L_for_true_index,L_rand_first_term, L_rand_sec_term, L_true_first_term, L_true_sec_term] = l2_match_with_cauchy_schwarz_bound1(img,pat,n_p,m_p,box_arr,w_arr,threshold,true_index,max_iter)


box_arr=int32(box_arr);

int_img = integral_image(img);


wins_norm2 = box_corr(img.^2,[1,n_p,1,m_p],1,n_p,m_p);

wins_norm = sqrt(wins_norm2);

r2=norm(pat,'fro')^2;

L=r2 +  wins_norm2-2*norm(pat,'fro')*wins_norm;
residual_r2=r2;

%rand_mask=zeros(size(L));


Lsec=norm(pat,'fro')*wins_norm;


C=zeros(size(L));
C2=zeros(size(L));

threshold0 = threshold;
threshold = threshold*n_p*m_p;

t=0;

percent_ind_left=zeros(1,length(box_arr));
ind_left=zeros(1,min(max_iter,length(box_arr)));



%disp(['numbox=',num2str(length(box_arr)),',num of indices = ',num2str(numel(L))]);
ind=L<threshold;
r=find(ind);


rand_indexes=r(randperm(length(r),10));

rand_index=r(randperm(length(r),1));

L_for_rand_indexes=zeros(length(rand_indexes),min(max_iter,length(box_arr)));
L_for_true_index=zeros(1,min(max_iter,length(box_arr)));

L_rand_first_term=zeros(1,min(max_iter,length(box_arr)));
L_rand_sec_term=zeros(1,min(max_iter,length(box_arr)));
L_true_first_term=zeros(1,min(max_iter,length(box_arr)));
L_true_sec_term=zeros(1,min(max_iter,length(box_arr)));

while t<min(max_iter,size(box_arr,1));
    
    ind=L<threshold;
    
    [i,j]=find(ind);
    t=t+1;
    
    disp(['numbox=',num2str(length(box_arr)),' t=',num2str(t),' n=',num2str(length(i))]);
    ind_left(t)=length(i);
    L_for_rand_indexes(:,t)=L(rand_indexes)';
    L_for_true_index(t)=L(true_index(1),true_index(2));
    
    L_true_first_term(t)=C2(true_index(1),true_index(2));
    L_true_sec_term(t)=Lsec(true_index(1),true_index(2));
    
    L_rand_first_term(t)=C2(rand_index);
    L_rand_sec_term(t)=Lsec(rand_index);
    
    percent_ind_left(t)=100*length(i)/numel(wins_norm);
    
    C=C+box_corr_given_indices_c(int_img,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));
    
    C2=C2+box_corr_c(img,box_arr(t,:),w_arr(t),n_p,m_p);
    
    f_r=box_arr(t,1);
    t_r=box_arr(t,2);
    f_c=box_arr(t,3);
    t_c=box_arr(t,4);
    
    support_sz=(t_r-f_r+1)*(t_c-f_c+1);
    f=zeros(n_p,m_p);
    f(f_r:t_r,f_c:t_c)=w_arr(t);
    s=sum(sum(pat(f_r:t_r,f_c:t_c)));
    pat=pat-f;
    
    
    residual_r2=residual_r2 + (w_arr(t)*w_arr(t)*double(support_sz)) - (2*w_arr(t)*s) ;
    
    L1=r2 + wins_norm2(ind)-2*C(ind);
    
    Lsec=sqrt(residual_r2)*wins_norm;
    
    L2=-2*sqrt(residual_r2)*wins_norm(ind);
    
    L(ind)=L1+L2;
end
%figure(), plot(1:t,L_for_rand_indexes), title('L for different random indices');
%input('')

end


