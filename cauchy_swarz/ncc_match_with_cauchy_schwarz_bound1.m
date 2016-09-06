
function [U,ind_left,i,j,U_for_rand_indexes,U_given_indexes,U_for_true_index,U_rand_first_term, U_rand_sec_term, U_true_first_term, U_true_sec_term] = ncc_match_with_cauchy_schwarz_bound1(img,max_iter,n_p,m_p,box_arr,w_arr,threshold,residual_pat_norms,true_index,given_indexes)


box_arr=int32(box_arr);

int_img = integral_image(img);

wins_sums= box_corr(img,[1,n_p,1,m_p],1,n_p,m_p);

wins_mu=wins_sums/(n_p*m_p);

%pat_mu=sum(pat(:))/(n_p*m_p);


wins_norm2 = box_corr(img.^2,[1,n_p,1,m_p],1,n_p,m_p);

%std_dev_wins=sqrt(wins_norm2/(n_p*m_p) - wins_mu.*wins_mu);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));

% if ~isreal(norms_of_wins_minus_mu)
%     if max(max( imag(sqrt(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)))))>0.001
%         max(max( abs(imag(sqrt(wins_norm2-(wins_mu.*wins_mu*n_p*m_p))))))
%         error('norms of wins minu mu is not real')
%     else
%        norms_of_wins_minus_mu=real(norms_of_wins_minus_mu); 
%     end
% end

ss=norms_of_wins_minus_mu*residual_pat_norms(1);
U=ones(size(ss));
C=zeros(size(U));
C2=zeros(size(U));
W=zeros(size(U));


ind_left=zeros(1,min(max_iter,size(box_arr,1)));

ind=U>threshold;
r=find(ind);
rand_indexes=r(randperm(length(r),1));

U_for_rand_indexes=zeros(length(rand_indexes),min(max_iter,size(box_arr,1)));
U_given_indexes=zeros(length(given_indexes),min(max_iter,size(box_arr,1)));
U_for_true_index=zeros(1,min(max_iter,size(box_arr,1)));

U_true_first_term=zeros(1,min(max_iter,size(box_arr,1)));
U_true_sec_term=zeros(1,min(max_iter,size(box_arr,1)));

U_rand_first_term=zeros(1,min(max_iter,size(box_arr,1)));
U_rand_sec_term=zeros(1,min(max_iter,size(box_arr,1)));


t=0;
Uall=U;

U_true_first_term(1)=C(true_index(1),true_index(2));
U_true_sec_term(1)=W(true_index(1),true_index(2));
U_rand_first_term(:,1)=C(rand_indexes)';
U_rand_sec_term(:,1)=W(rand_indexes)';
    
while t<min(max_iter,size(box_arr,1));
    
    ind=U>threshold;    
    [i,j]=find(ind);
    t=t+1;
    
    %disp(['numbox=',num2str(length(box_arr)),' t=',num2str(t),' n=',num2str(length(i))]);
    ind_left(t)=length(i);
   
    U_for_rand_indexes(:,t)=Uall(rand_indexes)';
    U_given_indexes(:,t)=Uall(given_indexes)';
    
    %U_for_true_index(t)=U(true_index(1),true_index(2));
    U_for_true_index(t)=U(true_index(1),true_index(2));
    
    percent_ind_left(t)=100*length(i)/numel(norms_of_wins_minus_mu);
    
    
    C=C+box_corr_given_indices_c(int_img,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));
    
%     b=box_corr_c(img,box_arr(t,:),w_arr(t),n_p,m_p);
%     'b'
%     size(b)
%     size(C2)
    C2=C2+box_corr_c(img,box_arr(t,:),w_arr(t),n_p,m_p);
    
    cur_b=box_arr(t,:);
    support_sz=(cur_b(2)-cur_b(1)+1)*(cur_b(4)-cur_b(3)+1);
    W=W+double(support_sz)*w_arr(t)*wins_mu;
    
    U_true_first_term(t+1)=C2(true_index(1),true_index(2));
    U_true_sec_term(t+1)=W(true_index(1),true_index(2));
    
    U_rand_first_term(:,t+1)=C2(rand_indexes)';
    U_rand_sec_term(:,t+1)=W(rand_indexes)';
     
    Uall=(C2-W+residual_pat_norms(t+1)*norms_of_wins_minus_mu)./ss;
    
    L1=(C(ind)-W(ind)+residual_pat_norms(t+1)*norms_of_wins_minus_mu(ind))./ss(ind);
  
    U(ind)=L1;
end
%figure(), plot(1:t,L_for_rand_indexes), title('L for different random indices');
%input('')

end


