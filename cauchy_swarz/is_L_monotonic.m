
function [L,arr] = is_L_monotonic(img,n_p,m_p,box_arr,w_arr,threshold , reconst_filter,debug)
%function [ arr,foundloci,foundlocj, m ] = l2_match_with_cauchy_schwarz_bound1(img,pattern, n_p,m_p,box_arr,w_arr,threshold , reconst_filter,debug)

box_arr=int32(box_arr);
% construct integral image + zeros pad (for boundary problems)
int_img = integral_image(img);

%[n_p,m_p]=size(pattern);
wins_norm2 = box_corr(img.^2,[1,n_p,1,m_p],1,n_p,m_p);
%wins_norm2(1:6,1:6)
wins_norm = sqrt(wins_norm2);

%reconst_filter=reconstruct_filter(box_arr,w_arr,n_p,m_p,0);
%figure, imshow(r)
r2=norm(reconst_filter,'fro')^2;

L=r2 +  wins_norm2-2*norm(reconst_filter,'fro')*wins_norm;%L is a matrix with the bound for wach window
%L=r2 +  wins_norm2+wins_norm;%L is a matrix with the bound for wach window

C=zeros(size(L));

threshold0 = threshold;
threshold = threshold*n_p*m_p;
%disp(['threshold ',num2str(threshold)]);

if debug
    D=r2+wins_norm2-2*conv2(img,rot90(reconst_filter,2),'valid');
    m=min(min(D));
    [ti,tj]=find(D==m);
end

t=0;

if debug
    %figure(998), imagesc(D/(n_p*m_p)), title('D/(n_p*m_p) - true distance map'); colorbar; pause(0.3);
end
arr=zeros(1,length(box_arr));


while t<length(box_arr);
    
    %disp(['iteration #',num2str(t),' out of ',num2str(length(box_arr))])
    
    if debug
        disp(['true minimum: ',num2str(m),' L at minimum: ',num2str(L(ti,tj)),' threshold ',num2str(threshold)]);
        figure(999), imagesc([D/(n_p*m_p) L/(n_p*m_p)], [0 threshold0]), title('L/(n_p*m_p) -  cauchy-schwarz lower bound to distance map'); colorbar; pause(0.1);
        %input('');
    end
    if(debug && ~isempty(find(L>D,1)))
        error('L>D')
    end
    ind=L<threshold; 
    
    [i,j]=find(ind);
     
    t=t+1;
    arr(t)=100*length(i)/numel(wins_norm);
    
    C=C+box_corr_given_indices_c(int_img,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));
    
    f=zeros(n_p,m_p);
    f(box_arr(t,1):box_arr(t,2),box_arr(t,3):box_arr(t,4))=w_arr(t);
    reconst_filter=reconst_filter-f;
    
    L1=r2 + wins_norm2(ind)-2*C(ind);
    L2=-2*norm(reconst_filter,'fro')*wins_norm(ind); 
    Lbefore=L;
    L(ind)=L1+L2;%VERY IMPORTANT TO UPDATE ONLY WHERE L<THRESHOLD
    ll=length(find(Lbefore>L));
    if(ll>0)
        disp(['L decreas, in ',num2str(ll),' indices']);
        input('')
    end
        
end

end

function wind = index(w,ind)
wind=zeros(size(w));
wind(ind)=w(ind);
end



