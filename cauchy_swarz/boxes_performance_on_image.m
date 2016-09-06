
%%%%%%%%%%%%%%%%%%%%%%%%
%Given an image and a pattern which may not be in the image, return an array that the i'th element is the percentage of indexes after the i'th box,   
%%%%%%%%%%%%%%%%%%%%%%%%
function [pref] = boxes_performance_on_image(img,filter,n_p,m_p,box_arr,w_arr,threshold , reconst_filter,debug, exact)
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

if exact
    L=r2 +  wins_norm2-2*norm(filter,'fro')*wins_norm;%L is a matrix with the bound for wach window
else
    L=r2 +  wins_norm2-2*norm(reconst_filter,'fro')*wins_norm;%L is a matrix with the bound for wach window
end



%L=r2 +  wins_norm2+wins_norm;%L is a matrix with the bound for wach window



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
arr=zeros(1,size(box_arr,1));
ind_left=zeros(1,size(box_arr,1));

ind=L<threshold; 
[i0,j0]=find(ind);
    
pref=zeros(1,length(w_arr));

reconst_filter_o=reconst_filter;

%disp(['numbox=',num2str(length(w_arr)),',num of indices = ',num2str(numel(L))]);
if(length(i0)==0)
    return
end
i=i0;
j=j0;
for t=1:length(w_arr);

    if exact
        L=r2 +  wins_norm2-2*norm(filter,'fro')*wins_norm;%L is a matrix with the bound for wach window
    else
        L=r2 +  wins_norm2-2*norm(reconst_filter_o,'fro')*wins_norm;%L is a matrix with the bound for wach window
    end


    if(w_arr(t)==0)
        continue;
    end
    f=zeros(n_p,m_p);
    f(box_arr(t,1):box_arr(t,2),box_arr(t,3):box_arr(t,4))=w_arr(t);
    reconst_filter=reconst_filter_o-f;
    

    C=box_corr_given_indices_c(int_img,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));
    L1=r2 + wins_norm2(ind)-2*C(ind);

    if exact
        L2=-2*norm(filter,'fro')*wins_norm(ind); 
    else
         L2=-2*norm(reconst_filter,'fro')*wins_norm(ind); 
    end 
    L(ind)=L1+L2;%VERY IMPORTANT TO UPDATE ONLY WHERE L<THRESHOLD
    
    ind=L<threshold;
    [i,j]=find(ind);
    %t
    %length(i)
    %if(length(i0)==0)
    %    error('BAAAAA')
    %end
    pref(t)=abs((length(i0)-length(i)))/length(i0);      
end
%plot(ind_left)
%plot(1:length(z),z)
end

function wind = index(w,ind)
wind=zeros(size(w));
wind(ind)=w(ind);
end



