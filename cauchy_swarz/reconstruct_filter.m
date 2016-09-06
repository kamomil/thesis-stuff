function [ filter] = reconstruct_filter(box_arr,w_arr,n_p,m_p,show_process)

%figure(x)
filter=zeros(n_p,m_p);
for i=1:size(box_arr,1)
    f=zeros(n_p,m_p);
    f(box_arr(i,1):box_arr(i,2),box_arr(i,3):box_arr(i,4))=w_arr(i);
    filter = filter+f;
    if show_process && i<20
        imshow(mat2gray(filter));
        title(['box num=',num2str(i)])
        input('')
    end
end

end
