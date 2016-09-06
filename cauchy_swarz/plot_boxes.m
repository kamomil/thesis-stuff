function [ boxes_to_plot ] = plot_boxes( box_arr,w_arr,n_p,m_p )

boxes_to_plot=zeros(n_p,0,3);
redstrip=zeros(n_p,4,3);
redstrip(:,:,1)=1;  
for t=1:length(w_arr)
    b=zeros(n_p,m_p,3);
    f_r=box_arr(t,1);
    t_r=box_arr(t,2);
    f_c=box_arr(t,3);
    t_c=box_arr(t,4);
    b(f_r:t_r,f_c:t_c,:)=w_arr(t);
    b=mat2gray(b);
%     brgb=zeros(n_p,m_p,3);
%     brgb(:,:,1)=b;
%     brgb(:,:,2)=b;
%     brgb(:,:,3)=b;    
    boxes_to_plot=[boxes_to_plot ,redstrip, b];
    
end

end

