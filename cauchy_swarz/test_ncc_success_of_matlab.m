function [ bad_pats ] = test_ncc_success_of_matlab(points)


npoints=length(points);

%sd=[2 7 12 13 16 17 23 24 25 29 31 36 45];
next=0;
bad_pats=zeros(0,2);
for i=1:npoints
    
%     if sum(sd==i)>0
%         continue;
%     end
   
    I=rgb2gray(imread(points(i).im_name));
    I=mat2gray(I);

    %I=imnoise(I,'gaussian',0,0.05);
    
    for k=1:length(points(i).pats)
        
        next=next+1;         
        
        sz=points(i).pats(k).sz;
        disp([num2str(i),' ',num2str(k),' ',num2str(sz(1))]);
        tl=points(i).pats(k).top_left;
        br=tl+sz-1;
        
        pattern=I(tl(1):br(1),tl(2):br(2));
        
        [n_p,m_p]=size(pattern);          
%         boxes=filters_boxes1(i,k).boxes;
%         box_arr=boxes(:,1:4);
%         w_arr=boxes(:,5);
        pattern=pattern-sum(pattern(:))/(numel(pattern));
%         box_arr=box_arr(2:end,:);
%         w_arr=w_arr(2:end);
%          
%          rf1=reconstruct_filter(box_arr,w_arr,n_p,m_p,0);

%         if debug
%            % figure(995), imshow([pattern,rf1]), title('the pattern')
%         end
 
        max_iter=15;
        frac_for_direct=0.0001;
        
        
        
%         t01=tic;
        % [Lm,ith,jth,ind_leftf]=ncc_match_cauchy_schwarz_most_efficient(I,n_p,m_p,box_arr,w_arr,threshold,get_residual_pat_norms(pattern,box_arr,w_arr),max_iter,frac_for_direct);       
%         toc_t01=toc(t01);
        %disp([num2str(k),' ',num2str(toc_t01)]);
        
       % Lm(83:86,154:157)
       
       t01=tic;
        [U,ith2,jth2,vals,ind_num] = ncc_match_cauchy_with_mex(I,n_p,m_p,box_arr,w_arr,threshold,get_residual_pat_norms(pattern,box_arr,w_arr),max_iter,pattern,frac_for_direct);
        toc_t01=toc(t01);
%         
%         ith2=ith2(1:ind_num)+1;
%         jth2=jth2(1:ind_num)+1;
%         vals=vals(1:ind_num);
        %[ind_num,length(ith)]
        %[jth(end-100:end)' ; jth2(end-100:end)']
        
        %[jth(1:100)' ; jth2(1:100)']
%         figure(245)
%         plot(1:length(ind_leftf),ind_leftf/numel(Lm));
%         axis([1 length(ind_leftf) 0 1])
%         legend('ind left efficient');
%        
         t3=tic;
%         %%%%%%% MATLABS NCC %%%%%%%%%%%%%%%%%%
         wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
         wins_mu=wins_sums/(n_p*m_p);
         wins_norm2 = box_corr(I.^2,[1,n_p,1,m_p],1,n_p,m_p);
         norms_of_wins_minus_mu = abs(sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p))));
         
         norms_of_wins_minus_mu = norms_of_wins_minus_mu+100*eps;
         norms_of_wins_minus_mu(norms_of_wins_minus_mu<0.0001)=0.0001;
         cor=conv2(I,rot90(pattern,2),'valid');
         ncc_map=cor./(norms_of_wins_minus_mu*norm(pattern,'fro'));
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         toc_t3=toc(t3);
%      
        %min(min(norms_of_wins_minus_mu))   
        %min(min(ncc_map))
%         t4=tic;
%         ncc_map = normxcorr2(pattern, I);
%        
%         toc_t4=toc(t4);
%         ncc_map=ncc_map(sz(1):end-sz(1)+1,sz(2):end-sz(2)+1);
        
        [mati,matj]=find(ncc_map==max(max(ncc_map)),1,'first');
    
%             
        if(find(ncc_map==max(max(ncc_map))) ~=sub2ind(size(ncc_map), tl(1),tl(2)))
            '###error##'
            bad_pats=[ bad_pats ; [i,k] ];
%             [i,k,sz(1)]
%             [xu,yu,xn,yn,tl(1),tl(2)]
            %input('')
        end
            
        %input('')
      
    end
    %input('')
end
%x
end

function debug_direct_ncc(ith,jth,vals,ncc_map,norms_of_wins_minus_mu,pattern,cor)
    for i=1:length(ith)
        if (abs(vals(i)-ncc_map(ith(i),jth(i)))>0.005)
           disp(['id=',num2str(i),' i=',num2str(ith(i)),' j=',num2str(jth(i)),' v=',num2str(vals(i)), ' true=',num2str(ncc_map(ith(i),jth(i)))]);
           disp([' d=',num2str(norms_of_wins_minus_mu(ith(i),jth(i))*norm(pattern,'fro')),' corr=',num2str(cor(ith(i),jth(i)))]);
         
            input('')
            break;
        end
    end
    
end
function [r] = get_residual_pat_norms(pat,box_arr,w_arr)

r=zeros(1,length(w_arr)+1);
for i=1:length(w_arr)
    r(i)=norm(pat,'fro');
    
    f_r=box_arr(i,1);
    t_r=box_arr(i,2);
    f_c=box_arr(i,3);
    t_c=box_arr(i,4);
    
    pat(f_r:t_r,f_c:t_c)=pat(f_r:t_r,f_c:t_c)-w_arr(i);
end
r(end)=norm(pat,'fro');

end

function [b,w]=sort_boxes_by_size(boxes,w_arr)

s=zeros(length(w_arr),1);
for i=1:length(w_arr)
    s(i)=(boxes(i,2)-boxes(i,1)+1)*(boxes(i,4)-boxes(i,3)+1);
end

[s,ix]=sort(s,'descend');

b=boxes(ix,:);

w=w_arr(ix);
end


