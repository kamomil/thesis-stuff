function [ output_args ] = measure_distance_quality( points )

rand_points=randperm(53,10)
rand_pats=randperm(20,10)

%for i=1:length(rand_points)

ncc=zeros(10,2000);
sdd=zeros(10,2000);
next=1;
for i=1:length(rand_points)
    i
    I=rgb2gray(imread(points(rand_points(i)).im_name));
    
    
    I=mat2gray(I);
    %     if debug
    %         figure(996), imshow(I), title('the image')
    %     end
    
    
    k=rand_pats(i);
    
    %for k=1:length(points(i).pats)
    %for k=1:length(points(i).pats)
    
    %k=rand_pats(i);
    sz=points(rand_points(i)).pats(k).sz;
    tl=points(rand_points(i)).pats(k).top_left;
    br=tl+sz-1;
    
    pattern=I(tl(1):br(1),tl(2):br(2));
    
    
    wins_norm2 = box_corr(I.^2,[1,sz(1),1,sz(2)],1,sz(1),sz(2));
    
    ssd_map=sum(sum(pattern.^2))+wins_norm2-2*conv2(I,rot90(pattern,2),'valid');
    ssd_sorted=sort(ssd_map(:));
    
    ssd(next,:)=(ssd_sorted(1:2000)')/(sz(1)*sz(2));
    
    %figure(110)
    %imagesc(ssd_map);
    
    
    
    ncc_map = normxcorr2(pattern, I);
    ncc_map=ncc_map(sz(1):end-sz(1)+1,sz(2):end-sz(2)+1);
    ncc_sorted=sort(ncc_map(:),'descend');
    
    ncc(next,:)=ncc_sorted(1:2000)';
    
    %     figure(111)
    %     imagesc(ncc_map);
    
    %     figure(112)
    %     %ncc_sorted' ; ssd_sorted'])
    %     figure(1)
    %     plot(1:2000,ncc_sorted(1:2000)');
    %     figure(2)
    %     plot(1:2000,(ssd_sorted(1:2000)')/(sz(1)*sz(2)));
    
    % legend('ncc','ssd');
    % input('')
    next=next+1;
    
end

figure(1)
plot(1:2000,ncc), title('ncc best 2000');
figure(2)
plot(1:2000,ssd),title('ssd best 2000');

figure(3)
plot(1:100,ncc(:,1:100)), title('ncc best 100');
figure(4)
plot(1:100,ssd(:,1:100)),title('ssd best 100');


end


