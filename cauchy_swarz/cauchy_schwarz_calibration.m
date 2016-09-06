function [ ] = cauchy_schwarz_calibration(points,iminds,filters_boxes,debug,thresholds, noises)


%assuming there is no pattern with more than 1000 boxes
av_ind_num32=zeros(length(thresholds)*length(noises),1000);
av_ind_num96=zeros(length(thresholds)*length(noises),1000);

pat_num32=zeros(1,1000);
pat_num96=zeros(1,1000);


thresholds0=thresholds;

success_rate32=zeros(length(thresholds),length(noises));
success_rate96=zeros(length(thresholds),length(noises));

leg={};
for noiseidx=1:length(noises)
    disp(['NOISE IDX ',num2str(noiseidx)])
     thresholds=thresholds0+9*noises(noiseidx);
     for tidx=1:length(thresholds0)
                idx=tidx+(noiseidx-1)*length(thresholds);
                    leg{idx}=['noise = ',num2str(noises(noiseidx)),', threshold=',num2str(thresholds(tidx))];
     end
            
    for i=iminds%1:length(points)
        
        I=mat2gray(rgb2gray(imread(points(i).im_name)));
        
        J=imnoise(I,'gaussian',0,noises(noiseidx));
        if debug
            figure(996), imshow(I), title('the image')
        end
        
        for k=1:length(points(i).pats)
            
            sz=points(i).pats(k).sz(1);
            if(sz ~=32 && sz ~= 96)
                continue;
            end
            tl=points(i).pats(k).top_left;
            
            disp(['-----------------------computing for ',num2str(i),'th image: ',points(i).im_name,'----------']);
            disp(['-----------------------',num2str(k),'th filter ----------']);
            
            disp(' ');
            
            pat=I(tl(1):tl(1)+sz-1,tl(2):tl(2)+sz-1 );
            
            if debug
                figure(995), imshow(pat), title('the pattern')
            end
            
            boxes=filters_boxes(i,k).boxes;
            box_arr=boxes(:,1:4);
            w_arr=boxes(:,5);
            
            
           
            
            
            if(noiseidx==1)
                if sz==32
                    pat_num32=pat_num32+[ones(1,min(1000,length(w_arr))),zeros(1,1000-length(w_arr))];
                else
                    pat_num96=pat_num96+[ones(1,min(1000,length(w_arr))),zeros(1,1000-length(w_arr))];
                end
            end
                
            [sn,sm]=size(I);
            total_ind=(sn-sz+1)*(sm-sz+1);
            for tidx=1:length(thresholds)
               % 'THRESHOLD'
               % thresholds(tidx)
                [L,ind_left,ics,jcs]=l2_match_with_cauchy_schwarz_bound1(J,pat,sz,sz,box_arr,w_arr,thresholds(tidx));
                
                idx=tidx+(noiseidx-1)*length(thresholds);
                if sz==32
                    av_ind_num32(idx,1:length(ind_left))=av_ind_num32(idx,1:length(ind_left))+ind_left/total_ind;
                else
                    av_ind_num96(idx,1:length(ind_left))=av_ind_num96(idx,1:length(ind_left))+ind_left/total_ind;
                end
             
                Lia = ismember([ics,jcs],tl,'rows');
                if sum(Lia)
                    disp('YES') 
                    if sz==32
                        success_rate32(tidx,noiseidx)=success_rate32(tidx,noiseidx)+1;
                    else
                        success_rate96(tidx,noiseidx)=success_rate96(tidx,noiseidx)+1;
                    end
                else
                   disp( 'NO')
                end
                
                
            end
        end
        save('threshold_calibration.mat','av_ind_num32','av_ind_num96','success_rate32','success_rate96','pat_num32','pat_num96','thresholds0','noises','iminds');
    end
    
    
end
        
        av_ind_num32=av_ind_num32./repmat(pat_num32,length(thresholds)*length(noises),1);
        av_ind_num96=av_ind_num96./repmat(pat_num96,length(thresholds)*length(noises),1);
        
        success_rate32=success_rate32/pat_num32(1);
        success_rate96=success_rate96/pat_num96(1);
        
%         leg={};
%         for noiseidx=1:length(noises)
%             for tidx=1:length(thresholds0)
%                 idx=tidx+(noiseidx-1)*length(thresholds);
%                     leg{idx}=['noise = ',num2str(noises(noiseidx)),', threshold=',num2str(thresholds(tidx))];
%             end
%         end
        figure,
        plot(1:1000,av_ind_num32)
        legend(leg)
        figure,
        plot(1:1000,av_ind_num96)
        legend(leg)
        
        figure,
        plot(1:length(noises),success_rate32);
         figure,
        plot(1:length(noises),success_rate96);
        
        
        
        
        
    end
    
    
    
    
