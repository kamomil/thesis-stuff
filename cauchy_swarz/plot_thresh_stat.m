load('thresh_stat.mat');

sz=[16,32,48,64,96,128];
nm=[151,192,153,158,89,57];
for szidx=1:length(sz)
   
    colors={'b','r','g','c','k'};
    lm={'o','+','*','.','x'};
    figure(szidx)
    r=res(szidx,:);
    for noiseidx=1:length(noises_vars)
        for thidx=1:length(thresholds)
            to_plot=[];
            for iter=1:length(r)
                if nm(szidx)-r(iter).nm(noiseidx,thidx)<4
                    to_plot=[to_plot,r(iter).ind_av_per(noiseidx,thidx)/r(iter).nm(noiseidx,thidx)];
                else
                    break;
                end
            end
            if ~isempty(to_plot)
                hold on
                plot(to_plot,'Color',colors{noiseidx},'Marker',lm{thidx});
                hold off
            end
        end
    end
    title(num2str(sz(szidx)));
    text(1.5,1,'blue=0 noise, red=0.01 noise, green=0.05 noise cyan=0.1 noise , black=0.15 noise' )
    text(1.5,0.9,'o=0.99 thresh, +=0.9 , *=0.85 .=0.8 , x=0.75 noise' )
    
end