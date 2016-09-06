load('ncc_vals_stat.mat');

sizes=[16,32,48,64,96,128];

figure(123652)
c={'b','g','r','c','k','m'};
next=1;

means=zeros(length(noises_vars),length(sizes));
stds=zeros(length(noises_vars),length(sizes));
noises_vars=sqrt(noises_vars);
for szidx=1:length(sizes)
    to_plot=zeros(1,length(noises_vars));
    sigma=zeros(1,length(noises_vars));
    for nidx=1:length(noises_vars)
        to_plot(nidx)=mean(measures(nidx,szidx).vals);
        sigma(nidx)=std(measures(nidx,szidx).vals);
        means(nidx,szidx)=mean(measures(nidx,szidx).vals);
        stds(nidx,szidx)=std(measures(nidx,szidx).vals);
    end
    %mean(measures(nidx,szidx).vals)
    %std(measures(nidx,szidx).vals);
    
    errorbar(sqrt(noises_vars),to_plot,sigma,c{next});
    hold on
    next=next+1;
end
title('mean and standart deviation of ncc','FontSize',20);
xlabel('Standart deviation','FontSize',15);
ylabel('mean ncc value','FontSize',15);

legend('16','32','48','64','96','128');
means 
stds
means-stds
    