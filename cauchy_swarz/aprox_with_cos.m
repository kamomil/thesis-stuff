function [m,fs]=x(sig)

r=0:length(sig)-1;
d=dic(length(sig));

% for i=1:length(d)
%     co=cos(r*2*pi*d(i));
%     norm(co)
% end
% return



fs=[];

next=1;
while(norm(sig)>0.01)
    norm(sig)
    plot(sig)
    a=input('');
    f=d(1);
    co=cos(r*2*pi*d(1));
    m=norm(sig-(norm(sig)/norm(co))*co);
    a=max(abs(sig));
    for i=2:length(d)
        co=cos(r*2*pi*d(i));
        co=co/norm(co);
        c=norm(sig-norm(sig)*co);
        if c<m
            m=c;
            f=d(i);
        end
    end
    m
    fs(next)=f;
    next=next+1;
    sig=sig-(norm(sig)/norm(cos(r*2*pi*f)))*cos(r*2*pi*f);
end

%plot(r,[sig ; cos(r*2*pi*f) ; sig-cos(r*2*pi*f) ])
%legend('sig','app','diff');
    

% 
% 
% l=zeros(length(sig),1);
% 
% for i=1:length(l)
% f=fft(sig);
% f=dct(sig);
% l(i)=f(i);
% sig=sig-f(i)/length(sig);
% end


end

function freq=dic(n)

freq=[];
freq(1)=0;
next=2;
for i=1:n-1
    j=1;
    while(j/i<2)
        w=j/i;
        if(j/i>1)
            w=w-2;
        end
        freq(next)=w;
        next=next+1;
        j=j+1;
    end
    
end
freq=unique(freq);
end