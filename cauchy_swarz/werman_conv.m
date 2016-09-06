%f- the long signal
%p the pattern
%c = the (length d)-set of coefficients such that p(k)=sum(c.*(p(k+1:k+d))) (the homogeneous equation  
function y=werman_conv(f,p,c)

  n=length(f);
  r=length(p);
  d=length(c);

  p_with_tail=zeros(1,d+r);
  p_with_tail(1:r)=p;

  for i=1:d
    p_with_tail(r+i)=p_with_tail(r+i-d)-dot(p_with_tail(r+i-d+1:r+i-1),c(1:end-1));
  end
  ptail=p_with_tail(r+1:end);

  y=zeros(1,n-r+1);
  F=zeros(1,d);

  %---INTILIZATION--%
  for i=1:d
    s1=dot(p(1:d-i+1),f(i:d));
    y(i)=s1+sum(p(d-i+2:end).*f(d+1:i+r-1));
    s2=dot(ptail(1:d-i+1),f(i+r:d+r));
    F(i)=y(i)-s1+s2;
  end
  %---RUNNING--%
  for i=d+1:n-r
    y(i)=dot(rot90(c,2),F);
    for j=1:d-1
      F(j)=F(j+1)-p(d-j+1)*f(i)+ptail(d-j+1)*f(i+r);
    end
    F(d)=y(i)-p(1)*f(i)+ptail(1)*f(i+r);
  end
  y(end)=dot(rot90(c,2),F);

end