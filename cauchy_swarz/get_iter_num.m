function [ n ] = get_iter_num(alpha,alpha0,a,b,a0,psz)

b=b*psz;
M1=(b-a0-alpha0*a+alpha*a0);
M2=alpha0*(b-a-b*alpha);
M1/M2

n=(log(M1/M2)/log(alpha))+1;

end

