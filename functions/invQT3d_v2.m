function Tinv=invQT3d_v2(opt)
% function to calculate the explicit inverse of: p1*R'R + p1*I + p3*I
%
% z = argmin p1*|| TV(z) - v1 - e11 ||^2 + p1*|| LLR(z) - v2 - e12 ||^2  + p3*|| x - z - e3 ||^2
% z = ( p1*R'R + p1*I + p3*I )^(-1) * [ p1*R'(v1+e11) + p1*(v2+e12) + p3*(x-e3) ];

pI=0;
if opt.lambda3~=0
    pI=pI+opt.p1;
end
pI = pI + +opt.p3;

kx=opt.size(1);
ky=opt.size(2);
nt=opt.size(3);

if (opt.lambda1==0)
    T1=0;
    T2=0;
else    
    T1z=zeros(kx,ky,nt);
    T2z=T1z;

    t1=[1;-1];
    T1z(1:2,1,1)=t1;
    T2z(1,1:2,1)=t1;

    T1=fftn(T1z);
    T2=fftn(T2z);
end

if (opt.lambda2==0)
    T3=0;
else    
    T3z=zeros(kx,ky,nt);

    t1=[1;-1];
    T3z(1,1,1:2)=t1;
    
    T3=fftn(T3z);
end

% make phase shift in Fourier domain to move circushift inside
[ph2,ph1,ph3]=meshgrid(0:ky-1,0:kx-1,0:nt-1);

shift1=exp(1j*2*pi*ph1/kx);
shift2=exp(1j*2*pi*ph2/ky);
shift3=exp(1j*2*pi*ph3/nt);

Tinv=opt.p1*(-T1.*shift1.*T1-T2.*shift2.*T2-T3.*shift3.*T3)+pI;

end
