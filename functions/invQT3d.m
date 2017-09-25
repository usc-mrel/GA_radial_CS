function Tinv=invQT3d(opt)
% function to calculate the explicit inverse for TV

p1=0;
if opt.lambda3~=0
    p1=p1+opt.p1;
end

if opt.xf~=0
    p1=p1+opt.p1;
end

kx=opt.size(1);
ky=opt.size(2);
nt=opt.size(3);

T1z=zeros(kx,ky,nt);
T2z=T1z;
T3z=T2z;

t1=[1;-1];
% three direction TV
T1z(1:2,1,1)=t1;
T2z(1,1:2,1)=t1;
T3z(1,1,1:2)=t1;


T1=fftn(T1z);
T2=fftn(T2z);
T3=fftn(T3z);

% make phase shift in Fourier domain to move circushift inside
[ph2,ph1,ph3]=meshgrid(0:ky-1,0:kx-1,0:nt-1);

shift1=exp(1j*2*pi*ph1/kx);
shift2=exp(1j*2*pi*ph2/ky);
shift3=exp(1j*2*pi*ph3/nt);

if opt.lambda1~=0
Tinv=opt.p1*(-T1.*shift1.*T1-T2.*shift2.*T2-T3.*shift3.*T3)+p1+opt.p3;
else
Tinv=p1+opt.p3;
end

end
