function [imgR,opt,cost1,cost2,cost3]=ADMM_recon_radial_LLR_TV(param,opt)
% function to realize reconstruction of 2D radial dynamic data, 
% using parallel imaging (SENSE), compressed sensing and locally low rank constrained recon
%
% input: kU: undersampled kspace, only the sampled points
%        sMaps: sensitivity Map, nt dimension should be 1
%        opt: recon parameters
% output: imgR: reconstructed image(nx*ny*kz)
%
% Cost function:
% f(x) = || NUFFT( sMaps * x ) - kU ||^2 + lambda1* ||TV(x)||^1 
% + lamda3* ||LLR(x) ||^p
%
% Algorithm: ADMM with variable splitting
% Cost function f(x,v1,v2) = 
% || NUFFT( sMaps * x ) - kU ||^2 + lambda1*|| v1 ||^1 + lamda3*|| v2 ||^p
% + p1*|| TV(z) - v1 - e11 ||^2 + p1*|| LLR(z) - v2 - e12 ||^2  + p3*|| x - z - e3 ||^2
%
% Iteratively solve z, v1, v2 and x:
%       z = argmin p1*|| TV(z) - v1 - e11 ||^2 + p1*|| LLR(z) - v2 - e12 ||^2  + p3*|| x - z - e3 ||^2
%       z = ( p1*R'R + p1*I + p3*I )^(-1) * [ p1*R'(v1+e11) + p1*(v2+e12) + p3*(x-e3) ];
%       v1 = argmin lambda1*|| v1 ||^1 + p1*|| TV(z) - v1 - e11 ||^2
%       v2 = argmin lamda3*|| v2 ||^p + p1*|| LLR(z) - v2 - e12 ||^2
%       x = argmin || NUFFT( sMaps * x ) - kU ||^2 + p3*|| x - z - e3 ||^2
%
%
% When using the code, please cite 
% X. Miao et al. Magnetic Resonance Imaging 2016;34(6):707 - 714

% Author: Xin Miao, Yi Guo
% 07/12/2015

nx=opt.size(1);ny=opt.size(2); nt=opt.size(3);

%% set initial values
x_k=param.recon_nufft;
[v1_k,LL]=compTx3d(x_k,opt);
v2_k=x_k;

e11_k=zeros(LL,1,opt.class); %residue variable for TV
e12_k=zeros(nx,ny,nt,opt.class); % residue variable for LLR
e3_k=zeros(nx,ny,nt,opt.class); % size of image

%% set different weighting for spatial and temporal TV
if ((opt.lambda1==0) &&  (opt.lambda2~=0))
    TVlambda=[ones([nx*ny*nt,1],opt.class)*opt.lambda2];
    disp('temporal TV');
else
    TVlambda=[ones([nx*ny*nt*2,1],opt.class)*opt.lambda1;ones([nx*ny*nt,1],opt.class)*opt.lambda2];
end

Tinv=invQT3d_v2(opt);

%% Begin iteration
iter=1;

while iter<opt.itermax

%% update z(k)

zftemp=opt.p3*(x_k-e3_k);
if opt.lambda1~=0||opt.lambda2~=0
    zftemp=zftemp+opt.p1*compThx3d(v1_k+e11_k,opt);
end

if opt.lambda3~=0
    zftemp=zftemp+opt.p1*(v2_k+e12_k);
end

z_k=ifftn(fftn(zftemp)./Tinv); 

%% update v1(k) using shrink
if opt.lambda1~=0||opt.lambda2~=0 
    v1_k=shrink1(compTx3d(z_k,opt)-e11_k,TVlambda/opt.p1); 
end

%% update v2(k)
if opt.lambda3~=0 % shrink LLR coeffient with svt
    if opt.lambda3<1000000
        v2_k=LLRp_v3(z_k-e12_k, opt.B, opt.lambda3/opt.p1, opt.p, opt.overlap, opt.do_shift, opt.useParrallel);     
    else
        v2_k = GLR_xm(z_k-e12_k, opt.lambda3/opt.p1, opt.p);
        disp('ktslr');
    end
end

%% update x(k): hard to solve explicitly, therefore use CG
tic;
x_k = CG_solver(x_k, param, opt.p3, (z_k+e3_k));
toc;

%% update residue variables: e11, e12 and e3
if opt.lambda1~=0||opt.lambda2~=0
    e11_k=e11_k-compTx3d(z_k,opt)+v1_k;
end

if opt.lambda3~=0
    e12_k=e12_k-z_k+v2_k;
end

e3_k=e3_k-x_k+z_k;

%% plot inter-iteration results
k_res= abs(param.E*x_k - param.y).^2;
k_res=sum(k_res(:));

if opt.lambda1~=0||opt.lambda2~=0
    obj11=compTx3d(x_k,opt);
else
    obj11=0;
end

if (opt.lambda3>0)&&(opt.lambda3<100000)
    obj12 = norm_evaluate(x_k,opt.p,opt.lambda3,opt.B,0);
elseif opt.lambda3>100000
    obj12 = norm_evaluate_glr(x_k,opt.p,opt.lambda3);
    disp('kt-SLR');
else
    obj12 = 0;
end

cost1(iter)=k_res;
cost2(iter)=sum(TVlambda.*abs(obj11));
cost3(iter)=abs(obj12);

cost(iter)=cost1(iter)+cost2(iter)+cost3(iter);

if iter>1
    if abs((cost(iter)-cost(iter-1))/cost(iter))<1e-6
        break;
    end
end

if mod(iter-1,opt.update)==0
    if opt.plot
%         tmp = abs(x_k(118:188,141:211,:));
        subplot(1,5,1);imshow(abs(x_k((-80:79)+nx/2+1,(-80:79)+nx/2+1,1)),[0 0.7*max(abs(x_k(:)))]);colormap('gray'); title(num2str(iter));axis equal; title(num2str(opt.lambda2));
%         subplot(1,5,2);imshow(repmat(squeeze(tmp(28:65,21,:)),[1 3]),[]);colormap('gray'); title(num2str(iter));axis equal; 
        subplot(1,5,3);plot(double(log10(cost1(1:iter)))); 
        subplot(1,5,4);plot(double(cost2(1:iter))); 
        subplot(1,5,5);plot(double(cost3(1:iter)));drawnow;

    end
end

iter=iter+1;

end

imgR=reshape(x_k,[nx,ny,nt]);
clear x_k;
