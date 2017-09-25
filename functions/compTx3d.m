function [Tx,LL]=compTx3d(x,opt)
% Calculate 2d spatial Total variation for 3d data set
x=reshape(x,[opt.size(1),opt.size(2),opt.size(3)]);

if ((opt.lambda1==0) &&  (opt.lambda2~=0))
    
TV3=x-circshift(x,[0 0 1]);
Tx=[TV3(:)];
LL=numel(x);
    
else

TV1=x-circshift(x,[1 0 0]);
TV2=x-circshift(x,[0 1 0]);
TV3=x-circshift(x,[0 0 1]);
Tx=[TV1(:);TV2(:);TV3(:)];
LL=numel(x)*3;

end

end