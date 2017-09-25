function [csm]=est_coilmaps_2coil(kloc,kdata,w_cal,n1,n2)
nt = size(kdata,4);
nc = size(kdata,3);

im_mch_gr = zeros(n1,n2,nt,nc);

for t = 1:nt
    GFFT_us = NUFFT( kloc(:,:,t), ones(size(kloc(:,:,t))), [0,0], [n1 n2]);clc;

    for c =1:nc
        im_mch_gr(:,:,t,c) = (GFFT_us'*(kdata(:,:,c,t).*((w_cal(:,:,t)))));
    end
end

im_mch_gr=squeeze(mean(im_mch_gr,3));
[csm,~] = ismrm_estimate_csm_walsh(im_mch_gr, 20,0.0);

figure;montage(reshape(abs(csm),[n1 n2 1 nc]),'DisplayRange', [0 1]); title('Magnitude of coil sensitivity maps');
figure;montage(reshape(angle(csm),[n1 n2 1 nc]),'DisplayRange', [-pi pi]); title('Phase of coil sensitivity maps');


