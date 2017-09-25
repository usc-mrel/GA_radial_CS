clear all;clc;
% This script realizes reconstruction of 2D radial dynamic data, 
% using parallel imaging (SENSE), compressed sensing and locally low rank constrained recon
%
% Cost function:
% f(x) = || NUFFT( sMaps * x ) - kU ||^2 + lambda1* ||TV(x)||^1 
% + lamda3* ||LLR(x) ||^p
%
% Algorithm: ADMM with variable splitting

% When using the code, please cite 
% X. Miao et al. Magnetic Resonance Imaging 2016;34(6):707 - 714

% Author: Xin Miao, Yi Guo
% 07/12/2015

addpath('nufft_toolbox')
addpath('functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
data_file='test_data.mat';
overSampling = 2;

lambdas_LLR =[0]; % regularization weight for locally low rank constraint (LLR)
lambdas_TV = [0.002]; % regularization weight for total variation constraint (TV)
schatten_p =[0.5]; % Schatten-p norm used on the LLR constraint

usePar = 0; % parallel computing 
overlap = 0; % overlapping patches for LLR
do_shift = 0; % shifting patches for LLR
do_plot = 1;
max_iter = 50;

% Initalize options from parsed inputs
opt = struct('useParrallel', usePar, 'overlap', overlap, 'do_shift', ...
    do_shift, 'plot', do_plot, 'itermax', max_iter);

if opt.useParrallel
    %Clear any open interactive parrallel computing session
    try
        delete(gcp);
    catch
    end
    parpool(12);            %Setup matlab parrallel pool
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data

load(data_file);
k_rad = data; clear data;
k_rad = k_rad/max(abs(k_rad(:)));
[nFE,nP,nC]=size(k_rad);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate trajectory

[kloc] = Goldenratiosampling(nFE,1,nP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate sampling density

H = designFilter('ram-lak', nFE, 1);
bn = repmat(H,1,nP);
w = fftshift(bn,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporally segment k_rad

nspokes=13;
nt=floor(nP/nspokes);% number of frames

% crop the data according to the number of spokes per frame
k_rad=k_rad(:,1:nt*nspokes,:);
kloc=kloc(:,1:nt*nspokes);
w=w(:,1:nt*nspokes);

% sort the data into a time-series
for ii=1:nt
    k_radu(:,:,:,ii)=k_rad(:,(ii-1)*nspokes+1:ii*nspokes,:);
    klocu(:,:,ii)=kloc(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coil compression and smap computation
disp('Coil compression and smap computation');

[b1]=est_coilmaps_2coil(klocu,k_radu,wu,nFE/overSampling,nFE/overSampling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare NUFFT operator
% multicoil NUFFT operator
param.E=MCNUFFT(klocu,wu,b1);
% undersampled data
param.y=k_radu;
clear k_rad k_radu kloc klocu wu w
param.recon_nufft=param.E'*param.y;

%% %%%%%%%%%%%%%%%%%%   ADMM recon  %%%%%%%%%%%%%%%%%%%%%%%%
%% Recon parameters setting
opt.B = 5;					% Patch Size
opt.lambda1=0;			  % Spatial TV penalty
opt.lambda2=lambdas_TV;			% Temporal TV penalty
opt.lambda3=lambdas_LLR;				% LLR penalty
opt.xf=0;                   % set to 1 to turn on xf

opt.p1=0.05; 	% ADMM: dummy variable penalty for two splitting
opt.p3=0.05;
opt.p=schatten_p;
opt.class='double';
opt.update=1; % plot update number
opt.size = [size(b1,1),size(b1,2),nt];
param.max_iter = 10;

%% recon
tic,
[imgR,opt,cost1,cost2,cost3]=ADMM_recon_radial_LLR_TV(param,opt); 
imgR=double(abs(imgR));    

if opt.useParrallel
    %Clear any open interactive parrallel computing session
    delete(gcp);
end

% filename = sprintf('rest_LLR_%s_tTVlambda_%s_newNUFFT.mat',num2str(opt.lambda3),num2str(opt.lambda2));
% save(filename,'imgR','cost1','cost2','cost3');

