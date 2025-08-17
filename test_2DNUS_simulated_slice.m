%% This demo is used to reconstruct synthetic data using the HERO, MF, and SVD algorithms.
%% Authored by Enping Lin enpinglin@qq.com
%%             Ze Fang    zefang23@qq.com

clear all,close all,clc
addpath('LowRank_Toolbox');
addpath data
addpath function

%% 
load simulated_data_100_100.mat
load mask_100_100.mat

% load simulated_data_150_150.mat
% load mask_150_150.mat

% load simulated_data_200_200.mat
% load mask_200_200.mat

% load simulated_data_250_250.mat
% load mask_250_250.mat

% load simulated_data_300_300.mat
% load mask_300_300.mat

[N1,N2] = size(RR);

FID_Ideal = zeros(2*N1, 2*N2);
FID_Ideal(1:2:end,1:2:end) = RR;
FID_Ideal(1:2:end,2:2:end) = RI;
FID_Ideal(2:2:end,1:2:end) = IR;
FID_Ideal(2:2:end,2:2:end) = II;

%% mask generation
N3 = 4;

Y = [RR RI IR II];
Y = reshape(Y,N1,N2,N3);

%%
YM = Y.*repmat(mask,[1,1,N3]);
YM_max = max(YM(:));
YM = YM./YM_max;
lambda = 0.01;
Wyes = 1;

InArg.YM = YM;
InArg.mask = mask;
InArg.lambda = lambda;
InArg.Wyes = 1;

% 100: HERO:lambda = 0.01, st = 120
%      MF:  lambda = 1   , st = 120
%      SVD: lambda = 0.1

% 150: HERO:lambda = 0.01, st = 180
%      MF:  lambda = 1   , st = 180
%      SVD: lambda = 0.1


% 200: HERO:lambda = 0.01, st = 240
%      MF:  lambda = 1   , st = 240
%      SVD: lambda = 0.1


% 250: HERO:lambda = 0.01, st = 325
%      MF:  lambda = 1   , st = 325
%      SVD: lambda = 0.1


% 300: HERO:lambda = 0.01, st = 405
%      MF:  lambda = 1   , st = 405
%      SVD: lambda = 0.1

InArg.st = 120; 
InArg.maxloop = 100; 
    
[ X_temp,Xdiff,OV,RecTime ] = NUS2D_HERO( InArg ); % HERO
% [ X_temp,Xdiff,OV,RecTime ] = NUS2D_MF( InArg ); % MF
% [ X_temp,Xdiff,RecTime ] = NUS2D_SVD( InArg ); % SVD

X_temp = X_temp.*YM_max;
X1_temp = X_temp(:,:,1);
X2_temp = X_temp(:,:,2);
X3_temp = X_temp(:,:,3);
X4_temp = X_temp(:,:,4);

FID_Rec = zeros(2*N1, 2*N2);
FID_Rec(1:2:end,1:2:end) = X1_temp;
FID_Rec(1:2:end,2:2:end) = X2_temp;
FID_Rec(2:2:end,1:2:end) = X3_temp;
FID_Rec(2:2:end,2:2:end) = X4_temp;

FID_Ide_temp = zeros(2*N1, N2);
FID_Rec_temp = zeros(2*N1, N2);

FID_Ide_temp(1:2:end,:) = RR + 1i * RI;
FID_Ide_temp(2:2:end,:) = IR + 1i * II;
spec_Ide_temp = fft(FID_Ide_temp, [], 2);
spec_Ide_temp = real(spec_Ide_temp);

spec_Ide = spec_Ide_temp(1:2:end, :) + 1i*spec_Ide_temp(2:2:end, :);
spec_Ide = fft(spec_Ide, [], 1);
spec_Ide = real(spec_Ide);
spec_Ide = spec_Ide/max(spec_Ide,[],'all');

FID_Rec_temp(1:2:end,:) = X1_temp + 1i * X2_temp;
FID_Rec_temp(2:2:end,:) = X3_temp + 1i * X4_temp;
spec_Rec_temp = fft(FID_Rec_temp, [], 2);
spec_Rec_temp = real(spec_Rec_temp);

spec_Rec = spec_Rec_temp(1:2:end, :) + 1i*spec_Rec_temp(2:2:end, :);
spec_Rec = fft(spec_Rec, [], 1);
spec_Rec = real(spec_Rec);
spec_Rec = spec_Rec/max(spec_Rec,[],'all');

rlne = norm((spec_Ide-spec_Rec),'fro')/norm(spec_Ide,'fro');
