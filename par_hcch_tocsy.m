%% This demo is used to reconstruct the apcE data using the HERO, MF, and SVD algorithms.
%% Authored by Enping Lin enpinglin@qq.com
%%             Ze Fang    zefang23@qq.com
clear all,close all,clc
addpath('LowRank_Toolbox');
addpath data
addpath function
%% 
load d2o_hcch_tocsy.mat
load mask_tocsy.mat 

left = 70;
AB_st = 55;
length = 160;
FID_Ideal = data;

%% 
[N1,N2,num_slice] = size(data);
N1 = N1/2;
N2 = N2/2;
N3 = 4;
X_temp = zeros(N1,N2,N3);
X1_temp = zeros(N1,N2,num_slice);
X2_temp = zeros(N1,N2,num_slice);
X3_temp = zeros(N1,N2,num_slice);
X4_temp = zeros(N1,N2,num_slice);
T1 = tic;

InArg = cell(length,1);
parfor (iter = 1:length,8)
    fprintf('parfor Iteration: %d\n',iter)
    FID_HC = squeeze(FID_Ideal(:,:,iter+left));
    RR = FID_HC(1:2:end,1:2:end);
    RI = FID_HC(1:2:end,2:2:end);
    IR = FID_HC(2:2:end,1:2:end);
    II = FID_HC(2:2:end,2:2:end);
    Y = [RR RI IR II];
    Y = reshape(Y,N1,N2,N3);
    %%
    YN = Y + 0*(randn(N1,N2)+1i*randn(N1,N2));
    YM = YN.*repmat(mask,[1,1,N3]);
    YM_max = max(YM(:));
    YM = YM./YM_max;
    lambda = 0.01;
    Wyes = 1;
    InArg{iter}.YM = YM;
    InArg{iter}.mask = mask;
    InArg{iter}.lambda = lambda;
    InArg{iter}.Wyes = 1;
    InArg{iter}.st = AB_st;
    InArg{iter}.maxloop = 100; 
    % [ X_temp,~,~,~ ] = NUS2D_HERO( InArg{iter} );
    [ X_temp,~,~,~ ] = NUS2D_MF( InArg{iter} );
    % [ X_temp,~,~ ] = NUS2D_SVD( InArg{iter} );
    X_temp = X_temp.*YM_max;
    X1_temp(:,:,iter+left) = X_temp(:,:,1);
    X2_temp(:,:,iter+left) = X_temp(:,:,2);
    X3_temp(:,:,iter+left) = X_temp(:,:,3);
    X4_temp(:,:,iter+left) = X_temp(:,:,4);
end

fprintf('Consume %f mins \n',toc(T1)/60)

for iter2 = 1:271

    if iter2>=71 && iter2<=230
        continue;
    end

    FID_HC = squeeze(FID_Ideal(:,:,iter2));
    RR = FID_HC(1:2:end,1:2:end);
    RI = FID_HC(1:2:end,2:2:end);
    IR = FID_HC(2:2:end,1:2:end);
    II = FID_HC(2:2:end,2:2:end);

    Y = zeros(N1,N2,N3);
    Y(:,:,1) = RR;
    Y(:,:,2) = RI;
    Y(:,:,3) = IR;
    Y(:,:,4) = II;

    YM = Y.*repmat(mask,[1,1,N3]);

    X1_temp(:,:,iter2) = YM(:,:,1);
    X2_temp(:,:,iter2) = YM(:,:,2);
    X3_temp(:,:,iter2) = YM(:,:,3);
    X4_temp(:,:,iter2) = YM(:,:,4);    

end

FID_Rec = zeros(size(data));
FID_Rec(1:2:end,1:2:end,:) = X1_temp;
FID_Rec(1:2:end,2:2:end,:) = X2_temp;
FID_Rec(2:2:end,1:2:end,:) = X3_temp;
FID_Rec(2:2:end,2:2:end,:) = X4_temp;

save('tocsy_recon_MF.mat','FID_Rec')
