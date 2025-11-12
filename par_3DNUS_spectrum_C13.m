%% This demo is used to reconstruct the IR24 data using the HERO, MF, and SVD algorithms.
%% Authored by Enping Lin enpinglin@qq.com
%%             Ze Fang    zefang23@qq.com
clear all,close all,clc
addpath('LowRank_Toolbox');
addpath data
addpath function
%% 
load C13_4D.mat
load mask_C13_4D.mat

%% parameter setting
left = 219;
AB_st = 240;
length_block = 240;
FID_Ideal = data;

%% 
[N1,N2,N3,num_block] = size(data);
N1 = N1/2;
N2 = N2/2;
N3 = N3/2;
N4 = 2^(length(size(data))-1);

X_temp = zeros(N1,N2,N3,N4);
X1_temp = zeros(N1,N2,N3,num_block);
X2_temp = zeros(N1,N2,N3,num_block);
X3_temp = zeros(N1,N2,N3,num_block);
X4_temp = zeros(N1,N2,N3,num_block);
X5_temp = zeros(N1,N2,N3,num_block);
X6_temp = zeros(N1,N2,N3,num_block);
X7_temp = zeros(N1,N2,N3,num_block);
X8_temp = zeros(N1,N2,N3,num_block);
F1 = squeeze(FID_Ideal(1,1,1,:));
F1 = abs(F1/max(F1(:)));
T1 = tic;
InArg = cell(length_block,1);
parfor (iter = 1:length_block,8)
    fprintf('parfor Iteration: %d\n',iter)

    FID_HC = squeeze(FID_Ideal(:,:,:,iter+left));
    RRR = FID_HC(1:2:end,1:2:end,1:2:end);
    RIR = FID_HC(1:2:end,2:2:end,1:2:end);
    IRR = FID_HC(2:2:end,1:2:end,1:2:end);
    IIR = FID_HC(2:2:end,2:2:end,1:2:end);
    RRI = FID_HC(1:2:end,1:2:end,2:2:end);
    RII = FID_HC(1:2:end,2:2:end,2:2:end);
    IRI = FID_HC(2:2:end,1:2:end,2:2:end);
    III = FID_HC(2:2:end,2:2:end,2:2:end);

    Y = zeros(N1,N2,N3,N4);
    Y(:,:,:,1) = RRR;
    Y(:,:,:,2) = RIR;
    Y(:,:,:,3) = IRR;
    Y(:,:,:,4) = IIR;
    Y(:,:,:,5) = RRI;
    Y(:,:,:,6) = RII;
    Y(:,:,:,7) = IRI;
    Y(:,:,:,8) = III;

    %%
    YM = Y.*repmat(mask,[1,1,1,N4]);
    YM_max = max(YM(:));
    YM = YM./YM_max;
    lambda = 0.01;  % SVD: 0.1, MF and HERO: 0.01
    Wyes = 1;

    InArg{iter}.YM = YM;
    InArg{iter}.mask = mask;
    InArg{iter}.lambda = lambda;
    InArg{iter}.Wyes = 1;
    InArg{iter}.st = AB_st;
    InArg{iter}.maxloop = 2; 

    % [ X_temp,~,~,~ ] = NUS3D_HERO( InArg{iter} );
    [ X_temp,~,~,~ ] = NUS3D_MF( InArg{iter} );
    % [ X_temp,~,~ ] = NUS3D_SVD( InArg{iter} );

    X_temp = X_temp.*YM_max;
    X1_temp(:,:,:,iter+left) = X_temp(:,:,:,1);
    X2_temp(:,:,:,iter+left) = X_temp(:,:,:,2);
    X3_temp(:,:,:,iter+left) = X_temp(:,:,:,3);
    X4_temp(:,:,:,iter+left) = X_temp(:,:,:,4);
    X5_temp(:,:,:,iter+left) = X_temp(:,:,:,5);
    X6_temp(:,:,:,iter+left) = X_temp(:,:,:,6);
    X7_temp(:,:,:,iter+left) = X_temp(:,:,:,7);
    X8_temp(:,:,:,iter+left) = X_temp(:,:,:,8);
end

fprintf('Consume %f mins \n',toc(T1)/60)

for iter2 = 1:512

    if iter2>=220 && iter2<=459
        continue;
    end


    FID_HC = squeeze(FID_Ideal(:,:,:,iter2));
    RRR = FID_HC(1:2:end,1:2:end,1:2:end);
    RIR = FID_HC(1:2:end,2:2:end,1:2:end);
    IRR = FID_HC(2:2:end,1:2:end,1:2:end);
    IIR = FID_HC(2:2:end,2:2:end,1:2:end);
    RRI = FID_HC(1:2:end,1:2:end,2:2:end);
    RII = FID_HC(1:2:end,2:2:end,2:2:end);
    IRI = FID_HC(2:2:end,1:2:end,2:2:end);
    III = FID_HC(2:2:end,2:2:end,2:2:end);

    Y = zeros(N1,N2,N3,N4);
    Y(:,:,:,1) = RRR;
    Y(:,:,:,2) = RIR;
    Y(:,:,:,3) = IRR;
    Y(:,:,:,4) = IIR;
    Y(:,:,:,5) = RRI;
    Y(:,:,:,6) = RII;
    Y(:,:,:,7) = IRI;
    Y(:,:,:,8) = III;
    YM = Y.*repmat(mask,[1,1,1,N4]);

    X1_temp(:,:,:,iter2) = YM(:,:,:,1);
    X2_temp(:,:,:,iter2) = YM(:,:,:,2);
    X3_temp(:,:,:,iter2) = YM(:,:,:,3);
    X4_temp(:,:,:,iter2) = YM(:,:,:,4); 
    X5_temp(:,:,:,iter2) = YM(:,:,:,5);
    X6_temp(:,:,:,iter2) = YM(:,:,:,6);
    X7_temp(:,:,:,iter2) = YM(:,:,:,7);
    X8_temp(:,:,:,iter2) = YM(:,:,:,8);      

end

FID_Rec = zeros(size(data));
FID_Rec(1:2:end,1:2:end,1:2:end,:) = X1_temp;
FID_Rec(1:2:end,2:2:end,1:2:end,:) = X2_temp;
FID_Rec(2:2:end,1:2:end,1:2:end,:) = X3_temp;
FID_Rec(2:2:end,2:2:end,1:2:end,:) = X4_temp;
FID_Rec(1:2:end,1:2:end,2:2:end,:) = X5_temp;
FID_Rec(1:2:end,2:2:end,2:2:end,:) = X6_temp;
FID_Rec(2:2:end,1:2:end,2:2:end,:) = X7_temp;
FID_Rec(2:2:end,2:2:end,2:2:end,:) = X8_temp;

FID_Rec = real(FID_Rec);
save('C13_4D_recon_MF.mat','FID_Rec')