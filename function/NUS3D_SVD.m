function [ X,Xdiff,RecTime ] = NUS3D_SVD( InArg  )

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% NUS3D_SVD to solve
%                     arg(X) min 1/2 * (||Mask.*(X-Y)||_2).^2 + lambda*||HX||_*
%                      where, H is the operator of block Hankel matrix pertutation
% Input :
%               Y : 3D Undersampled data
%            Mask : sampled mask
%          lambda : regularized parameter
% Output :
%               X : reconstructed data
%           Xdiff : the differentiation between neighbor iteration X
%              OV : Objective function value through iteration
% Authored by Enping Lin enpinglin@qq.com
%             Ze Fang    zefang23@qq.com

    Y = InArg.YM;
    Mask = InArg.mask;
    lambda = InArg.lambda;
    if ~isfield(InArg,'maxloop')
       maxloop = 100;
    else
       maxloop= InArg.maxloop;
    end 
    if ~isfield(InArg,'Wyes')
       Wyes = 1;
    else
       Wyes= InArg.Wyes;
    end    

    tic;
    [N1,N2,N3,N4] = size(Y);
    sr = ceil(size(Y,2)/2);
    sc = size(Y,2) - sr+1;
    slr = ceil(size(Y,3)/2);
    slc = size(Y,3) - slr + 1;
    nr = ceil(size(Y,1)/2);
    n2 = sr;
 %   Generate HhH weight   
    nc = size(Y,1) - nr + 1;
    t1 = min(nr,nc) - 1;
    t2 = max(nr,nc) - t1;
    rw = ([1:1:t1, (t1+1)*ones(1,t2),t1:-1:1]).';
    
    t3 = min(sr,sc)-1;
    t4 = max(sr,sc)-t3;
    cw = ([1:1:t3, (t3+1)*ones(1,t4),t3:-1:1]);

    t5 = min(slr,slc)-1;
    t6 = max(slr,slc)-t5;
    lw = ([1:1:t5, (t5+1)*ones(1,t6),t5:-1:1]);
    HhH = zeros(length(rw),length(cw),length(lw));
    for it_HhH = 1:length(lw)
        HhH(:,:,it_HhH) = rw*cw*lw(it_HhH);
    end
    HhH = repmat(HhH,[1,1,N4]);
    Xlast = Y;

    Mask = repmat(Mask,[1,1,N4]);
    if Wyes
        WMask = HhH .* Mask;
    else
        WMask = Mask;
    end
    mu = 0.1;

    % ---------------- Initialization --------------------
    Z = zeros(sr*nr*slr,sc*nc*slc,N4);
    for c_temp = 1:N4
        Z_temp = Matrix2BHankel_3D(Y(:,:,:,c_temp),nr,n2,slr); 
        Z(:,:,c_temp) = Z_temp;
    end
    Y = reshape(Y,[N1,N2,N3*N4]); 
    D = zeros(size(Z));
    HX = zeros(size(Z));
    [b1,b2,~] = size(Z);
    Xdown = 1*WMask+mu*HhH;

    HhmuZ_D = zeros(N1,N2,N3,N4);
    for c_temp = 1:N4
        muZ_D_part = mu*Z(:,:,c_temp)-D(:,:,c_temp);
        HhmuZ_D(:,:,:,c_temp) = BHankel2Matrix_3D( muZ_D_part,sr,sc,slr,slc );
    end
    HhmuZ_D = reshape(HhmuZ_D,[N1,N2,N3*N4]);
    Xdiff=zeros(maxloop-1,1);

    for itloop = 1:maxloop 
        fprintf('Iteration: %d\n',itloop)
        % ------------ Update X ------------------
        Xup = WMask.*Y + HhmuZ_D;
        X = Xup./Xdown;
        if itloop >= 2
           Xdiff(itloop-1) = norm(X-Xlast,'fro')/norm(X,'fro') ;
        end 
        Xlast = X;
        % ------------ Update Z ------------------        
        for c_temp = 1:N4
            HX_temp = Matrix2BHankel_3D(X(:,:,(c_temp-1)*N3+1:c_temp*N3),nr,n2,slr);
            HX(:,:,c_temp) = HX_temp;
        end

        for c_temp = 1:N4
            Z_temp = SVT(HX(:,:,c_temp)+D(:,:,c_temp)/mu ,lambda/mu );
            Z(:,:,c_temp) = Z_temp;
        end

        % ------------ Update D ------------------
        D = D + mu*(HX-Z);

        for c_temp = 1:N4
            muZ_D_part = mu*Z(:,(c_temp-1)*b2+1:c_temp*b2) - D(:,(c_temp-1)*b2+1:c_temp*b2);
            HhmuZ_D(:,:,(c_temp-1)*N3+1:c_temp*N3) = BHankel2Matrix_3D( muZ_D_part,sr,sc,slr,slc );
        end
        if itloop >= 2 && Xdiff(itloop-1)<1e-3
             if itloop<maxloop 
                 Xdiff(itloop:end) = [];
             end    
             break;        
        end 
    end
    X = reshape(X,[N1,N2,N3,N4]);
    RecTime = toc/60;
    fprintf('Finish Iteration: %d, Time lapse:%5.2f min \n',itloop,RecTime)
end
