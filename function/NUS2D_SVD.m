function [ X,Xdiff,RecTime ] = NUS2D_SVD( InArg  )

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% NUS2D_SVD to solve
%                     arg(X) min 1/2 * (||Mask.*(X-Y)||_2).^2 + lambda*||HX||_*
%                      where, H is the operator of block Hankel matrix pertutation
% Input :
%               Y : 2D Undersampled data
%            Mask : sampled mask
%          lambda : regularized parameter
% Output :
%               X : reconstructed data
%           Xdiff : the differentiation between neighbor iteration X
%              OV : Objective function value through iteration
% Authored by Enping Lin enpinglin@qq.com

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
    [N1,N2,N3] = size(Y);
    sr = ceil(size(Y,2)/2);
    sc = size(Y,2) - sr+1;
    n1 = ceil(size(Y,1)/2);
    n2 = sr;
 %   Generate HhH weight   
    n3 = size(Y,1) - n1 + 1;
    t1 = min(n1,n3) - 1;
    t2 = max(n1,n3) - t1;
    rw = ([1:1:t1, (t1+1)*ones(1,t2),t1:-1:1]).';
    
    t1 = min(sr,sc)-1;
    t2 = max(sr,sc)-t1;
    cw = ([1:1:t1, (t1+1)*ones(1,t2),t1:-1:1]);
    HhH = rw*cw;
    HhH = repmat(HhH,[1,N3]);
    Xlast = Y;

    Mask = repmat(Mask,[1,N3]);
    if Wyes
        WMask = HhH .* Mask;
    else
        WMask = Mask;
    end
    mu = 0.1;

    % ---------------- Initialization --------------------

    Z = zeros(sr*n1,sc*n3,N3);
    for c_temp = 1:N3
        Z_temp = Matrix2BHankel(Y(:,:,c_temp),n1,n2); 
        Z(:,:,c_temp) = Z_temp;
    end
    Y = reshape(Y,[N1,N2*N3]); 
    D = zeros(size(Z));
    HX = zeros(size(Z));
    [b1,b2,~] = size(Z);
    Xdown = 1*WMask+mu*HhH;

    HhmuZ_D = zeros(N1,N2,N3);
    for c_temp = 1:N3
        muZ_D_part = mu*Z(:,:,c_temp)-D(:,:,c_temp);
        HhmuZ_D(:,:,c_temp) = BHankel2Matrix( muZ_D_part,sr,sc );
    end
    HhmuZ_D = reshape(HhmuZ_D,[N1,N2*N3]);
    
    Xdiff = zeros(maxloop-1,1);

    for itloop = 1:maxloop 
        fprintf('Iteration: %d\n',itloop)
        % tic;
        % ------------ Update X ------------------
        Xup = WMask.*Y + HhmuZ_D;
        X = Xup./Xdown;
        if itloop >= 2
           Xdiff(itloop-1) = norm(X-Xlast,'fro')/norm(X,'fro') ;
        end 
        Xlast = X;
        % ------------ Update Z ------------------        
        for c_temp = 1:N3
            HX_temp = Matrix2BHankel(X(:,(c_temp-1)*N2+1:c_temp*N2),n1,n2);
            HX(:,:,c_temp) = HX_temp;
        end

        for c_temp = 1:N3
            Z_temp = SVT(HX(:,:,c_temp)+D(:,:,c_temp)/mu ,lambda/mu );
            Z(:,:,c_temp) = Z_temp;
        end

        % ------------ Update D ------------------
        D = D + mu*(HX-Z);

        for c_temp = 1:N3
            muZ_D_part = mu*Z(:,(c_temp-1)*b2+1:c_temp*b2) - D(:,(c_temp-1)*b2+1:c_temp*b2);
            HhmuZ_D(:,(c_temp-1)*N2+1:c_temp*N2) = BHankel2Matrix( muZ_D_part,sr,sc );
        end

        % ------------ Estimate Objective Value --------------
        if itloop >= 2 && Xdiff(itloop-1)<1e-3
             if itloop<maxloop 
                 Xdiff(itloop:end) = [];
             end    
             break;        
        end        
    end
    X = reshape(X,[N1,N2,N3]);
    RecTime = toc/60;
    fprintf('Finish Iteration: %d, Time lapse:%5.2f min \n',itloop,RecTime)
end
