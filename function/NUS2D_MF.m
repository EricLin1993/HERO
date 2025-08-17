function [ X,Xdiff,OV,RecTime ] = NUS2D_MF(InArg )

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% NUS2D_MF to solve
%                     arg(X) min 1/2 * (||Mask.*(X-Y)||_2).^2 + lambda/2*(||A||F2+||B||F2)
%                         s.t. HX = A*B    
%                      where, H is the operator of block Hankel matrix
%                      pertutation, 1/2*(||A||F2+||B||F2) = ||HX||_*. 
% Input :
%               Y : 2D Undersampled data;
%            Mask : sampled mask;
%          lambda : regularized parameter
%   Output :
%               X : reconstructed data;
%           Xdiff : the differentiation between neighbor iteration X;
%              OV : Objective function value through iteration;
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
    if ~isfield(InArg,'st')
       InArg.st = 0.5*min(size(YM));
    else
       st= InArg.st;
    end   

    tic;
    [N1,N2,N3] = size(Y);
    sr = ceil(size(Y,2)/2);
    sc = size(Y,2) - sr+1;
    n1 = ceil(size(Y,1)/2);
    n2 = sr;
 %   Generate HhH weight   
    n3 = size(Y,1) - n1+1;
    t1 = min(n1,n3)-1;
    t2 = max(n1,n3)-t1;
    rw = ([1:1:t1, (t1+1)*ones(1,t2),t1:-1:1]).';
    
    t1 = min(sr,sc)-1;
    t2 = max(sr,sc)-t1;
    cw = ([1:1:t1, (t1+1)*ones(1,t2),t1:-1:1]);
    HhH = rw*cw;
    HhH = repmat(HhH,[1,N3]);
    Y = reshape(Y,[N1,N2*N3]); 
    Xlast = Y;

    Mask = repmat(Mask,[1,N3]);
    if Wyes
        WMask = HhH .* Mask;
    else
        WMask = Mask;
    end
    mu = 1;

    % ---------------- Initialization --------------------
    A=rand(sr*n1,st) - 0.5;
    B=rand(st,sc*n3*N3) - 0.5;

    D = zeros(sr*n1,sc*n3*N3);
    Xdown = 1*WMask+mu*HhH;

    HhmuAB_D = zeros(N1,N2*N3);
    for c_temp = 1:N3
        c1 = (c_temp-1)*sc*n3 + 1;
        c2 = c_temp*sc*n3;
        muAB_D_part = mu*A*B(:,c1:c2) - D(:,c1:c2);
        HhmuAB_D(:,(c_temp-1)*N2+1:c_temp*N2) = BHankel2Matrix( muAB_D_part,sr,sc );
    end

    Xdiff = zeros(maxloop-1,1);
    OV = zeros(maxloop,1);

    for itloop = 1:maxloop 
        fprintf('Iteration: %d\n',itloop)
        % ------------ Update X ------------------
        Xup = WMask.*Y + HhmuAB_D;
        X = Xup./Xdown;
        if itloop >= 2
           Xdiff(itloop-1) = norm(X-Xlast,'fro')/norm(X,'fro') ;
        end
        Xlast = X;
        
        % ------------ Update A,B ------------------
        Bh = B';
        BBh = zeros(st,st);
        for Bc = 1:N3
            Temp_BBh = B(:,(Bc-1)*sc*n3+1:Bc*sc*n3)*Bh((Bc-1)*sc*n3+1:Bc*sc*n3,:); 
            BBh = BBh + Temp_BBh;
        end   
 
        A_up = zeros(sr*n1,st);
        for c_temp = 1:N3
            c1 = (c_temp-1)*sc*n3 + 1;
            c2 = c_temp*sc*n3;
            HX(:,c1:c2) = Matrix2BHankel(X(:,(c_temp-1)*N2+1:c_temp*N2),n1,n2);
            X_part(:,c1:c2) = mu*HX(:,c1:c2)+D(:,c1:c2);
            A_up_part = X_part(:,c1:c2)*Bh((c_temp-1)*sc*n3+1:c_temp*sc*n3,:);
            A_up = A_up + A_up_part;
        end 

        A = A_up/(lambda*eye(size(BBh))+mu*BBh);

        % B new
        Ah = A';
        AhA = Ah*A;
        Bup = Ah*X_part;
        B = (lambda*eye(size(AhA))+mu*AhA) \ Bup;
        
        % ------------ Update D ------------------
        D = D+mu*(HX-A*B);  

        for c_temp = 1:N3
            c1 = (c_temp-1)*sc*n3 + 1;
            c2 = c_temp*sc*n3;
            muAB_D_part = mu*A*B(:,c1:c2) - D(:,c1:c2);
            HhmuAB_D(:,(c_temp-1)*N2+1:c_temp*N2) = BHankel2Matrix( muAB_D_part,sr,sc );
        end        

        % ------------ Estimate Objective Value --------------
        
        Nu = 1/2*( norm(B,'fro')+norm(A,'fro') );
        OV(itloop) = 0.5*( norm(Mask.*(X-Y),'fro') ).^2 + lambda*Nu;
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

