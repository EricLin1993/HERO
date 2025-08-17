function [ X,Xdiff,OV,RecTime ] = NUS2D_HERO( InArg  )

%---------------------- Block Hankel MF Memory saved --------------------------------
%----------------------------------------------------------------------
% NUS2D_HERO to solve
%                     arg(X) min 1/2 * (||Mask.*(X-Y)||_2).^2 + lambda/2*(||A||F2+||B||F2)
%                         s.t. HX = A*B    
%                      where, H is the operator of block Hankel matrix
%                      pertutation, 1/2*(||A||F2+||B||F2) = ||HX||_*. 
% Input :
%               Y : 2D Undersampled data
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
    Xlast = Y;

    Mask = repmat(Mask,[1,N3]);
    if Wyes
        WMask = HhH .* Mask;
    else
        WMask = Mask;
    end
    mu = 1;
    
    % ---------------- Initialization --------------------

%% generate A and B by random matrix
    A = rand(sr*n1,st) - 0.5;
    B = rand(st,sc*n3*N3) - 0.5;

%%
    Y = reshape(Y,[N1,N2*N3]);
	d = zeros(size(Y));

    for Bc = 1:N3
        B_part = B(:,(Bc-1)*sc*n3+1:Bc*sc*n3);
        HhAB(:,:,Bc) = HhAB_Generate(sr,sc,n1,n3,A,B_part);
    end
    HhAB = reshape(HhAB,[N1,N2*N3]);    
    Xdiff=zeros(maxloop-1,1);
    OV = zeros(maxloop,1);
    Xdown = 1*WMask+mu*HhH;
    for itloop = 1:maxloop 
        fprintf('Iteration: %d\n',itloop)
        % ------------ Update X ------------------
        Xup =  WMask.*Y + mu*HhAB - HhH.*d; % HhH.*d ;%    BHankel2Matrix( mu*A*B-D,sr,sc );
        X = Xup./Xdown;
        X = reshape(X,N1,N2*N3); 

        if itloop >= 2
            Xdiff(itloop-1) = norm(X-Xlast,'fro')/norm(X,'fro') ;
        end 
        
        Xlast = X;
        % ------------ Update A,B ------------------
        % A new
        Bh = B';
        BBh = zeros(st,st);
        for Bc = 1:N3
            Temp_BBh = B(:,(Bc-1)*sc*n3+1:Bc*sc*n3)*Bh((Bc-1)*sc*n3+1:Bc*sc*n3,:); 
            BBh = BBh + Temp_BBh;
        end   

        Temp1 = zeros(sr*n1,st);
        
        for Bc = 1:N3
            X_part = mu*X(:,(Bc-1)*N2+1:Bc*N2)+d(:,(Bc-1)*N2+1:Bc*N2);
            B_part = Bh((Bc-1)*sc*n3+1:Bc*sc*n3,:);
            Temp1_part = HXBh_Generate_pro(sr,sc,n1,n3,X_part,B_part);
            Temp1 = Temp1 + Temp1_part;                                                     
        end 

        A = Temp1/(lambda*eye(size(BBh))+mu*BBh);

        % B new
        Ah = A';
        AhA = Ah*A;
        Temp2 = [];
        
        for Ar = 1:N3
            X_part = mu*X(:,(Ar-1)*N2+1:Ar*N2)+d(:,(Ar-1)*N2+1:Ar*N2);
            Temp2_part = AhHX_Generate_pro(sr,sc,n1,n3,X_part,Ah);
            Temp2 = [Temp2 Temp2_part];
        end

        B = (lambda*eye(size(AhA))+mu*AhA) \ Temp2 ;
        
        % ------------ Update D ------------------
        HhAB = reshape(HhAB,[N1,N2,N3]);       
        for Bc = 1:N3
            B_part = B(:,(Bc-1)*sc*n3+1:Bc*sc*n3);
            HhAB(:,:,Bc) = HhAB_Generate(sr,sc,n1,n3,A,B_part);
        end
        HhAB = reshape(HhAB,[N1,N2*N3]);
        d = d + mu*(X - HhAB./HhH);

        % ------------ Estimate Objective Value --------------
        OV(itloop) = 0.5*( norm(Mask.*(X-Y),'fro') ).^2 + 0.5*lambda*(norm(A(:)).^2+norm(B(:)).^2);
        X = reshape(X,[N1,N2,N3]);  
        if itloop >= 2 && Xdiff(itloop-1)<1e-3
             if itloop<maxloop 
                 Xdiff(itloop:end) = [];
             end    
             break;        
        end 
    end

    RecTime = toc/60;
    fprintf('Finish Iteration: %d, Time lapse:%5.2f min \n',itloop,RecTime)
end
