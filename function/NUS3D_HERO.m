function [ X,Xdiff,OV,RecTime ] = NUS3D_HERO( InArg  )

%---------------------- Block Hankel MF Memory saved --------------------------------
%----------------------------------------------------------------------
% NUS3D_MF_improved to solve
%                     arg(X) min 1/2 * (||Mask.*(X-Y)||_2).^2 + lambda/2*(||A||F2+||B||F2)
%                         s.t. HX = A*B    
%                      where, H is the operator of block Hankel matrix
%                      pertutation, 1/2*(||A||F2+||B||F2) = ||HX||_*. 
% Input :
%               Y : 3D Undersampled data
%            Mask : sampled mask
%          lambda : regularized parameter
% Output :
%               X : reconstructed data;
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
    [N1,N2,N3,N4] = size(Y);
    sr = ceil(size(Y,2)/2); % s:block, r:row, c:column, l:layer
    sc = size(Y,2) - sr + 1;
    slr = ceil(size(Y,3)/2);
    slc = size(Y,3) - slr + 1;
    nr = ceil(size(Y,1)/2); % n:sub
    nc = size(Y,1) - nr + 1;
    
    sr_nr_slr = sr * nr * slr;
    sc_nc_slc = sc * nc * slc;

 %   Generate HhH weight   
    t1 = min(nr,nc)-1;
    t2 = max(nr,nc)-t1;
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
    mu = 1;
    
    % ---------------- Initialization --------------------

%% generate A and B by random matrix
    A = rand(sr_nr_slr,st) - 0.5;
    B = rand(st,sc_nc_slc*N4) - 0.5;
    
%%
    Y = reshape(Y,[N1,N2,N3*N4]);
	d = zeros(size(Y));
    HhAB = zeros(N1, N2, N3, N4);

    for Bc = 1:N4
        B_part = B(:,(Bc-1)*sc_nc_slc+1:Bc*sc_nc_slc);
        HhAB(:,:,:,Bc) = HhAB_Generate_3Dpro(sr,sc,nr,nc,slr,slc,A,B_part);
    end

    HhAB = reshape(HhAB,[N1,N2,N3*N4]);
    Xdiff=zeros(maxloop-1,1);
    OV = zeros(maxloop,1);
    Xdown = 1*WMask+mu*HhH;
    Temp2 = zeros(st, sc_nc_slc*N4);

    for itloop = 1:maxloop 
        fprintf('Iteration: %d\n',itloop)
        % ------------ Update X ------------------
        Xup =  WMask.*Y + mu*HhAB - HhH.*d;
        X = Xup./Xdown;
        X = reshape(X,N1,N2,N3*N4); 

        if itloop >= 2
            Xdiff(itloop-1) = norm(X-Xlast,'fro')/norm(X,'fro') ;
        end 
        
        Xlast = X;

        % ------------ Update A,B ------------------
        % A new
        Bh = B';
        BBh = 0;
        Temp1 = zeros(sr*nr*slr,st);
        for Bc = 1:N4
            B_part = Bh((Bc-1)*sc_nc_slc+1:Bc*sc_nc_slc,:);
            Temp_BBh = B(:,(Bc-1)*sc_nc_slc+1:Bc*sc_nc_slc)*B_part; 
            BBh = BBh + Temp_BBh;

            X_part = mu*X(:,:,(Bc-1)*N3+1:Bc*N3)+d(:,:,(Bc-1)*N3+1:Bc*N3);
            Temp1_part = HXBh_Generate_3D(sr_nr_slr,sc,nc,slc,X_part,B_part);
            Temp1 = Temp1 + Temp1_part;
        end   

        A = Temp1/(lambda*eye(size(BBh))+mu*BBh);

        % B new
        Ah = A';
        AhA = Ah*A;

        for Ar = 1:N4
            X_part = mu*X(:,:,(Ar-1)*N3+1:Ar*N3)+d(:,:,(Ar-1)*N3+1:Ar*N3);
            Temp2_part = AhHX_Generate_3D(sc_nc_slc,sr,nr,slr,X_part,Ah);
            Temp2(:, (Ar-1)*sc_nc_slc+1:Ar*sc_nc_slc) = Temp2_part;
        end

        B = (lambda*eye(size(AhA))+mu*AhA) \ Temp2;
        
        % ------------ Update D ------------------
        HhAB = reshape(HhAB,[N1,N2,N3,N4]);       
        for Bc = 1:N4
            B_part = B(:,(Bc-1)*sc_nc_slc+1:Bc*sc_nc_slc);
            HhAB(:,:,:,Bc) = HhAB_Generate_3Dpro(sr,sc,nr,nc,slr,slc,A,B_part);
        end
        HhAB = reshape(HhAB,[N1,N2,N3*N4]);
        d = d + mu*(X - HhAB./HhH);

        % ------------ Estimate Objective Value --------------
        OV(itloop) = 0.5*( norm(Mask.*(X-Y),'fro') ).^2 + 0.5*lambda*(norm(A(:)).^2+norm(B(:)).^2);
        X = reshape(X,[N1,N2,N3,N4]);  
        if itloop >= 2 && Xdiff(itloop-1)<1e-3 % 1e-4
             if itloop<maxloop 
                 Xdiff(itloop:end) = [];
             end    
             break;        
        end 
    end
    RecTime = toc/60;
    fprintf('Finish Iteration: %d, Time lapse:%5.2f min \n',itloop,RecTime)
end
