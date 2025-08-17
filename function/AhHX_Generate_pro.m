function [AhHX] = AhHX_Generate_pro(Cr,Cc,Rr,Rc,X,Ah)
% Authored by Enping Lin enpinglin@qq.com
%             Ze Fang    zefang23@qq.com
    
    BHcn = Cc*Rc;
    AhHX = zeros(size(Ah,1),BHcn);

    r0 = Rr + size(X,1) - 1;
    c0 = Cr + size(X,2) - 1;

    if mod(r0,2) ~= 0
        r0 = r0 + 1;
    end

    if mod(c0,2) ~= 0
        c0 = c0 + 1;
    end    

    X0 = fft2(X,r0,c0);
    Ind1 = Rr : size(X,1);
    Ind2 = Cr : size(X,2);

    for it = 1:size(Ah,1)
        K = flip( flip(reshape(Ah(it,:),Rr,Cr),1),2 );
        Temp0 = ifft2(X0.*fft2(K,r0,c0));
        Temp1 = Temp0(Ind1, Ind2);
        Temp1 = Temp1(:);
        AhHX(it,:) = Temp1.';
    end
end

