function [HXBh] = HXBh_Generate_pro(Cr,Cc,Rr,Rc,X,Bh)
% Authored by Enping Lin enpinglin@qq.com
%             Ze Fang    zefang23@qq.com
    
    BHrn = Cr*Rr;
    HXBh = zeros(BHrn,size(Bh,2));
    
    r0 = Rc + size(X,1) - 1;
    c0 = Cc + size(X,2) - 1;

    if mod(r0,2) ~= 0
        r0 = r0 + 1;
    end

    if mod(c0,2) ~= 0
        c0 = c0 + 1;
    end    
    
    X0 = fft2(X,r0,c0);
    Ind1 = Rc : size(X,1);
    Ind2 = Cc : size(X,2);

    for it = 1:size(Bh,2)
        K = flip( flip(reshape(Bh(:,it),Rc,Cc),1),2 );        
        Temp0 = ifft2(X0.*fft2(K,r0,c0));
        Temp1 = Temp0(Ind1, Ind2);
        HXBh(:,it) = Temp1(:);
    end 
end



  