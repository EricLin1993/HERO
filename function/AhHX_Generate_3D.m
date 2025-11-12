function [AhHX] = AhHX_Generate_3D(sc_nc_slc,sr,nr,slr,X,Ah)
% Authored by Enping Lin  enpinglin@qq.com
% 

    AhHX = zeros(size(Ah,1),sc_nc_slc);

    r0 = nr + size(X,1) - 1;
    c0 = sr + size(X,2) - 1;
    l0 = slr + size(X,3) - 1;

    if mod(r0,2) ~= 0
        r0 = r0 + 1;
    end

    if mod(c0,2) ~= 0
        c0 = c0 + 1;
    end   

    if mod(l0,2) ~= 0
        l0 = l0 + 1;
    end 
    
    X0 = fftn(X, [r0,c0,l0]);

    Ind1 = nr : size(X,1);
    Ind2 = sr : size(X,2);
    Ind3 = slr : size(X,3);
    
    for it = 1:size(Ah,1)
        K = flip(flip(flip(reshape(Ah(it,:),nr,sr,slr),1),2),3);
        Temp0 = ifftn( X0.* fftn(K,[r0,c0,l0] ) );
        Temp1 = Temp0(Ind1, Ind2 , Ind3);
        Temp1 = Temp1(:);
        AhHX(it,:) = Temp1.';
    end
end
