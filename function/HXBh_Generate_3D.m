function [HXBh] = HXBh_Generate_3D(sr_nr_slr,sc,nc,slc,X,Bh)
% Authored by Enping Lin  enpinglin@qq.com
% 

    HXBh = zeros(sr_nr_slr,size(Bh,2));

    r0 = nc + size(X,1) - 1;
    c0 = sc + size(X,2) - 1;
    l0 = slc + size(X,3) - 1;

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
    
    Ind1 = nc : size(X,1);
    Ind2 = sc : size(X,2);
    Ind3 = slc : size(X,3);

    for it = 1:size(Bh,2)
        K = flip(flip(flip(reshape(Bh(:,it),nc,sc,slc),1),2),3);
        Temp0 = ifftn( X0.* fftn(K,[r0,c0,l0]) );
        Temp1 =  Temp0(Ind1, Ind2, Ind3);
        HXBh(:,it) = Temp1(:);
    end
end
