function [HhAB] = HhAB_Generate(Cr,Cc,Rr,Rc,A,B)
% Authored by Enping Lin enpinglin@qq.com
%             Ze Fang    zefang23@qq.com
 
    cn = Cr+Cc-1;
    rn = Rr+Rc-1;
    HhAB = zeros(rn,cn);

    t1 = min(Cr,Cc);
    t2 = max(Cr,Cc);
    rw = ([1:1:t1-1, (t1)*ones(1,t2-t1+1), t1-1:-1:1]).';

    for it1 = 1:cn
        ri = min(it1,Cr);
        Temp = 0;
        countm = 0;%
        while ri>0 
            ci = it1+1-ri; % the ci-th block anti-diagonal
            if ci>Cc
                break;
            end    
            if countm <= 3 % this is the maximum number of Hankel matrix used to estimate the column of HhAB (default 3)
                Temp = Temp + A((ri-1)*Rr+1:ri*Rr,:)*B(:,(ci-1)*Rc+1:ci*Rc);  
                countm = countm +1; 
            else            
                break;
            end 
            ri = ri-1;
        end
        Temp = rw(it1)/countm *Temp;
        HhAB(:,it1) = Hankel2vec( Temp );
    end


      
end

