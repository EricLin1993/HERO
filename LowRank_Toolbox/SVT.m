function [ Y ] = SVT( X,tau )
% Authored by Enping Lin  enpinglin@qq.com
% 2020.6.29

   [U,S,V] = svd(X,'econ');
   sv = diag(S);
   thr = tau  ;   
   Y = U*diag(max(sv-thr,0))*V';

end

