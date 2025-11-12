function [ TBH ] = Matrix2BHankel_3D(X,n1,n2,n3)
% Transform Matrix X into a Block Hankel matrix BH  
% Input :  
%         X : the matirx to be permuted in Block Hankel matrix
%         n1 : the row number of the sub Hankel matrix
%         n2 : the number of the sub Hankel matrix in the row of the Block Hankel matrix
% Authored by Ze Fang  zefang23@qq.com
% 2024.8.30
%%--------------------------------------------------------------
   [r,c,l] = size(X);
   sl = 1:l;

   hcl = r - n1 + 1; 
   shcl = c - n2 + 1;
   hi = l - n3 + 1;
   TBH = zeros(n1*n2*n3,hcl*shcl*hi);
   for it = 1:l
       Hc(:,:,it) = Matrix2BHankel(X(:,:,it),n1,n2);
   end  

   Hsl = hankel(sl(1:n3),sl(n3:end));
   for it1 = 1:size(Hsl,1)
      for it2 = 1:size(Hsl,2) 
         TBH( (it1-1)*n1*n2+1:it1*n1*n2,(it2-1)*hcl*shcl+1:it2*hcl*shcl ) = squeeze(Hc(:,:,Hsl(it1,it2)));
      end
   end  
end

