function [ X ] = BHankel2Matrix_3D( TBH,r,c,tsl,tn3)
% Transform the Block Hankel matrix into a  matrix    
%  Input: BH : the Block Hankel matrix
%         r : the number of the sub Hankel matrix in the row of the Block Hankel matrix
%         c : the number of the sub Hankel matrix in the column of the Block Hankel matrix
%  Output: X : the constructed matrix

% Authored by Enping Lin  enpinglin@qq.com
% 2021.3.29
%%--------------------------------------------------------------

   a = min(tsl,tn3);   % r=tsl, c=tn3
   b = max(tsl,tn3);
   [shs] = size(TBH)./[tsl,tn3];
   sr = shs(1); % n1
   sc = shs(2); % hcl
   
   for it1 = 1:a 
       s = zeros(sr,sc);
       for it2 =1:it1
          br = it1+1-it2; 
          bc = it2;
          s =s + TBH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
       end
       Hc(:,:,it1) = s;
   end  
 
  for it1 = a+1:b 
       s = zeros(sr,sc);
       if tsl <= tn3
           for it2 =it1-tsl+1:it1
              br =  it1+1-it2;
              bc = it2;
              s = s + TBH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
           end
       else
           for it2 =1:tn3
              br =  it1+1-it2;
              bc = it2;
              s = s + TBH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
           end
       end
       Hc(:,:,it1) = s;
%        Hc(:,:,it1) = s/a;
      
  end   
  for it1 = b+1:tn3+tsl-1 
       s = zeros(sr,sc);
       for it2 =it1-tsl+1:tn3
          br =  it1+1-it2;
          bc = it2;
          s = s + TBH((br-1)*sr+1:br*sr,(bc-1)*sc+1:bc*sc);
       end
       Hc(:,:,it1) = s;
%        Hc(:,:,it1) = s/(c+r-it1);

  end    
% ------------ Sub-Hankel to the collumn of Matirx --------------
   for it = 1:size(Hc,3)
       X(:,:,it) = BHankel2Matrix(Hc(:,:,it),r,c);
   end

end

