function [HhAB] = HhAB_Generate_3Dpro(sr,sc,nr,nc,slr,slc,A,B)
% Authored by Ze Fang

    N1 = nr+nc-1;    
    N2 = sr+sc-1;
    N3 = slr+slc-1;
    HhAB = zeros(N1,N2,N3);

    t1 = min(nr, nc);
    t2 = max(nr, nc);
    st1 = min(sr, sc);
    st2 = max(sr, sc);

    w = ([1:1:t1-1, (t1)*ones(1,t2-t1+1), t1-1:-1:1]).';   
    sw = ([1:1:st1-1, (st1)*ones(1,st2-st1+1), st1-1:-1:1]).';

    snc = sc*nc;
    snr = sr*nr;
    slr_snr = (slr-1)*snr;
    slc_snc = (slc-1)*snc;

    for row_iter = 1:N1
        for column_iter = 1:N2
            count_block = 0;
            minr = min(nr, row_iter);
            minc = min(sr, column_iter);
            ri = (minc-1)*nr + minr;
            ci = (column_iter - minc)*nc + row_iter - minr + 1;
            AB_part = zeros(slr, slc);

            row_array = ri : snr : ri+slr_snr;
            column_array = ci : snc : ci+slc_snc;

            while row_array(1)>0 && column_array(end)<=slc*snc && count_block<3
                A_part = A(row_array,:);
                B_part = B(:, column_array);
                AB_temp = A_part * B_part;
                AB_part = AB_part + AB_temp;
                row_array = row_array - nr;
                column_array = column_array + nc;
                count_block = count_block + 1;
            end

            AB_part = AB_part*w(row_iter)*sw(column_iter)/count_block;
            HhAB(row_iter, column_iter, :) = Hankel2vec(AB_part);

        end
    end
end

