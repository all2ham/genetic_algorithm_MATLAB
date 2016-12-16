 function [dec]=bintodec_multi(bin,dim,bits)
    % modified version of bintodec to handle multivariable bitstrings
    global range;
    % Length of the string without signs
    % nn=length(bin)-dim;
    % get the binaries
    for i = 1:dim
        dtemp(i,:) = bin((bits*(i-1)+i+1):(bits*i+i));
        signs(i) = 1-2.*bin(bits*(i-1)+i);
    end
    % Sign=+1 if bin(l)=0; Sign=-l if bin(l)=l.
    dec=zeros(1,dim);
    % floating point/decimal place in a binary string
    dp=floor(log2(max(abs(range))));
    for j = 1:length(dec)
        for i=1:bits,
            dec(j)=dec(j)+dtemp(j,i)*2^(dp-i);
        end
    end
    dec=dec.*signs;
 end

