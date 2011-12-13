function [ bark ] = FreqToBark( f )
%FreqToBark Converts from frequency to bark

bark = 13*arctanApprox(0.00076*f) + 3.5*arctanApprox((f/7500)^2);

end

