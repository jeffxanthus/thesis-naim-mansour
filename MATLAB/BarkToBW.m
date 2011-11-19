function [ freq ] = BarkToBW( b )
%BarkToFreq converts from bark to bandwidth (Hz)

freq = 52548/(b^2-52.56*b+690.39);
end

