function H = mpd_H(x,Prob);

% H = [ 2 0 0 0 0 ; 0 8 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 0 ];

H = zeros(5,5);
H(1,1) = 2.0;
H(2,2) = 8.0;