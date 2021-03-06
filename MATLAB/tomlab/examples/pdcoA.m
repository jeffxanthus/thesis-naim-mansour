function y = pdcoA( mode, m, n, x )

%        y = pdcoA( mode, m, n, x )
%    computes y = Ax (mode=1) or A'x (mode=2)
%    for a matrix A defined by the array in global pdAAA
%    (which is generated by pdcotest2.m).
%
%    Aname will always be the string 'pdcoA'.

%-----------------------------------------------------------------------
% 16 Oct 2002: First version of pdcoA.m for use by pdcotest2.m
%              to test defining A as a function.
%-----------------------------------------------------------------------


%function y = pdcoA( Aname, mode, m, n, x )
  global pdAAA

  if mode==1,  y = pdAAA*x;  else  y = pdAAA'*x;  end
  
%-----------------------------------------------------------------------
% End function pdcoA
%-----------------------------------------------------------------------