% For information on how to use TOMLAB /CONOPT, see help conoptTL.m.
%
% Call to MEX file conopt:
%
% [modsta,solsta,iter,f_k,x_k,xmar,xbas,xsta,yval,ymar,ybas,ysta,vrsn] = ...
%    conopt(n,mL,mNL,x_0,xl,xu,A,b,btype,nza,c,ctype,ConsPattern,nzc,...
%           nHess,hessrow,hesscol,DebugFV,PriLev,PrintFile,StatFile,OptFile,...
%           options,Prob) ;

%# mex

help conopt