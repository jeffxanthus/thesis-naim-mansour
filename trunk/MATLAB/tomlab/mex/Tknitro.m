% For information on how to use TOMLAB /KNITRO, see help knitroTL.m.
%
% Call to MEX file Tknitro:
%
% [Inform,iter,x_k,f_k,g_k,c_k,v_k,cJac,H_k] = ...
%     Tknitro(n,mL,mNL,xl,xu,x_0,A,cl,cu,ConsPattern,nz,hesscol, ...
%             hessrow,f_Low,PriLev,PrintFile,options,cb,Prob);

%# mex

help Tknitro