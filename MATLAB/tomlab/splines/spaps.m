%SPAPS Smoothing spline.
%
%   [SP,VALUES] = SPAPS(X,Y,TOL)  returns the B-form and, if asked, the values 
%   at  X , of the smoothest d-valued function  f , in the sense that
%
%       F(D^M f) := sum integral { (D^M f(t))^2 : X(1) < t < X(n) }
%
%   is smallest, for which
%
%       E(f) :=  sum sum_j {W(j)*(Y(:,j) - f(X(j)))^2 : j=1,...,n }
%
%   is no bigger than TOL, with the weights W chosen so that  E(f)  is
%   the Composite Trapezoid Rule approximation to  F(y - f), and with 
%   M = 2,   i.e.,  f  is a  c u b i c  smoothing spline. This is so because
%   f is constructed as the unique minimizer of
%                     rho*E(f) + F(D^2 f) ,
%   with the smoothing parameter RHO (optionally returned as a third output
%   argument) so chosen that  E(f) == TOL . Hence, after conversion to ppform
%   (and up to roundoff), the result should be the same as obtained by
%   cpaps(x,y,RHO/(1+RHO)).
%
%   Here, the default for d is 1 and, for that case, the extra summation in
%   the definition of F and E(f) does nothing. However, if Y is not just a 
%   vector, then Y is expected to have length(X) columns, and then d is the 
%   number of rows in Y, and, correspondingly, the output is a d-vector-valued 
%   function, and the sum in the definition of F and the first sum in the 
%   definition of E(f) is the sum over the d rows.  
%   See below for the case of  g r i d d e d  data.
%
%   SPAPS(X,Y,TOL,W)
%   SPAPS(X,Y,TOL,M)
%   SPAPS(X,Y,TOL,W,M)
%   SPAPS(X,Y,TOL,M,W)
%   provide W and/or M explicitly, with M presently restricted to
%   M = 1 (linear) and M = 3 (quintic).
%
%   If X is not increasing, both X and Y will be reordered in unison to
%   make it so. After that, X must be  s t r i c t l y  increasing.
%   Also, Y must be of size d-by-length(X).
%
%   If X is a cell-array, of length  r , then Y is expected to supply
%   correspondingly gridded data, with Y of size [length(X{i}): i=1:r]
%   for scalar data, and of size [ d, [length(X{i}): i=1:r] ]  for 
%   d-vector-valued data. In this case of gridded data, the optional
%   argument M is either an integer or else an r-vector of integers 
%   (from the set  {1,2,3} ), and the optional argument W is expected
%   to be a cell array of length r with W{i} either empty (to get the
%   default choice) or else a vector of the same length as X{i}, i=1:r.
%
%   For example, the statements
%
%      w = ones(size(x)); w([1 end]) = 100;
%      sp = spaps(x,y,1.e-2,w,3);
%
%   give a quintic smoothing spline approximation to the given data which
%   close-to-interpolates the first and last datum, while being within about
%   1.e-2 of the rest.
%
%   See also SPAPI, SPAP2, CSAPS.


% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written June 20, 1999. Last modified Aug 28, 2000.
%

function [sp,values,rho] = spaps(x,y,tol,varargin)

if nargin <3
   tol=[];
   if nargin <2
      y=[];
      if nargin <1
         x=[];
      end
   end
end

disp('Dummy routine spaps for TOMLAB');

sp=[];
values=[];
rho=[];

