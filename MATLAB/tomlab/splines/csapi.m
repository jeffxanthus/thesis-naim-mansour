%CSAPI Cubic spline interpolant with not-a-knot end condition.
%
%   VALUES  = CSAPI(X,Y,XX)
%
%   returns the values at XX of the cubic spline interpolant to the
%   given data (X,Y) (using the not-a-knot end condition).  The
%   ordinates may be vector-valued, in which case Y(:,j) is the j-th
%   ordinate.
%
%   The alternative call
%
%   PP  = CSAPI(X,Y)
%
%   returns the ppform of the cubic spline interpolant instead, for
%   later use with FNVAL, FNDER, FNPLT, etc.
%
%   For example, 
%
%      values = csapi([-1:5]*(pi/2),[0 1 0 -1 0], linspace(0,2*pi));
%
%   gives a surprisingly good fine sequence of values for the sine over its
%   period.
%
%   It is also possible to interpolate to data on a rectangular grid and
%   obtain values of the interpolant on a rectangular grid or at scattered
%   points, by the calls
%
%   VALUES = CSAPI({X1, ...,Xm},Y,XX)
%   or
%   PP = CSAPI({X1, ...,Xm},Y)
%
%   in which Y is expected to have size [d,length(X1),...,.length(Xm)]
%   (or [length(X1),...,.length(Xm)] if the function is to be scalar-valued),
%   and XX is either a list of m-vectors XX(:,j) or else a cell-array 
%   {XX1,...,XXm} specifying the m-dimensional grid at which to evaluate
%   the interpolant.
%
%   See also CSAPE, SPAPI, SPLINE.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written June 20, 1999. Last modified Aug 28, 2000.
%

function output = csapi(x,y,xx)

if nargin <3
   xx=[];
   if nargin <2
      y=[];
      if nargin <1
         x=[];
      end
   end
end

disp('Dummy routine csapi for TOMLAB');
pause

output=[];

