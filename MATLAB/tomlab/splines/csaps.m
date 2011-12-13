%CSAPS Cubic smoothing spline.
%
%   CSAPS(X,Y,P)  returns the ppform of the cubic smoothing spline for the
%   given data (X,Y)  and depending on the smoothing parameter P from 
%   [0 .. 1].  This smoothing spline  f  minimizes
%
%   P * sum_i W(i)(Y(i) - f(X(i)))^2  +  (1-P) * integral (D^2 f)^2
%
%   with W = ones(length(X),1) the default value for  W .
%   For  P=0, the smoothing spline is the least-squares straight line fit to
%   the data, while, at the other extreme, i.e., for  P=1, it is
%   the `natural' or variational cubic spline interpolant.  The
%   transition region between these two extremes is usually only a
%   rather small range of values for P and its location strongly
%   depends on the data.
%   If you have difficulty choosing  P , but have some feeling for the size
%   of the noise in  Y , consider using instead  spaps(X,Y,tol)  which, in
%   effect, chooses  P  in such a way that
%           sum_i W(i)(Y(i) - f(X(i)))^2  =  tol  .
%
%   CSAPS(X,Y,P,XX) returns instead the value(s) at XX of the cubic smoothing
%   spline, unless XX is empty, in which case the ppform of the cubic smoothing
%   spline is returned. This latter option is important when the user wants
%   to supply the weights W, as in the following:
%
%   CSAPS(X,Y,P,XX,W) returns, depending on whether or not XX is empty, the
%   ppform, or the values at XX,  of the cubic smoothing spline for the 
%   specified weights W.
%   
%   For example, 
%
%      x = linspace(0,2*pi,21); y = sin(x)+(rand(1,21)-.5)*.1;
%      pp = csaps(x,y,.4,[],[ones(1,11), repmat(5,1,10)]);
%
%   returns a smoothed version of the data that is much closer to the data
%   in the right half, because of the much larger weight there.
%
%   It is in general difficult to choose the parameter P without
%   experimentation. For that reason, use of SPAPS is recommended instead
%   since there P is chosen so as to produce the smoothest spline within a
%   specified tolerance of the data.
%
%   It is also possible to smooth data on a rectangular grid and
%   obtain smoothed values on a rectangular grid or at scattered
%   points, by the calls
%   VALUES = CSAPS({X1, ...,Xm},Y,P,XX)
%   VALUES = CSAPS({X1, ...,Xm},Y,P,XX,W)
%   or the ppform of the corresponding cubic smoothing spline by the calls
%   PP = CSAPS({X1, ...,Xm},Y,P)
%   PP = CSAPS({X1, ...,Xm},Y,P,[],W)
%   in which Y is expected to have size [d,length(X1),...,.length(Xm)]
%   (or [length(X1),...,.length(Xm)] if the function is to be scalar-valued),
%   and P is either a scalar or an m-vector,
%   and XX is either a list of m-vectors XX(:,j) or else a cell-array 
%   {XX1, ..., XXm} specifying the m-dimensional grid at which to evaluate
%   the interpolant, and, correspondingly, W, if given, is a cell array of
%   weight sequences for the m dimensions (with an empty Wi indicating the
%   default choice).
%
%   See also SPAPS, CSAPSDEM.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written June 20, 1999. Last modified Aug 28, 2000.
%

function output = csaps(x,y,p,xx,w)

if nargin <5
   w=[];
   if nargin <4
      xx=[];
      if nargin <3
         p=[];
         if nargin <2
            y=[];
            if nargin <1
               x=[];
            end
         end
      end
   end
end

disp('Dummy routine csaps for TOMLAB');

output=[];


