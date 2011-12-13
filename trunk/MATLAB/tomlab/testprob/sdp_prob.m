% sdp_prob:
%
% Defines Linear Semi-Definite Program-Integer Programming problems
%
% function [probList, Prob] = sdp_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function [probList, Prob] = sdp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'sdp.ps example 2'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

if P == 1
    Name = 'sdp.ps example 2';
    % Objective function
    c = [1 2 3]';
    % Two linear constraints
    A = [0 0 1; 5 6 0];
    b_L = [-Inf;-Inf];
    b_U = [ 3  ; -3 ];
    x_L = -1000*ones(3,1);
    x_U =  1000*ones(3,1);
    % Two linear matrix inequality constraints. It is OK to give only
    % the upper triangular part.
    SDP = [];
    % First constraint
    SDP(1).Q{1} = [2 -1 0; 0 2 0; 0 0 2];
    SDP(1).Q{2} = [2 0 -1; 0 2 0; 0 0 2];
    SDP(1).Qidx = [1; 3];
    % Second constraint
    SDP(2).Q{1} = diag([0  1]);
    SDP(2).Q{2} = diag([1 -1]);
    SDP(2).Q{3} = diag([3 -3]);
    SDP(2).Qidx = [0; 1; 2];
    x_0 = [];
else
    error('sdp_prob: Illegal problem number')
end

n = max( [length(c),size(A,2),length(x_L),length(x_U)] );
if isempty(x_0)
    x_0 = zeros(n,1);
end

Prob = sdpAssign(c, SDP, A, b_L, b_U, x_L, x_U, x_0, Name);
Prob.P = P;

% MODIFICATION LOG:
%
% 020701 hkh  Written (mip_prob)
% 030123 ango Major revision, remove mip_prob problems
% 030124 ango Add check for x_0
% 030127 ango Change LMI(i).Q0 --> LMI(i,1).Q0
% 041117 med  xxx_prob removed and code added
% 041214 frhe Converted problem to the new LMI format
% 080603 med  Switched to conAssign, cleaned