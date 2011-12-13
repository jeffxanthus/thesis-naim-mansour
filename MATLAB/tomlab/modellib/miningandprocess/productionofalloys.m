% function Prob = productionofalloys(compsize, mincomp, maxcomp,
% rawcompmat, rawavail, rawcost)
%
% Creates a TOMLAB LP problem for production design of alloys
%
% INPUT PARAMETERS
% compsize     Size of component order (in generic unit 1).
% mincomp      Lower limit on specification. Vector. (in generic unit 2).
% maxcomp      Upper limit on specification. Vector. (in generic unit 2).
%
% rawcompmat   Component matrix for raw materials (in generic unit 2).
% rawavail     Availability (in general unit 1).
% rawcost      Cost (in generic currency / general unit 1)
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type LP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Dec 7, 2005.

function Prob = productionofalloys(compsize, mincomp, maxcomp, rawcompmat,...
   rawavail, rawcost)

if nargin < 6
   error('The function requires 6 inputs');
end

if isempty(compsize) | isempty(mincomp) | isempty(maxcomp) | isempty(rawcompmat) |...
       isempty(rawavail) | isempty(rawcost)
    error('One of the inputs are empty');
end

mincomp = mincomp(:); % lower bounds on components
maxcomp = maxcomp(:); % upper bounds on components

m = length(mincomp);

if m~=length(maxcomp) | m~=size(rawcompmat,2)
   error('Incorrect sizes for either maxcomp or rawcompmat (columns should equal rows in maxcomp');
end

rawavail = rawavail(:);
rawcost = rawcost(:);

n = size(rawcompmat,1);

if n~=length(rawavail) | n~=length(rawcost)
   error('Incorrect sizes for either rawavail or rawcost');
end

% FORMULATE PROBLEM

c   = rawcost;
x_L = zeros(n,1);
x_U = rawavail;

% Production constraint

% A1   = ones(1,n);
% b_L1 = compsize;
% b_U1 = compsize;

% Component constraint

% A2   = rawcompmat';
% b_L2 = mincomp*compsize;
% b_U2 = maxcomp*compsize;

% Production and component constraint together

A = [sparse(ones(1,n));rawcompmat'];
b_L = [compsize;mincomp*compsize];
b_U = [compsize;maxcomp*compsize];

Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production of Alloys');

% MODIFICATION LOG
%
% 051007 med   Created.
% 051207 med   isscalar removed.