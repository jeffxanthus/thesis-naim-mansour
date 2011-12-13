%% UIMP - Production Scheduling Problem
% TomSym implementation of GAMS Example (UIMP,SEQ=11)
%
% A company manufactures nuts, bolts and washers using three different
% machines that can be operated in normal or overtime production mode. The
% company needs to plan operations for the next two periods.
%
% Ellison, E F D, and Mitra, P, UIMP - User Interface for Mathematical
% Programming. ACM Transactions on Mathematical Software 8, 2 (1982).
%
% i: time periods (summer, winter)
% j: production mode (normal, overtime)
% k: products (nuts, bolts, washers)
% l: machines (m1, m2, m3)

% Machine hours (hours per unit), mh(l,k)
mh = [ 4  4  6
       7  6  6
       3  0  0 ];

% Addfactors for mh(i,j) 
mhadd = [ 0  -1 
          1   0 ];

%  Availability  (hours)
av = [ 100  80 
       100  90 
        40  30 ];

% Machine hours required
t = zeros(2,2,3,3);
for i=1:2
    for j=1:2
        for k=1:3
            for l=1:3
                t(i,j,k,l) = mh(l,k);
                if mh(l,k) > 0
                    t(i,j,k,l) = t(i,j,k,l)+mhadd(i,j);
                end
            end
        end
    end
end
t(2,2,3,1) = 5;

% Machine hours available;
a = zeros(2,2,3);
a(1,:,:) = av';
a(2,:,:) = av'+10;

% Production cost data
tc = [ 2 3 4 
       4 3 2
       1 0 0 ];

% Addfactors for tc
tcadd = [ 0 1
          1 2 ];

% Production cost
c = zeros(2,2,3,3);
for i=1:2
    for j=1:2
        for k=1:3
            for l=1:3
                c(i,j,k,l) = tc(l,k);
                if tc(l,k) > 0
                    c(i,j,k,l) = c(i,j,k,l)+tcadd(i,j);
                end
            end
        end
    end
end

% Selling price
p = [ 10 10  9
      11 11 10 ];

% Demand
d = [ 25 30 30
      30 25 25 ];

% Storage cost
s = ones(3,1);

% Storage capacity
h = [20;20;0];

% Production
toms 3x3 x11 x12 x21 x22

% Products stored
toms 2x3 y

% Products sold
toms 2x3 z

%  Positive Variables: x, y
cbnd = {x11 >= 0, x12 >= 0, x21 >= 0, x22 >= 0, y >= 0};

% m3 can only produce nuts
cbnd = {cbnd{:}, x11(2:3,3) <= 0, x12(2:3,3) <= 0, ...
    x21(2:3,3) <= 0, x22(2:3,3) <= 0};

cost = 0;
for i=1:2
    for k=1:3
        cost = cost + s(k)*y(i,k);
        if i == 1
            cost = cost + sum(squeeze(c(1,1,k,:)).*x11(k,:)'+...
                squeeze(c(1,2,k,:)).*x12(k,:)');
        else
            cost = cost + sum(squeeze(c(2,1,k,:)).*x21(k,:)'+...
                squeeze(c(2,2,k,:)).*x22(k,:)');
        end
    end
end

revenue = sum(sum(p.*z));

profit = revenue - cost;

eq1 = {};
for l=1:3
    eq1 = {eq1{:}, sum(squeeze(t(1,1,:,l)).*x11(:,l)) <= a(1,1,l)};
    eq1 = {eq1{:}, sum(squeeze(t(1,2,:,l)).*x12(:,l)) <= a(1,2,l)};
    eq1 = {eq1{:}, sum(squeeze(t(2,1,:,l)).*x21(:,l)) <= a(2,1,l)};
    eq1 = {eq1{:}, sum(squeeze(t(2,2,:,l)).*x22(:,l)) <= a(2,2,l)};
end

eq2 = {};
for i=1:2
    for k=1:3
        if i>1
            eq2 = {eq2; sum(x21(k,:)+x22(k,:))...
                + y(i-1,k) == z(i,k) + y(i,k)};
        else
            eq2 = {eq2; sum(x11(k,:)+x12(k,:))...
                == z(i,k) + y(i,k)};
        end
    end
end

eq3 = {z >= d; y <= [h'; h']};

options = struct;
options.scale = 'manual';
options.name = 'max_revenue';
solution1 = ezsolve(-revenue,{cbnd, eq1, eq2, eq3},[],options);
options.name = 'max_profit';
solution = ezsolve(-profit,{cbnd, eq1, eq2, eq3},[],options);

disp(' ');
disp('Maximum revenue:');
disp(subs(revenue,solution1));

disp('Maximum profit:');
disp(subs(profit,solution));