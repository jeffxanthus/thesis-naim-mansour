% function f = lgo2_f(x, Prob)
%
% Test functions for global optimization. Two or more dimensions.

% Reference:
% Pintér, J.D., Bagirov, A., and Zhang, J. (2003) An Illustrated Collection of
% Global Optimization Test Problems. Research Report, Pintér Consulting Services,
% Inc. Halifax, NS, Canada; and CIAO-ITMS, University of Ballarat, Ballarat, Vic.,
% Australia.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function f = lgo2_f(x, Prob)

x=x(:);

P=Prob.P;
x1=x(1);
x2=x(2);

if P == 1 % M1
   f = sin(x1)*(cos(x1)-sin(x1))^2+sin(x2)*(cos(x2)-sin(x2))^2;
elseif P == 2 % M2
   f = sin(3*x1)+sin(5*x1)+sin(7*x1)+sin(11*x1)+sin(3*x2)+sin(5*x2)+sin(7*x2)+sin(11*x2);
elseif P == 3 % M3
   f = -20*exp(-0.02*sqrt(0.5*(x1^2+x2^2)))-exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2)))+20+exp(1);
elseif P == 4 % M4
   f = x2^3-x1^3+50*x1^2+25*x2^2-25*x1*x2+100*x2;
elseif P == 5 % M5
   f = (abs(x1)-5)^2+(abs(x2)-5)^2;
elseif P == 6 % M6
   f = x1^2+2*x2^2-3*cos(2*pi*x1)-4*cos(4*pi*x2)+7;
elseif P == 7 % M7
   f = x1^2+2*x2^2-5*cos(2*pi*x1)*cos(4*pi*x2)+5;
elseif P == 8 % M8
   f = x1^2+2*x2^2-cos(3*pi*x1+4*pi*x2)+1;
elseif P == 9 % M9
   f = (x2-5.1*x1^2/4/pi^2+5*x1/pi-6)^2+10*(1-1/8/pi)*cos(x1)+10;
elseif P == 10 % M10
   f = (-x1-2*x2+1/20*sin(4*pi*x2)+1)^2+(x2-1/2*sin(2*pi*x1))^2;
elseif P == 11 % M11
   f = 1/6*x1^6-1.05*x1^4+2*x1^2+x1*x2+x2^2;
elseif P == 12 % M12
   f = 1/3*x1^6-2.1*x1^4+4*x1^2+x1*x2+4*x2^4-4*x2^2;
elseif P == 13 % M13
   f = 10*(x1^2-x2)^2+(x1-1)^2;
elseif P == 14 % M14
   f = -cos(x1)*cos(x2)*exp(-(x1-pi)^2-(x2-pi)^2);
elseif P == 15 % M15
   f = -sin(x1*cos(pi/6)-x2*sin(pi/6))*(sin((x1*cos(pi/6)-x2*sin(pi/6))^2/pi))^2-sin(x2)*(sin(2*x2^2/pi))^20;
elseif P == 16 % M16
   f = max(abs(x1+x2/2),abs(x1/2+x2/3));
elseif P == 17 % M17
   f = ((3*x1^2+6*x2*x1-14*x1+3*x2^2-14*x2+19)*(x1+x2+1)^2+1) * ((12*x1^2-36*x2*x1-32*x1+27*x2^2+48*x2+18)*(2*x1-3*x2)^2+30);
elseif P == 18 % M18
   f = (x1^2+x2^2)/4000-cos(x1)*cos(x2/sqrt(2))+1;
elseif P == 19 % M19
   f1=0;
   f2=0;
   for i=1:5
      f1 = f1 + i*cos(i+(i-1)*x1);
      f2 = f2 + i*cos(i+(i-1)*x2);
   end
   f = f1*f2;
elseif P == 20 % M20
   f = (x1^4/4-7*x1^3/3+7*x1^2-8*x1+1)*x2^2*exp(-x2);
elseif P == 21 % M21
   f = 1/2*pi*((10*(sin(pi*(x2-1)/4+1))^2+1)*(((x1-1)/4+1)-1)^2 + ((x2-1)/4+1) + 10*sin(pi*((x1-1)/4+1))*sin(pi*((x1+1)/4+1))-1)^2;
elseif P == 22 % M22
   f = sin(3*pi*x2)^2 + (x1-1)^2*sin(3*pi*x2)^2 + (x2-1)^2*(sin(2*pi*x2)^2+1);
elseif P == 23 % M23
   f = 1/2*pi*((10*(sin(pi*((x2+1)/4+1))^2+1)+1) * (((x1+1)/4+1)-1)^2 + (((x2+1)/4+1)-1)^2 + 10*(sin(pi*((x1+1)/4+1))^2));
elseif P == 24 % M24
   f = -(x1*sin(4*pi*x1)+x2*sin(20*pi*x2)+21.5);
elseif P == 25 % M25
   f = -(0.5*exp(-(x1^2+x2^2)/0.1^2) + 1.2*exp(-((x2-1)^2+x1^2)/0.5^2) + 1.2*exp(-((x1-1)^2+x2^2)/0.5^2) + 1*exp(-(x2^2+(x1+0.5)^2)/0.5^2) + 1*exp(-(x1^2+(x2+0.5)^2)/0.5^2));
elseif P == 26 % M26
   f = (x1-1)^2+(x2-1)^2-x1*x2;
elseif P == 27 % M27
   f = -(0.2*sqrt((x2-1.3)^2+(x1-1)^2)/(max(abs(x1-1),abs(x2-1.3)) * sqrt(2) + 0.1) + 1.0) * cos(pi*sqrt(2)*max(abs(x1-1),abs(x2-1.3))) * exp(-sqrt(2)*max(abs(x1-1),abs(x2-1.3))/2/pi);
elseif P == 28 % M28
   f = log(x1-2)^2 + log(x2-2)^2 + log(10-x1)^2 + log(10-x2)^2 - (x1*x2)^0.2;
elseif P == 29 % M29
   f = 1 + sin(x1)^2 + sin(x2)^2 - exp(-x1^2-x2^2);
elseif P == 30 % M30
   f = 0;
   for i=1:5
      f = f + ((x1-1)+0.01*(x2-1)*i)^2;
   end
elseif P == 31 % M31
   f = max([x1^2+x2^4,(2-x1)^2+(2-x2)^2,2*exp(x2-x1)]);
elseif P == 32 % M32
   f = max([x1^4+x2^2,(2-x1)^2+(2-x2)^2,2*exp(x2-x1)]);
elseif P == 33 % M33
   f = max([x1^2+x2^2, x1^2+x2^2+10*(-4*x1-x2+4), x1^2+x2^2+10*(-x1-2*x2+6)]);
elseif P == 34 % M34
   x3 = x(3); 
   x4 = x(4);
   f1 = x1^2-5*x1+x2^2+2*x3^2+x4^2-5*x2-21*x3+7*x4;
   f2 = x1^2+x1+x2^2+x3^2+x4^2-x2+x3-x4-8;
   f3 = x1^2-x1+2*x2^2+x3^2+2*x4^2-x4-10;
   f4 = x1^2+2*x1+x2^2+x3^2-x2-x4-5;
   f = max([f1,f1+10*f2,f1+10*f3,f1+10*f4]);
elseif P == 35 % M35
   f = 0;
   for i = 1:100
      f = f + abs(x1-0.5)+0.01*i*abs(x2-0.5);
   end
elseif P == 36 % M36
   f = abs(x1-1)+100*abs(x2-abs(x1));
elseif P == 37 % M37
   x3 = x(3);
   x4 = x(4);
   f = abs(x1-1)+abs(x3-1)+10.1*(abs(x2-1)+abs(x4-1))+4.95*(abs(x2+x4-2)-abs(x2-x4))+100*abs(x2-abs(x1))+90*abs(x4-abs(x3));
elseif P == 38 % M38
   f = 0;
   for i=1:10
       f = f + abs(x1-0.5+0.01*i*(x2-0.5));
   end
   fe = max(abs(x1-0.5+0.01*(1:100)*(x2-0.5)));
   f = f - fe;
elseif P == 39 % M39
   f = x1^2+x2^2-cos(18*x1)-cos(18*x2);
elseif P == 40 % M40
   f = x1^2-10*cos(2*pi*x1)+x2^2-10*cos(2*pi*x2)+20;
elseif P == 41 % M41
   f = (x1+x2-2)^2/(x1+3*x2)+(3*x1-x2)^2/(2*x1+x2);
elseif P == 42 % M42
   f = 2*sin(2*x1)^2+3*sin(3*x2)^2;
elseif P == 43 % M43
   f = -cos(2*pi*sqrt(x1^2+x2^2))+0.1*sqrt(x1^2+x2^2)+1;
elseif P == 44 % M44
   f = (sin(x1^2+x2^2)^2-0.5)/(0.001*(x1^2+x2^2)+1)^2+0.5;
elseif P == 45 % M45
   f = (x1^2+x2^2)^0.25*floor(sin(50*(x1^2+x2^2)^0.1)^2+1);
elseif P == 46 % M46
   f1 = 0;
   f2 = 0;
   for i=1:5
       f1 = f1 + i*cos(i+(i+1)*x1);
       f2 = f2 + i*cos(i+(i+1)*x2);
   end
   f = f1*f2;
elseif P == 47 % M47
   f = 0.5*((x2+0.80032)^2+(x1+1.42513)^2);
   f1 = 0;
   f2 = 0;
   for i=1:5
       f1 = f1 + i*cos(i+(i+1)*x1);
       f2 = f2 + i*cos(i+(i+1)*x2);
   end
   f = f + f1*f2;
elseif P == 48 % M48
   f1 = 0;
   f2 = 0;
   for i=1:5
       f1 = f1 + i*sin(i+(i+1)*x1);
       f2 = f2 + i*sin(i+(i+1)*x2);
   end
   f = -f1-f2;
elseif P == 49 % M49
   f = -(x1*sin(sqrt(abs(x1))) + x2*sin(sqrt(abs(x2))));
elseif P == 50 % M50
   f = 2.5*sin(x1-30)*sin(x2-30)+sin(5*(x1-30))*sin(5*(x2-30));
elseif P == 51 % M51
    f = 1/((x2-9.112)^2+(x1-4.887)^2+0.811);
    f = f + 1/((x2-8.858)^2+(x1-0.226)^2+0.513);
    f = f + 1/((x2-8.777)^2+(x1-8.074)^2+0.965);
    f = f + 1/((x2-8.645)^2+(x1-0.432)^2+0.817);
    f = f + 1/((x2-8.583)^2+(x1-6.306)^2+0.964);
    f = f + 1/((x2-8.057)^2+(x1-1.460)^2+0.332);
    f = f + 1/((x1-8.304)^2+(x2-7.559)^2+0.869);
    f = f + 1/((x2-7.549)^2+(x1-3.352)^2+0.369);
    f = f + 1/((x2-7.027)^2+(x1-0.652)^2+0.462);  
    f = f + 1/((x2-7.006)^2+(x1-2.132)^2+0.714);  
    f = f + 1/((x2-6.686)^2+(x1-2.440)^2+0.828);  
    f = f + 1/((x1-7.650)^2+(x2-5.658)^2+0.669);  
    f = f + 1/((x2-5.579)^2+(x1-4.707)^2+0.352);  
    f = f + 1/((x1-9.496)^2+(x2-4.830)^2+0.608);  
    f = f + 1/((x1-8.632)^2+(x2-4.409)^2+0.813);  
    f = f + 1/((x1-8.327)^2+(x2-3.897)^2+0.463);  
    f = f + 1/((x2-3.605)^2+(x1-1.256)^2+0.524);
    f = f + 1/((x2-3.516)^2+(x1-2.699)^2+0.491);    
    f = f + 1/((x2-2.800)^2+(x1-0.679)^2+0.632);    
    f = f + 1/((x1-4.138)^2+(x2-2.562)^2+0.326);    
    f = f + 1/((x2-2.343)^2+(x1-0.652)^2+0.789);    
    f = f + 1/((x1-8.314)^2+(x2-2.261)^2+0.902);    
    f = f + 1/((x1-7.305)^2+(x2-2.228)^2+0.876);    
    f = f + 1/((x1-9.400)^2+(x2-2.041)^2+0.517);    
    f = f + 1/((x1-5.558)^2+(x2-1.272)^2+0.360);    
    f = f + 1/((x1-4.263)^2+(x2-1.074)^2+0.883);    
    f = f + 1/((x1-8.798)^2+(x2-0.880)^2+0.992);    
    f = f + 1/((x1-9.681)^2+(x2-0.667)^2+0.806);    
    f = f + 1/((x1-2.196)^2+(x2-0.415)^2+0.908);    
    f = f + 1/((x2-9.152)^2+(x1-8.025)^2+0.100);
    f = -f;
elseif P == 52 % M52
   f = 1/4*(x1^2+x2^2)+exp(sin(50*x1))-sin(10*(x1+x2))+sin(60*exp(x2))+sin(70*sin(x1))+sin(sin(80*x2));
elseif P == 53 % M53
   f = sin(x1+x2)^2-cos(5*x1+3*x2)+sin(2*x1-5*x2)^2-cos(x1-7*x2);
end

% MODIFICATION LOG
%
% 040114  med  Created
% 040121  med  Added M37-M53
% 080603  med  Switched to glcAssign, cleaned