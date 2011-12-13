% glcQG_c - constraints for Hock-Schittkowski 59, Quickguide
%
% function cx = glcQG_c(x, Prob)

function cx = glcQG_c(x, Prob)

cx = [x(1)*x(2)-700;x(2)-x(1)^2/125;(x(2)-50)^2-5*(x(1)-55)];