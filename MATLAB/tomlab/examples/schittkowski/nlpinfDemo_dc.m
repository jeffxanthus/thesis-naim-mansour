function dcx = nlpinfDemo_dc(x, Prob)

if Prob.USER.P == 1
   dcx = [];
   disp(Prob.P);
elseif Prob.USER.P == 2
   t = Prob.LS.t;
   ix = [1 length(t)];
   dcx(:,1) = t(ix).*(t(ix) + x(2))./(t(ix).^2 + x(3)*t(ix) + x(4));
   dcx(:,2) = x(1)*t(ix)./(t(ix).^2 + x(3)*t(ix) + x(4));
   dcx(:,3) = -x(1)*(t(ix).^2).*(t(ix) + x(2))./(t(ix).^2 + x(3)*t(ix) + x(4)).^2;
   dcx(:,4) = -x(1)*t(ix).*(t(ix) + x(2))./(t(ix).^2 + x(3)*t(ix) + x(4)).^2;
end
