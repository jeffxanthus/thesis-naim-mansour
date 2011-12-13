% tomlic.m
%
% TOMLAB license information utility. Calls tomlablic(80). 
%

try
  tomlablic(80)
catch
  error('Failed to invoke tomlablic(80) to get license information')
end
