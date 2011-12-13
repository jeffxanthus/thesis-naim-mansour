% tomlab.m calls the new license file tomlablic

function [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12] = tomlab(varargin)

if length(varargin) == 0
   tomlablic;
elseif length(varargin) == 1
   if varargin{1} == 10
      [y1] = tomlablic(varargin{:});
   else
      [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12] = tomlablic(varargin{:});
   end
else
   if exist('tomlablic') == 3
      [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12] = tomlablic(varargin{:});
   else
      error('License file tomlablic does not exist');
   end
end