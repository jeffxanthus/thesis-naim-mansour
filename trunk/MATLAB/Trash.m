function [a,b] = Trash()
%TRASH 
frameOverlap=50;
numberOfSamples=882;
position=(1-(50/100))*882+1;
i=2;
while position+(numberOfSamples)<=9713
    position=position+(1-(frameOverlap/100))*numberOfSamples
    i=i+1
    pause
end
position

end

