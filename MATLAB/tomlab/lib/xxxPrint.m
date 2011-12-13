function xxxPrint(v,name)
for i=1:length(v)
    if length(v)==1
       fprintf('%s',name)
       fprintf(':')
    else
       fprintf('%s',name)
       fprintf(' #%2d:',i)
    end
    fprintf(' %40.20f',v(i));
    fprintf(' %bx',v(i));
    fprintf(' %x\n',v(i));
end