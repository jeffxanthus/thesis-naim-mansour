function checkarch

r = {'mexglx','mexa64','mexmaci','mexmaci64','mexsol','mexw64','mexw32','dll'};
v = {'Linux 32-bit','Linux 64-bit', 'Mac OS (32-bit)', 'Mac OS (64-bit)',...
  'Solaris 64-bit', 'Windows 64-bit', 'Windows 32-bit R2007A and newer', ...
  'Windows 32-bit R2006A and older'};
c = computer;

if length(which('tomlablic','-ALL')) > 2
  error('Error: Two versions of TOMLAB are in the Matlab path.');
end

t = which('tomlablic');
[p,f,e] = fileparts(t);

% Special case (replaces tlmxarch)
if isequal(c,'PCWIN')
  verstr = version;
  if ~isempty(findstr(verstr,'R201'))
    old = 0; % 201x is a "new" Matlab.
  else
    V = str2num(verstr(1:3));
    if V<=7.3
      old=1;
    else
      old=0;
    end
  end
  if isequal(e,'.mexw32') & old
    err = sprintf('The version of TOMLAB installed at\n\n %s \n\nis for Matlab version 7.4 and newer but you are running version\n\n%s.\n\n',p,verstr);
    err = [err 'Please go to http://tomopt.com/tomlab/ and download the correct version of TOMLAB' ];
    % Try to show error dialog (display might not exist)
    try
      errordlg(err,'TOMLAB Installation Error');
    catch
    end
    error(err);
  elseif isequal(e,'.dll') & ~old
    err=sprintf('The version of TOMLAB installed at\n\n%s\n\nis for Matlab version 7.3 and older but you are running version\n\n%s.\n\n',p,verstr);
    err = [err 'Please go to http://tomopt.com/tomlab/ and download the correct version of TOMLAB' ];
    
    % Try to show error dialog (display might not exist)
    try
      errordlg(err,'TOMLAB Installation Error');
    catch
    end
    error(err);
  end
end

%%
% Check if there's a tomlablic.* with wrong mex extension
if isempty(t)
  p = pwd;
  for k=1:length(r)
    F = ['tomlablic.' r{k} ];
    if exist(F,'file')
      str = sprintf(['TOMLAB''s startup routine has detected a possible installation error. \n' ...
        'The TOMLAB installation at\n\n%s\n\ncontains\n\n%s\n\nwhich indicates a ' ...
        'TOMLAB installation for:\n\n%s'], p,F,v{k});
      
      % Try to show error dialog (display might not exist)
      try
        errordlg(str,'TOMLAB installation error');
      catch
      end
      error(str);
    end
  end
end
