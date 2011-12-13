function testlsqcurvefit

fprintf('\n\n-----------------------------------------------------\n\n');   
fprintf('\n\nThree Least Squares test problems solved by lsqcurvefit\n\n');   
fprintf('\n\n-----------------------------------------------------\n\n');   
   
if 1
   Name='Walsh';
   fprintf('\n');   
   fprintf('Test of problem %s',Name);   
   fprintf('\n');   

   t=[ 2000; 5000; 10000; 20000; 30000; 50000];
   y=[0.9427; 0.8616; 0.7384; 0.5362; 0.3739; 0.3096];

   % C = askparam(ask, 'Give new C (Normal value 96.05) ', 90, 100, 96.05, uP);
   % uP(1)=C;

   % x_0=[-5;5];
   x_0=[-0.01;400];

   % Avoid division by 0
   % x_L=[-Inf; 0];
   % x_U=[   0; Inf];
   small = 1E-100;
   x_L=[-Inf;    small];
   x_U=[ -small;   Inf];
   if exist('optimset')
      options=optimset;
   else
      options=[];
   end
   % Define Tomlab problem structure to input the solver to be used
   % If the call is to the Opt Tbx lsqnonlin, it is simply not used
   global otxProb
   otxProb = ProbDef;
   otxProb.Solver.Tomlab = 'clsSolve';

   options.Display='iter';
   options.Diagnostics='on';

   % Define input needed by curve_rJ
   InProb.P  = 1;
   InProb.uP = []; % Use default value on parameter C (C=96.05)

   % Note that this structure is just additional user input
   % It is not recognized by Tomlab as a problem structure, because
   % the field TOMLAB is missing

   if nargout('lsqcurvefit')==8
      % Tomlab lsqcurvefit is returning the Result structure as well
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =lsqcurvefit...
         ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);
      % Using the Result structure the standard PrintResult routine is called
      disp(' ');
      disp('CALL TOMLAB PrintResult routine');
      disp(' ');
      PrintResult(Result,2);
   else
      options.Jacobian='on';
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k] =lsqcurvefit...
         ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);
   end

   Name='Gisela';
   fprintf('\n');   
   fprintf('\n\n-----------------------------------------------------\n\n');   
   fprintf('Test of problem %s',Name);   
   fprintf('\n');   
   fprintf('\n\n-----------------------------------------------------\n\n');   
   fprintf('\n');   


   t=[0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];
   y=[30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];
   % K = askparam(ask, 'Give new K (Normal value K=5) ', 1, 10, 5, uP);
   K=5;

   x_0=[6.8729,0.0108,0.1248]';
   x_L=[-Inf,-Inf,-Inf]';
   x_U=[Inf,Inf,Inf]';

   % Define input needed by lsq_rJ
   % Note that this structure is just additional user input
   % It is not recognized by Tomlab as a problem structure, because
   % the field TOMLAB is missing
   InProb.P=2;      % Set the problem number   
   InProb.uP(1)=K;  % Problem dependent parameter in standard place for user 
   InProb.Name=Name;% Set the problem name, nicer output in PrintResult later on

   % Define Tomlab problem structure to input the solver to be used
   % If the call is to the Opt Tbx lsqnonlin, it is simply not used
   global otxProb
   otxProb = ProbDef;
   otxProb.Solver.Tomlab = 'nlssol';

   % We can switch off the output from lsqnonlin, and instead use the
   % standard printing routine in TOMLAB (PrintResult)
   if exist('optimset')
      options=optimset;
   else
      options=[];
   end
   options.Diagnostics='on';
   options.Display='off';
   

   if nargout('lsqcurvefit')==8
      % Tomlab lsqnonlin is returning the Result structure as well
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =lsqcurvefit...
          ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);
      % Using the Result structure the standard PrintResult routine is called
      PrintResult(Result,1);
   else
      options.Display='iter';
      options.Jacobian='on';
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k] =lsqcurvefit...
          ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);
   end
   
   Name='Population problem';
   fprintf('\n');   
   fprintf('\n\n-----------------------------------------------------\n\n');   
   fprintf('Test of problem %s',Name);   
   fprintf('\n');   
   fprintf('\n\n-----------------------------------------------------\n\n');   
   fprintf('\n');   

   t=[-15; -10; -5; 0; 5; 10; 15];
   y=[60; 64; 71; 80; 90; 101; 116];
   x_0=[80,1]';
   x_L=[-Inf,-Inf]';
   x_U=[Inf,Inf]';
   
   % Define input needed by lsq_rJ
   % Note that this structure is just additional user input
   % It is not recognized by Tomlab as a problem structure, because
   % the field N is missing
   InProb.P=3;      % Set the problem number   
   InProb.Name=Name;% Set the problem name, nicer output in PrintResult later on
   % We send the t and y vector down to lsq_rJ:
   InProb.LS.t=t;
   InProb.LS.y=y;

   % Define Tomlab problem structure to input the solver to be used
   % If the call is to the Opt Tbx lsqnonlin, it is simply not used
   global otxProb
   otxProb = ProbDef;
   otxProb.Solver.Tomlab = 'nlssol';

   if exist('optimset')
      options=optimset;
   else
      options=[];
   end
   options.Diagnostics='on';
   options.Display='iter';

   if nargout('lsqcurvefit')==8
      % Tomlab lsqnonlin is returning the Result structure as well
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =lsqcurvefit...
         ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);

      PrintResult(Result,1);
   else
      options.Jacobian='on';
      [x, f_k, r_k, ExitFlag, Output, Lambda, J_k] =lsqcurvefit...
         ('curve_rJ', x_0, t, y, x_L, x_U, options, InProb);
   end

end
otxProb = []; % Clean global structure

