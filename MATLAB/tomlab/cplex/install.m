% ============================================================
% Installation of TOMLAB /CPLEX
% ============================================================
%
% *** Windows ***
%
% TOMLAB /CPLEX is installed using the same InstallShield
% self-extracting EXE file as the rest of TOMLAB.
%
% TOMLAB /CPLEX is installed in the following directories:
%
% tomlab\cplex - Main routines, M-files and MEX files
% tomlab\cplex\network  - Network interface
% tomlab\cplex\examples - A collection of example files
% tomlab\common - Common files
%
% In addition, if the user has selected the ILOG Shared Files in the
% installer dialog, the cplex*.dll file is installed in
% \tomlab\shared.
%
% Usually, the installer modifies the Windows PATH variable to include
% \tomlab\shared.
%
% If this has not happened after running the installer,
% please manually add the full path (e.g. c:\tomlab\shared)
% to the Windows PATH variable.
%
% --------------------------------------------------------------
%
% *** Linux ***
%
% TOMLAB /CPLEX is installed when extracting the tomlab-lnx-setup.tar.gz file.
%
% The ILOG CPLEX libcplex*.so library is included in
%
% ==============================================================
% Using TOMLAB /CPLEX with an existing ILOG CPLEX installation
% ==============================================================
%
% If you already have a correct installation of, and a valid license
% for ILOG CPLEX, you will only need to install the TOMLAB license
% file which has been sent to you as email attachment.
%
% The cplex*.dll (or libcplex*.so, Linux) must be present in the Windows
% PATH variable.
%
% ==============================================================
% Starting and running TOMLAB /CPLEX
% ==============================================================
%
% To start TOMLAB /CPLEX with TOMLAB Base Module, the user should
% execute only the \tomlab\startup.m file.
%
% >> cd d:/tomlab
% >> startup
%
% HINT: You may also set the TOMLAB and/or TOMLAB /CPLEX paths
% permanently in the Matlab system. To find out which PATHs are used,
% run the desired startup command as described above, then use the
% PATH command to see what paths TOMLAB created and set these
% permanently on your system. See the Matlab documentation on how to
% set Matlab paths.
%
% ==============================================================
% Usage hints
% ==============================================================
%
% Create a working directory for all your TOMLAB /CPLEX runs, e.g.
%
%    mkdir('d:\cprun')
%    cd('d:\cprun')
%
% -------------------------------------------------------------------------
% To check the TOMLAB /CPLEX installation, try out the example files, e.g.
%
%    x=cpxTest1;
%
% or
%
%    cpxKnaps
%
% Note that cpxKnaps takes two input argument that you can play around with.
%
% -------------------------------------------------------------------------
% If TOMLAB Base Module is installed, to run a first test using the
% test examples in the TOMLAB distribution, in \tomlab\testprob, run
%
%    cpxtomtest1
%
% -------------------------------------------------------------------------
% To see how to quickly define a problem in the TOMLAB format using the
% TQ format (TOMLAB Format), run
%
%    cpxtomtest2
%
% -------------------------------------------------------------------------
% To see how to set parameters that influence the runs in TOMLAB /CPLEX, run
%
%    cpxKnapsTL
%
% This file does the same as cpxKnaps, but is using the TQ format, and
% setting the relevant parameters into the Prob structure before calling the
% TOMLAB driver routine tomRun.
%
% -------------------------------------------------------------------------
% Contact sales@tomopt.com for license inquiries
%
% -------------------------------------------------------------------------

help install.m