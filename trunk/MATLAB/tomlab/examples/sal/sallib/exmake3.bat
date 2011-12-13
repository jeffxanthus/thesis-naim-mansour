@echo off

call mcc -B csharedlib:libcplexqp checkType.m cplexmex.dll cplexTL.m cpx2cbinfo.m cpx2retvec.m defblbu.m endSolve.m iniSolve.m LineParamDef.m LineParamSet.m mkbound.m optParamDef.m optParamSet.m PrintResult.m ProbCheck.m ProbDef.m qpAssign.m ResultDef.m snoptTL.m SOLGet.m tomFiles.m tomlablic.dll tomlabVersion.m tomRun.m xnargin.m xnargout.m solveqp.m DefPar.m qp_H.m qp_f.m qp_g.m cplex.m

if %ERRORLEVEL% geq 1 goto mccerror

call mbuild -v application.c libcplexqp.lib

if %ERRORLEVEL% geq 1 goto mbuilderror

exit /b 0

:mccerror

echo Failed to build shared library. Make sure that you have run the copyexfiles.m script. See the README file for details.
exit /b %ERRORLEVEL%

:mbuilderror

echo Failed to build client application. Make sure your dynamic library path is properly set. See the README file for details.
exit /b %ERRORLEVEL%