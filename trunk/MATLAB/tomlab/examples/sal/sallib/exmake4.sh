#!/bin/sh

mcc -B csharedlib:libcplexqp checkType.m cplexmex.mexglx cplexTL.m cpx2cbinfo.m cpx2retvec.m defblbu.m endSolve.m iniSolve.m LineParamDef.m LineParamSet.m mkbound.m optParamDef.m optParamSet.m PrintResult.m ProbCheck.m ProbDef.m qpAssign.m ResultDef.m snoptTL.m SOLGet.m tomFiles.m tomlablic.mexglx tomlabVersion.m tomRun.m xnargin.m xnargout.m solveqp.m DefPar.m qp_H.m qp_f.m qp_g.m cplex.m

if [ "$?" -eq "0" ]
then
    mbuild application.c -L. -lcplexqp -I. -DMCC4
    if [ "$?" -ne "0" ]
    then
	echo Failed to build client application. Make sure your dynamic library path is properly set. See the README file for details.
    fi
else
    echo Failed to build shared library. Make sure that you have run the copyexfiles.m script. See the README file for details.
fi
