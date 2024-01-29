/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

void write(word matrixName, label timeIndex,
	   fvVectorMatrix M,
	   volVectorField& U)
{ 
    std::stringstream fileName;
    char timeStep[8];

    sprintf(timeStep, "%05d", timeIndex);

    fileName << matrixName << timeStep << ".m";

    std::fstream file(fileName.str().c_str(), std::ios::out);
    file.precision(15);

    const scalarField& diag = M.diag();
    const vectorField& source = M.source();
    const labelUList& owner = M.lduAddr().lowerAddr();
    const labelUList& neighbor = M.lduAddr().upperAddr();
    const scalarField& lower = M.lower();
    const scalarField& upper = M.upper();


    file << "Data = [" << std::endl;
    for(label i=0; i<diag.size(); i++)
	{
	    file << i+1 << "\t" << i+1 << "\t" << diag[i] << std::endl;
	}
    for(label f=0; f<upper.size(); f++)
	{
	    file << owner[f]+1 << "\t" << neighbor[f]+1 << "\t" << lower[f] << std::endl;
	    file << neighbor[f]+1 << "\t" << owner[f]+1 << "\t" << upper[f] << std::endl;
	}
    file << "];\n " << matrixName << "Matrix = sparse(Data(:,1), Data(:,2), Data(:,3));\n";
    file << "Source = [" << std::endl;
    for(label i=0; i<source.size(); i++)
	{
	    file << i+1 << "\t" << source[i][0] << "\t" << source[i][1] << "\t" << source[i][2] << std::endl;
	}
    file << "];\n";
    file << "Solution = [" << std::endl;
    for(label i=0; i<U.size(); i++)
	{
	    file << i+1 << "\t" << U[i][0] << "\t" << U[i][1] << "\t" << U[i][2] << std::endl;
	}
    file << "];\n";
    file.close();











    //
    //
    // Output Boundary condition
    //
    //
    for (int k = 0; k < 3; k++) {
	char name[1000];
	char surfix[3][100] = {
	    "x", "y", "z"
	};
	sprintf(name, "Bdry%d_%s.m", timeIndex, surfix[k]);
	FILE *fp = fopen(name, "w");
	fprintf(fp, "bc%s = [\n", surfix[k]);
	    
	forAll(M.internalCoeffs(), patchI)
	    {
		const label*  fPtr = M.lduAddr().patchAddr(patchI).begin();
		Field<double> intF(M.internalCoeffs()[patchI].component(k));
		Field<double> bouF(M.boundaryCoeffs()[patchI].component(k));

		forAll(M.lduAddr().patchAddr(patchI), faceI)
		    {
			label fCell = fPtr[faceI];
			//diag(fCell) -= intF[faceI];
			fprintf(fp, "%5d %30.15e %30.15e\n",
				fCell+1, intF[faceI], bouF[faceI]
				);
		    }
	    }

	fprintf(fp, "];\n");
	fclose(fp);
    }

    
    return;

}	


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    //turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Generating random initial values
    //std::vector<double> F_arr = {};
    /*scalar LO=0.1;
    scalar HI=0.9;
    for (int ii = 0; ii < x.size(); ii++)
    {
        //F_arr.push_back(LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)))); 
        x[ii] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
    };*/

    //Info<< "Random Values: " << x[0] << " " << x[1] << " " << x[2] << " " << x[2000] << " " << x[2333] << endl;
    //alpha1=1;
    Info<< "\nStarting time loop\n" << endl;


    int iterr=0;

    //simple.loop();
    while (simple.loop())
    //while (!converged && iterr<1000)
    //while (true)
    {
        //Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "Time = " << iterr << nl << endl;
 	    iterr += 1;

	    vf1.correctBoundaryConditions();
	    vf2.correctBoundaryConditions();
        vf1Pen = vf1*(1 + q)/(vf1 + q);
        vf2Pen = vf2*(1 + q)/(vf2 + q);


        // --- Pressure-velocity SIMPLE corrector
        {
            #include "U1Eqn.H"
            #include "p1Eqn.H"
        }

        // --- Pressure-velocity2 SIMPLE corrector
        {
            #include "U2Eqn.H"
            #include "p2Eqn.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);

    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
