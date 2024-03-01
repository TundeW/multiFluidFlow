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

template <typename T, size_t N>
double array_size(const T (&array)[N]) {
    return static_cast<double>(N);
}

void write_input(std::ofstream &inputFile, double* x_init, double x_init_size){
  
    // Send column names to the stream
    for(int j = 0; j < x_init_size; ++j)
    {
        inputFile << x_init[j];
        if(j != x_init_size - 1) {
            inputFile << ","; // No comma at end of line
        }
    }
    inputFile << "\n";
    
    
}

void write_output(std::ofstream &labelFile, volVectorField& U1){

    // Send column names to the stream
    for(label i=0; i<U1.size(); i++)
	{
	    labelFile << U1[i][0] << "," << U1[i][1] << "," << U1[i][2];
        if(i != U1.size() - 1) {
            labelFile << ","; // No comma at end of line
        }
	}
    labelFile << "\n";
    
    
}

void write_velocity_distribution(std::ofstream &label2File, volVectorField& U1, volVectorField& U2){

    // Send column names to the stream
    for(label i=0; i<U1.size(); i++)
	{
	    label2File << U1[i][0]+U2[i][0] << "," << U1[i][1]+U2[i][1] << "," << U1[i][2]+U2[i][2];
        /*if(i != U1.size() - 1) {
            label2File << ","; // No comma at end of line
        }*/
        label2File << "\n";
	}
    label2File << "\n";  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    std::string filename = "input_data.csv";
    std::string labelfilename = "label_data.csv";
    std::string label2filename = "velocity_data.csv";
    // Create an output filestream object for both input and label data files
    std::ofstream inputFile(filename);
    std::ofstream labelFile(labelfilename);
    std::ofstream label2File(label2filename);

    int dataSize = 1; // Specify the size of the dataset
    srand (time(0)); // Seed random number generator with system time.
    
    for (int out_iter = 0; out_iter < dataSize; out_iter++) {

        #include "addCheckCaseOptions.H"
        #include "setRootCaseLists.H"
        #include "createTime.H"
        #include "createMesh.H"
        #include "createControl.H"
        #include "createFields.H"
        #include "initContinuityErrs.H"
        #include "CourantNo.H"
        Info << "\nStarting time loop\n" << endl;
    
        
        //double x_init[]={0.2225105,0.0238811,0.159695,0.0202648,0.118059,0.0742291,0.125606,0.0275261,0.220039,0.0220306,0.0914213,0.0419477,0.00560588,0.00851001,0.0140996,0.0112353,0.00490028,0.00966353,0.0122237,0.0117558,0.0116712,0.00967284,0.0150398,0.0106401};
        double x_init[]={0.025,0.025,0.175,0.025,0.225,0.075,0.125,0.075,0.075,0.075,0.075,0.025,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
        //x_init[0] = out_iter * 0.1;
        double x_init_size = array_size(x_init);

        //Update values of x with random values within a specified limit
        for (int ii = 0; ii < x_init_size; ii++)
        {   
            scalar LO=0;
            scalar HI;
            if ((ii%2 == 0) & (ii < x_init_size/2) ) {
                LO = 0.02;
                HI=0.23; // Pipe Locations x-coordinates
            }
            else if ((ii%2 != 0) & (ii < x_init_size/2) ) {
                LO = 0.02;
                HI=0.08; // Pipe Locations x-coordinates
            }
            else if ( (ii >= x_init_size/2) & (ii < x_init_size*0.75) ) {
                int ind = 2*(ii-x_init_size/2) + 1;
                LO = 0;
                HI = std::min(0.1 - x_init[ind], x_init[ind]); //Pipe Inner Radius (y-coord + inner radius < 0.1)
            }
            else if ( ii >= x_init_size*0.75 ) {
                int ind = 2*(ii-x_init_size*0.75) + 1;
                int ind2 = (x_init_size/2) + ii - (x_init_size*0.75);
                LO = 0;
                HI = std::min(0.1 - x_init[ind] - x_init[ind2], x_init[ind] - x_init[ind2]); //Pipe Thickness (y-coord + inner radius + thickness < 0.1)
            }
            
            //x_init[ii] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
            //Info << "out_iter: " << out_iter << " X_init[" << ii << "]: " << x_init[ii] << " LO: " << LO << " HI: " << HI << endl;
        };

        write_input(inputFile, x_init, x_init_size);

        int iterr=0;

        //simple.loop();
        //while (simple.loop())
        //while (!converged && iterr<1000)
        //while (true)
        simple.loop();
        Info<< "Recalculating vf1 and vf2\n" << endl;
        const volVectorField& Cell = mesh.C();
        int iiter = 0;
        forAll(Cell, cellI)
        {
            if (iiter < 60)
            {
                Info << "Element " << iiter << ": " << Cell[cellI] << endl;
            }
            iiter = iiter + 1;
            #include "densityProjection.H"
        }
        vf1.write();
        vf2.write();

        double qq = 0.01;
        for(label i=0; i<U1.size(); i++)
	    {   
            scalar vf1I = vf1[i];
            scalar vf2I = vf2[i];
            double vf1i = vf1I;
            double vf2i = vf2I;
            double vf1PenD = vf1i*(1 + qq)/(vf1i + qq);
            double vf2PenD = vf2i*(1 + qq)/(vf2i + qq);
	        U[i][0] = vf1PenD*(vf2PenD*(1) + (1 - vf2PenD)*1e-3) + (1 - vf1PenD)*1e-3; 
            U[i][1] = 1e-3; 
            U[i][2] = vf1PenD*(vf2PenD*1e-3 + (1 - vf2PenD)*(-1)) + (1 - vf1PenD)*1e-3;    
	    }
        U.write();
        phi = fvc::flux(U);

        do
        {
            //Info<< "Time = " << runTime.timeName() << nl << endl;
            Info<< "Time = " << iterr << nl << endl;
            iterr += 1;

            vf1.correctBoundaryConditions();
            vf2.correctBoundaryConditions();
            vf1Pen = vf1*(1 + q)/(vf1 + q);
            vf2Pen = vf2*(1 + q)/(vf2 + q);
            Info<< "Projection Completed!" << endl;


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

            dimensionedScalar Uzero("Uzero", dimensionSet(0,1,-1,0,0,0,0),0.0);
            dimensionedScalar Uone("Uone", dimensionSet(0,1,-1,0,0,0,0),1.0);  

            DT=vf1Pen*(vf2Pen*DT2 + (1 - vf2Pen)*DT1) + (1 - vf1Pen)*DTS;
	        DT.correctBoundaryConditions();

            fvScalarMatrix TEqn
            (
                fvm::div(phi1+phi2, T)
                //fvm::div(fvc::flux(U1)+fvc::flux(U2), T)
                - fvm::laplacian(DT, T)
            );

	        TEqn.relax();
            TEqn.solve();
	        T.correctBoundaryConditions();    

            runTime.write();

            runTime.printExecutionTime(Info);

        } while (simple.loop());

        //double qq = 0.01;
        /*for(label i=0; i<U1.size(); i++)
	    {   
            scalar vf1I = vf1[i];
            scalar vf2I = vf2[i];
            double vf1i = vf1I;
            double vf2i = vf2I;
            double vf1PenD = vf1i*(1 + qq)/(vf1i + qq);
            double vf2PenD = vf2i*(1 + qq)/(vf2i + qq);
	        //U1[i][0] = vf1PenD*(vf2PenD*(1) + (1 - vf2PenD)*0) + (1 - vf1PenD)*0; 
            U1[i][0] = 0;
            U1[i][1] = 0; 
            U1[i][2] = vf1PenD*(vf2PenD*0 + (1 - vf2PenD)*(-1)) + (1 - vf1PenD)*0;
            //U1[i][2] = 0;
            U2[i][0] = vf1PenD*(vf2PenD*(1) + (1 - vf2PenD)*0) + (1 - vf1PenD)*0; 
            //U2[i][0] = 0;
            U2[i][1] = 0; 
            U2[i][2] = 0;
            //U2[i][2] = vf1PenD*(vf2PenD*0 + (1 - vf2PenD)*(-1)) + (1 - vf1PenD)*0;
    
	    }
        U1.correctBoundaryConditions();
        U2.correctBoundaryConditions();

        phi1 = fvc::flux(U1);
        phi2 = fvc::flux(U2);

        runTime.write();

        runTime.printExecutionTime(Info);

        DT=vf1Pen*(vf2Pen*DT2 + (1 - vf2Pen)*DT1) + (1 - vf1Pen)*DTS;
	    DT.correctBoundaryConditions();

        fvScalarMatrix TEqn
        (
            fvm::div(phi1+phi2, T)
            //fvm::div(fvc::flux(U1)+fvc::flux(U2), T)
            - fvm::laplacian(DT, T)
        );

	    TEqn.relax();
        TEqn.solve();
	    T.correctBoundaryConditions();

        runTime.write();

        runTime.printExecutionTime(Info);*/

        write_output(labelFile, U1);
        write_velocity_distribution(label2File, U1, U2);
        Info<< "End\n" << endl;

        
    
    }

    // Close the file
    inputFile.close();
    labelFile.close();
    label2File.close();

    return 0;
}


// ************************************************************************* //
