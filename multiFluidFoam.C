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

void write_velocity_distribution(std::ofstream &labelFile, volVectorField& U1){

    // Send column names to the stream
    for(label i=0; i<U1.size(); i++)
	{
        double v0 = U1[i][0];
        double v1 = U1[i][1];
        double v2 = U1[i][2];
        //double U1_mag = pow(pow(v0, 2.0) + pow(v1, 2.0) + pow(v2, 2.0), 0.5);
        //double U1_mag = pow(pow(v0, 2.0) + pow(v1, 2.0) + pow(v2, 2.0), 0.5);
        double U1_mag = Foam::sqrt(Foam::pow(v0, 2.0) + Foam::pow(v1, 2.0) + Foam::pow(v2, 2.0));

	    labelFile << U1_mag;
        if(i != U1.size() - 1) {
            labelFile << ","; // No comma at end of line
        }
	}
    labelFile << "\n";
    
    
}

void write_density_distribution(std::ofstream &labelFile, volScalarField& vf){

    // Send column names to the stream
    for(label i=0; i<vf.size(); i++)
	{
	    labelFile << vf[i];
        if(i != vf.size() - 1) {
            labelFile << ","; // No comma at end of line
        }
	}
    labelFile << "\n";
    
    
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
    std::string velMagfilename = "velocity_magnitude_data.csv";
    std::string density1filename = "density1_data.csv";
    std::string density2filename = "density2_data.csv";
    // Create an output filestream object for both input and label data files
    std::ofstream inputFile(filename);
    std::ofstream labelFile(labelfilename);
    std::ofstream velMagFile(velMagfilename);
    std::ofstream density1File(density1filename);
    std::ofstream density2File(density2filename);

    int dataSize = 1000; // Specify the size of the dataset
    srand (time(0)); // Seed random number generator with system time.
    
    for (int out_iter = 0; out_iter < dataSize; out_iter++) {

        #include "addCheckCaseOptions.H"
        #include "setRootCaseLists.H"
        #include "createTime.H"
        #include "createMesh.H"
        #include "createControl.H"
        #include "createFields.H"
        #include "initContinuityErrs.H"
        Info << "\nStarting time loop\n" << endl;
    
        
        double x_init[]={0.2225105,0.0238811,0.159695,0.0202648,0.118059,0.0742291,0.125606,0.0275261,0.220039,0.0220306,0.0914213,0.0419477,0.00560588,0.00851001,0.0140996,0.0112353,0.00490028,0.00966353,0.0122237,0.0117558,0.0116712,0.00967284,0.0150398,0.0106401};
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
            
            x_init[ii] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
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
        forAll(Cell, cellI)
        {
            #include "densityProjection.H"
        }
        vf1.write();
        vf2.write();

        write_density_distribution(density1File, vf1);
        write_density_distribution(density2File, vf2);
            
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

            runTime.write();

            runTime.printExecutionTime(Info);

        } while (simple.loop());

        write_output(labelFile, U2);
        write_velocity_distribution(velMagFile, U2);
        Info<< "End\n" << endl;
    
    }

    // Close the file
    inputFile.close();
    labelFile.close();
    velMagFile.close();
    density1File.close();
    density2File.close();

    return 0;
}


// ************************************************************************* //
