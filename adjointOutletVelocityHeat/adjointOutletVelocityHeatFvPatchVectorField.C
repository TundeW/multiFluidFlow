/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityHeatFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    a_()
{}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
     a_( ptf.a_ )
{}


Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityHeatFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi1b");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U1");

    const fvPatchField<vector>& Uap =
    	patch().lookupPatchField<volVectorField, vector>("U1b");

    const fvsPatchField<scalar>& phip =
    	patch().lookupPatchField<surfaceScalarField, scalar>("phi1");
    	
    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
     //dimensionedScalar nu(transportProperties.lookup("nu1"));
    dimensionedScalar nu
    (
        "nu1",
        dimViscosity,
        transportProperties
    );

    const fvsPatchField<scalar>& phi2ap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi2b");

    const fvPatchField<vector>& U2p =
        patch().lookupPatchField<volVectorField, vector>("U2");

    const fvPatchField<vector>& U2ap =
    	patch().lookupPatchField<volVectorField, vector>("U2b");

    const fvsPatchField<scalar>& phi2p =
    	patch().lookupPatchField<surfaceScalarField, scalar>("phi2");
    	
     //dimensionedScalar nu_2(transportProperties.lookup("nu2"));

    dimensionedScalar nu_2
    (
        "nu2",
        dimViscosity,
        transportProperties
    );
   //const incompressible::RASModel& rasModel =
   // 	db().lookupObject<incompressible::RASModel>("RASProperties");

   // scalarField nueff = rasModel.nuEff()().boundaryField()[patch().index()];

    const scalarField& deltainv = 
    	patch().deltaCoeffs(); // dist^(-1) 

//Primal velocity, mag of normal component and tangential component
    scalarField Up_ns = phip/patch().magSf();

    vectorField Up_t = Up - (phip * patch().Sf())/(patch().magSf()*patch().magSf());

//Tangential component of adjoint velocity in neighbouring node
    vectorField Uaneigh = Uap.patchInternalField();
    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;

    vectorField Uap_t = (nu.value()*deltainv*Uaneigh_t) / (Up_ns+nu.value()*deltainv) ;

    vectorField Uap_n = (phiap * patch().Sf())/(patch().magSf()*patch().magSf());

 if (a_!=1) 
    {

        Up_ns = phi2p/patch().magSf();

        Up_t = U2p - (phi2p * patch().Sf())/(patch().magSf()*patch().magSf());

    //Tangential component of adjoint velocity in neighbouring node
        vectorField Uaneigh = U2ap.patchInternalField();
        vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();
        vectorField Uaneigh_t = Uaneigh - Uaneigh_n;

        Uap_t = (nu_2.value()*deltainv*Uaneigh_t) / (Up_ns+nu_2.value()*deltainv) ;

        Uap_n = (phi2ap * patch().Sf())/(patch().magSf()*patch().magSf());
    }   
    operator==(Uap_t+Uap_n);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityHeatFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityHeatFvPatchVectorField
    );
}


// ************************************************************************* //
