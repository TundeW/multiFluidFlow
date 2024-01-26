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

#include "adjointOutletPressureHeatFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletPressureHeatFvPatchScalarField::
adjointOutletPressureHeatFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    a_()
{}


Foam::adjointOutletPressureHeatFvPatchScalarField::
adjointOutletPressureHeatFvPatchScalarField
(
    const adjointOutletPressureHeatFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
     a_( ptf.a_ )
{}


Foam::adjointOutletPressureHeatFvPatchScalarField::
adjointOutletPressureHeatFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::adjointOutletPressureHeatFvPatchScalarField::
adjointOutletPressureHeatFvPatchScalarField
(
    const adjointOutletPressureHeatFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletPressureHeatFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi1");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi1b");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("U1b");

    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
     //dimensionedScalar nu(transportProperties.lookup("nu1"));
    dimensionedScalar nu
    (
        "nu1",
        dimViscosity,
        transportProperties
    );

    const fvsPatchField<scalar>& phi2p =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi2");

    const fvsPatchField<scalar>& phi2ap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi2b");

    const fvPatchField<vector>& U2ap =
        patch().lookupPatchField<volVectorField, vector>("U2b");

     //dimensionedScalar nu_2(transportProperties.lookup("nu2"));
    dimensionedScalar nu_2
    (
        "nu2",
        dimViscosity,
        transportProperties
    );



     
    scalarField Up_n = phip / patch().magSf();//Primal

    scalarField Uap_n = phiap / patch().magSf();//Adjoint



   // const incompressible::hahaha& rasModel =
//	 db().lookupObject<incompressible::hahaha>("RASProperties");

   // scalarField nueff = rasModel.nuEff()().boundaryField()[patch().index()];

    const scalarField& deltainv = patch().deltaCoeffs(); // distance^(-1)

    scalarField Uaneigh_n = (Uap.patchInternalField() & patch().nf());

   operator ==((Up_n * Uap_n) +2*nu.value()*deltainv*(Uap_n-Uaneigh_n));
 
    if (a_!=1) 
    {

        Up_n = phi2p / patch().magSf();//Primal

        Uap_n = phi2ap / patch().magSf();//Adjoint

        Uaneigh_n = (U2ap.patchInternalField() & patch().nf());

        operator ==((Up_n * Uap_n) +2*nu_2.value()*deltainv*(Uap_n-Uaneigh_n));
    }   
    

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletPressureHeatFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletPressureHeatFvPatchScalarField
    );
}

// ************************************************************************* //
