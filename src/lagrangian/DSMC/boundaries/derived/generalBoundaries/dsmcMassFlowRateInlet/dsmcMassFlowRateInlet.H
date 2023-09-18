/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Class
    Foam::dsmcMassFlowRateInlet

Description
    This class provides functions to insert the mass flow rate at a patch,
    given a temperature

    From Farbar and Boyd's publication referenced below, equation 8

        \f[
            \dot{n}_{s,i} = \frac{Y_s}{U}\frac{\dot{m}}{A}\frac{Av}{M_s}
        \f]

    where
    \vartable
        n_{s,i}     | Species inlet face number density [#/m^3]
        Y_s         | Species mass fraction
        U           | Inlet velocity [m/s]
        \dot{m}     | Mass flow rate [kg/s]
        A           | Inlet area [m^2]
        Av          | Avogadro's number [#/mol]
        M_s         | Species molecular weight [kg]
    \endvartable

    Use in the boundariesDict file, in conjunction with a deletion patch.

    Reference:
        E. Farbar and I.D. Boyd,  "Subsonic flow boundary conditions for the direct simulation Monte Carlo method", Computers & Fluids, 102, 99-110, 2014.

Usage
    \table
        Property         | Description                        | Required  | Default
        typeIds          | List of species                    | yes       |
        inletTemperature | Inlet temperature [K]              | yes       |
        massFlowRate     | Mass flow rate [kg/s]              | yes       |
        initialVelocity  | Guess for the inlet velocity [m/s] | yes       |
        moleFractions    | Mole fractions                     | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    generalBoundaryProperties
    {
        patchName                           flow;
    }

    boundaryModel   dsmcMassFlowRateInletPatch;

    dsmcMassFlowRateInletPatchProperties
    {
        typeIds                     (Ar);
        inletTemperature            300;
        massFlowRate                1e-11;
        initialVelocity             (1 0 0);
        moleFractions
        {
            Ar          1.0;
        }
    }
    \endverbatim

SourceFiles
    dsmcMassFlowRateInlet.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcMassFlowRateInlet_H
#define dsmcMassFlowRateInlet_H

#include "dsmcGeneralBoundary.H"
#include "polyMesh.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class dsmcMassFlowRateInlet Declaration
\*---------------------------------------------------------------------------*/

class dsmcMassFlowRateInlet
:
    public dsmcGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        scalarField moleFractions_;
        vectorField inletVelocity_;
        vectorField previousInletVelocity_;
        vectorField momentum_;
        List<scalarField> n_;
        scalarField mass_;
        scalar massFlowRate_;
        scalar inletTemperature_;
        vector initialVelocity_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("dsmcMassFlowRateInlet");

    //- Construct from components
    dsmcMassFlowRateInlet
    (
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~dsmcMassFlowRateInlet() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
