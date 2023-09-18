/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "Field.H"
#include "fvcGrad.H"
#include "polyMeshTetDecomposition.H"

#include "BasicMethod.H"
#include "DualMethod.H"
// #include "Moment.H"

// Scalar interpolation
// defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethod<Foam::scalar>, 0);
// namespace Foam
// {
//     defineTemplateRunTimeSelectionTable
//     (
//         mppicAveragingMethod<Foam::scalar>,
//         dictionary
//     );
// }
// 
// // Vector interpolation
// defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethod<Foam::vector>, 0);
// namespace Foam
// {
//     defineTemplateRunTimeSelectionTable
//     (
//         Foam::mppicAveragingMethod<Foam::vector>,
//         dictionary
//     );
// }
// 
// // Basic interpolation
// defineNamedTemplateTypeNameAndDebug
// (
//     Foam::mppicAveragingMethod<Foam::scalar>::BasicMethod<Foam::scalar>,
//     0
// );
// Foam::mppicAveragingMethod<Foam::scalar>::
// adddictionaryConstructorToTable<Foam::mppicAveragingMethod<Foam::scalar>::BasicMethod<Foam::scalar>>
//     addBasicMethodscalarConstructorToTable_;
// 
// defineNamedTemplateTypeNameAndDebug
// (
//     Foam::mppicAveragingMethod<Foam::vector>::BasicMethod<Foam::vector>,
//     0
// );
// Foam::mppicAveragingMethod<Foam::vector>::
// adddictionaryConstructorToTable<Foam::mppicAveragingMethod<Foam::vector>::BasicMethod<Foam::vector>>
//     addBasicMethodvectorConstructorToTable_;
// 
// // Dual interpolation
// defineNamedTemplateTypeNameAndDebug
// (
//     Foam::mppicAveragingMethod<Foam::scalar>::DualMethod<Foam::scalar>,
//     0
// );
// Foam::mppicAveragingMethod<Foam::scalar>::
// adddictionaryConstructorToTable<Foam::mppicAveragingMethod<Foam::scalar>::DualMethod<Foam::scalar>>
//     addDualMethodscalarConstructorToTable_;
// 
// defineNamedTemplateTypeNameAndDebug
// (
//     Foam::mppicAveragingMethod<Foam::vector>::DualMethod<Foam::vector>,
//     0
// );
// Foam::mppicAveragingMethod<Foam::vector>::
// adddictionaryConstructorToTable<Foam::mppicAveragingMethod<Foam::vector>::DualMethod<Foam::vector>>
//     addDualMethodvectorConstructorToTable_;



namespace Foam
{
    // Scalar interpolation
    defineNamedTemplateTypeNameAndDebug(mppicAveragingMethod<scalar>, 0);
    defineTemplateRunTimeSelectionTable
    (
        mppicAveragingMethod<scalar>,
        dictionary
    );

    // Vector interpolation
    defineNamedTemplateTypeNameAndDebug(mppicAveragingMethod<vector>, 0);
    defineTemplateRunTimeSelectionTable
    (
        mppicAveragingMethod<vector>,
        dictionary
    );

    namespace mppicAveragingMethods
    {
        // Basic interpolation
        defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethods::BasicMethod<scalar>, 0);
        mppicAveragingMethod<scalar>::
            adddictionaryConstructorToTable<BasicMethod<scalar> >
            addBasicMethodscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethods::BasicMethod<vector>, 0);
        mppicAveragingMethod<vector>::
            adddictionaryConstructorToTable<BasicMethod<vector> >
            addBasicMethodvectorConstructorToTable_;

        // Dual interpolation
        defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethods::DualMethod<scalar>, 0);
        mppicAveragingMethod<scalar>::
            adddictionaryConstructorToTable<DualMethod<scalar> >
            addDualMethodscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Foam::mppicAveragingMethods::DualMethod<vector>, 0);
        mppicAveragingMethod<vector>::
            adddictionaryConstructorToTable<DualMethod<vector> >
            addDualMethodvectorConstructorToTable_;
// 
//         // Moment interpolation
//         defineNamedTemplateTypeNameAndDebug(Moment<scalar>, 0);
//         AveragingMethod<scalar>::
//             adddictionaryConstructorToTable<Moment<scalar> >
//             addMomentscalarConstructorToTable_;
// 
//         defineNamedTemplateTypeNameAndDebug(Moment<vector>, 0);
//         AveragingMethod<vector>::
//             adddictionaryConstructorToTable<Moment<vector> >
//             addMomentvectorConstructorToTable_;
    }
}


// ************************************************************************* //
