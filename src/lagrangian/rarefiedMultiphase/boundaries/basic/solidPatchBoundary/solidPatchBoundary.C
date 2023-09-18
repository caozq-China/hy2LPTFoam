/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "solidPatchBoundary.H"
// #include "IFstream.H"
#include "graph.H"
#include "solidParticleCouplingCloud.H"
#include "Time.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidPatchBoundary, 0);

defineRunTimeSelectionTable(solidPatchBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidPatchBoundary::solidPatchBoundary
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
:
    solidBoundaryBase(mesh, spc, dict, "patch"),
    faces_(),
    patchSurfaceArea_(Zero),
    totalPatchSurfaceArea_(Zero),
    cells_(),
    velocity_(Zero),
    density_(Zero),
    temperature_(Zero),
    measurePropertiesAtWall_(false),
    calculateHeatConduction_(false),
    preIE_(0.0),
    preIMom_(Zero)
    
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    if (isA<cyclicPolyPatch>(patch))
    {
        FatalErrorInFunction
            << "Patch: " << patchName_ << " is a cyclic boundary. It should be a patch." 
            << nl << "in: "
            << mesh_.time().system()/solidBoundaries::dictName
            << exit(FatalError);
    }

    for(label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
//         nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]); //area on one processor
    }
    
    totalPatchSurfaceArea_ = patchSurfaceArea_;

    if(Pstream::parRun())
    {
        reduce(totalPatchSurfaceArea_, sumOp<scalar>());  //total area on all processors
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<solidPatchBoundary> solidPatchBoundary::New
(
    const polyMesh& mesh,
    solidParticleCouplingCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.get<word>("boundaryModel")
    );

    Info<< "Selecting solidPatchBoundaryModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "solidPatchBoundary",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<solidPatchBoundary>
    (
        cstrIter()(mesh, spc, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidPatchBoundary::measurePropertiesBeforeControl(const solidParticleCoupling& p)
{

    if(measurePropertiesAtWall_)
    {
        const label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());
    
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

        const auto& constSolidProps = spc_.constSolidProps(p.typeID());
    
        const scalar m = constSolidProps.massSphere();
    
        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);
    
        scalar U_dot_nw = p.U() & nw;
    
        const vector Ut = p.U() - U_dot_nw*nw;
    
        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);
        
        solidBoundaryMeasurements& bm = spc_.boundaryFluxMeasurements();
        bm.rhoNBF()[p.typeID()][wppIndex][wppLocalFace] += invMagUnfA;

        bm.rhoMBF()[p.typeID()][wppIndex][wppLocalFace] += m*invMagUnfA;
        bm.linearKEBF()[p.typeID()][wppIndex][wppLocalFace] += 0.5*m*(p.U() & p.U())*invMagUnfA;
        bm.volumeFractionBF()[p.typeID()][wppIndex][wppLocalFace] += (pi/6.0)*pow(p.D(),3.0)*invMagUnfA;
        bm.momentumBF()[p.typeID()][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
    }
}

void solidPatchBoundary::measurePropertiesAfterControl
(
    const solidParticleCoupling& p, 
    const scalar heatOfConduction
)
{
    if(measurePropertiesAtWall_)
    {       
        const label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());
    
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);
        
//         const scalar deltaT = mesh_.time().deltaTValue();
    
        const solidParticleCoupling::constantProperties& constSolidProps(spc_.constSolidProps(p.typeID()));
    
        const scalar m = constSolidProps.massSphere();
    
        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);
    
        scalar U_dot_nw = p.U() & nw;
    
        vector Ut = p.U() - U_dot_nw*nw;
    
        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);
        
        solidBoundaryMeasurements& bm = spc_.boundaryFluxMeasurements();
        
        bm.rhoNBF()[p.typeID()][wppIndex][wppLocalFace] += invMagUnfA;

        bm.rhoMBF()[p.typeID()][wppIndex][wppLocalFace] += m*invMagUnfA;
        bm.linearKEBF()[p.typeID()][wppIndex][wppLocalFace] += 0.5*m*(p.U() & p.U())*invMagUnfA;
        bm.volumeFractionBF()[p.typeID()][wppIndex][wppLocalFace] += (pi/6.0)*pow(p.D(),3.0)*invMagUnfA;
        bm.momentumBF()[p.typeID()][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

        scalar RWF = spc_.dsmcCloudReference()->axiRWF(wpp.faceCentres()[wppLocalFace]);
        scalar nParticle = spc_.nSolidParticles()*RWF;

        bm.numberFluxBF()[p.typeID()][wppIndex][wppLocalFace] += nParticle/fA;
        bm.massFluxBF()[p.typeID()][wppIndex][wppLocalFace] += m*nParticle/fA;
        
        if(calculateHeatConduction_)
        {
            bm.conductiveHeatFluxBF()[p.typeID()][wppIndex][wppLocalFace] += heatOfConduction*nParticle/fA;
        }
    }
}

void solidPatchBoundary::specularReflection(solidParticleCoupling& p)
{
    vector& U = p.U();

    vector nw = p.normal();
    nw.normalise();

    const scalar U_dot_nw = U & nw;

    if (U_dot_nw > 0.0)
    {
        U -= 2.0*U_dot_nw*nw;
    }
}

const labelList& solidPatchBoundary::controlPatch() const
{
    return faces_;
}

const labelList& solidPatchBoundary::controlZone() const
{
    return cells_;
}

const scalar& solidPatchBoundary::density() const
{
    return density_;
}

scalar& solidPatchBoundary::density()
{
    return density_;
}

const vector& solidPatchBoundary::velocity() const
{
    return velocity_;
}

vector& solidPatchBoundary::velocity()
{
    return velocity_;
}

const scalar& solidPatchBoundary::temperature() const
{
    return temperature_;
}

scalar& solidPatchBoundary::temperature()
{
    return temperature_;
}

} // End namespace Foam

// ************************************************************************* //
