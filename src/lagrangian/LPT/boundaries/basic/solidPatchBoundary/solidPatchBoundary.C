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
#include "graph.H"
#include "solidParcelCloud.H"
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
    solidParcelCloud& spc,
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
    preIkineticE_(0.0),
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
    solidParcelCloud& spc,
    const dictionary& dict
)
{
    word modelType
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting solidPatchBoundaryModel "
         << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalError
            << "solidPatchBoundary::New(const dictionary&) : " << endl
            << "    unknown solidPatchBoundary type "
            << modelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<solidPatchBoundary>
    (
        cstrIter()(mesh, spc, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidPatchBoundary::measurePropertiesBeforeControl(const solidParcel& p)
{
    const label wppIndex = patchId_;
    const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
    const label wppLocalFace = wpp.whichFace(p.face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar m = p.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = p.U() & nw;

    const vector Ut = p.U() - U_dot_nw*nw;

    scalar invMagUnfA = p.RWF()/max(mag(U_dot_nw)*fA, VSMALL);

    boundaryMeasurements& bm = spc_.boundaryFluxMeasurements();
    bm.rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfA;

    bm.rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*invMagUnfA;
    bm.momentumBF()[p.typeId()][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    //- pre-interaction energy
    preIkineticE_ = 0.5*m*(p.U() & p.U());

    //- pre-interaction momentum
    preIMom_ = m*p.U();
}

void solidPatchBoundary::measurePropertiesAfterControl(const solidParcel& p)
{
    const label wppIndex = patchId_;
    const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
    const label wppLocalFace = wpp.whichFace(p.face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = mesh_.time().deltaTValue();

    const scalar m = p.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = p.U() & nw;

    vector Ut = p.U() - U_dot_nw*nw;

    scalar invMagUnfA = p.RWF()/max(mag(U_dot_nw)*fA, VSMALL);

    boundaryMeasurements& bm = spc_.boundaryFluxMeasurements();

    bm.rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfA;
    bm.rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*invMagUnfA;
    bm.momentumBF()[p.typeId()][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    //- post-interaction energy
    scalar postIkineticE = 0.5*m*(p.U() & p.U());

    //- post-interaction momentum
    const vector postIMom = m*p.U();

    // scalar nParticle = spc_.nParticle();
    // Note: Do _not_ use the cloud._nParticles() command here because it
    // assumes the wrong RWF. Because this calculation happens during a
    // parcel move step we have to use the parcels RWF here. Since the
    // measurements before/after control are always done within one move
    // step the parcels RWF does not change, therefore we can use it to
    // calculate the pre- and post-interaction values below.
    // Note: this is not the case for parcels that are stuck to the wall,
    // but in that case a special function called
    //   measurePropertiesAfterDesorption
    // is used to calculate the post interaction values. In that case the
    // pre-interaction values (preIE and preIMom) are _not_ used.
    //const scalar nParticle = p.RWF() * spc_.coordSystem().dtModel().nParticles(wppIndex, wppLocalFace);
    scalar nParticle = spc_.nParticle();
    if(spc_.axisymmetric())
    {
        const vector fC = wpp.faceCentres()[wppLocalFace];
        
        scalar radius = fC.y();
        
//             scalar radius = sqrt((p.position().y()*p.position().y()) + (p.position().z()*p.position().z()));
        
        scalar RWF = 1.0;

        RWF = 1.0 + spc_.maxRWF()*(radius/spc_.radialExtent());
      
        nParticle *= RWF;
    }
    
    // Info<<"preIkineticE_="<<preIkineticE_<<endl;
    // Info<<"postIkineticE="<<postIkineticE<<endl;
    // Info<<"deltaT*fA="<<deltaT*fA<<endl;
    const scalar deltaKineticE = nParticle*(preIkineticE_ - postIkineticE)/(deltaT*fA);
    const vector deltaMomentumFlux = nParticle*(preIMom_ - postIMom)/(deltaT*fA);

    bm.fDBF()[p.typeId()][wppIndex][wppLocalFace] += deltaMomentumFlux;
    bm.qBF()[p.typeId()][wppIndex][wppLocalFace] += deltaKineticE;
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
