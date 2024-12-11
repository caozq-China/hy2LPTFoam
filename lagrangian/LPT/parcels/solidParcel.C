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

#include "solidParcel.H"
#include "solidParcelCloud.H"
#include "meshTools.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::solidParcel::setCellValues
(
    trackingData& td,
    const scalar dt,
    const label celli
)
{
    tetIndices tetIs = currentTetIndices();
    Uc_ = td.UInterp().interpolate(position(), tetIs);

    rhoc_ = td.rhoInterp().interpolate(position(), tetIs);

    muc_ = td.muInterp().interpolate(position(), tetIs);

    Tc_ = td.TInterp().interpolate(position(), tetIs);

    if (Tc_ < td.cloud().constProps(typeId_).TMin())
    {
        Tc_ = td.cloud().constProps(typeId_).TMin();
    }

    kappac_ = td.kappaInterp().interpolate(position(), tetIs);

    Cpc_ = td.CpInterp().interpolate(position(), tetIs);

    Mac_ = td.MaInterp().interpolate(position(), tetIs);

    gammac_ = td.gammaInterp().interpolate(position(), tetIs);

    omegac_ = td.omegaInterp().interpolate(position(), tetIs);
}

void Foam::solidParcel::cellValueSourceCorrection
(
    trackingData& td,
    const scalar dt,
    const label celli
)
{
    tetIndices tetIs = currentTetIndices();

    Uc_ += td.cloud().momentumTrans()[celli]/massCell(celli);

    Tc_ = td.TInterp().interpolate(position(), tetIs);

    if (Tc_ < td.cloud().constProps(typeId_).TMin())
    {
        Tc_ = td.cloud().constProps(typeId_).TMin();
    }
}


void Foam::solidParcel::calcSurfaceValues
(
    trackingData& td,
    const scalar T,
    scalar& Ts,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappas,
    const label celli
) const
{
    // Surface temperature using two thirds rule
    Ts = (2.0*T + Tc_)/3.0;

    if (Ts < td.cloud().constProps(typeId_).TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting parcel surface temperature to "
                << td.cloud().constProps(typeId_).TMin() <<  nl << endl;
        }

        Ts = td.cloud().constProps(typeId_).TMin();
    }

    // Assuming thermo props vary linearly with T for small d(T)
    const scalar TRatio = Tc_/Ts;

    rhos = rhoc_*TRatio;

    tetIndices tetIs = currentTetIndices();
    mus = td.muInterp().interpolate(position(), tetIs)/TRatio;
    kappas = td.kappaInterp().interpolate(position(), tetIs)/TRatio;

    Pr = Cpc_*mus/kappas;
    Pr = max(ROOTVSMALL, Pr);
}

void Foam::solidParcel::calc
(
    trackingData& td,
    const scalar dt,
    const label celli
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = td.cloud().nParticle()*RWF();

    const scalar mass0 = mass();

    // Calc surface values
    // ~~~~~~~~~~~~~~~~~~~
    scalar Ts = 0.0;
    scalar rhos = 0.0;
    scalar mus = 0.0;
    scalar Pr = 0.0;
    scalar kappas = 0.0;

    // Reynolds number
    scalar Reynolds = Re(U_, d_, rhoc_, muc_);

    if(td.cloud().heatTransfer().active())
    {
        calcSurfaceValues(td, T_, Ts, rhos, mus, Pr, kappas, celli);
        Reynolds = Re(U_, d_, rhoc_, mus);
    }
    
    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;
    scalar dworksTrans = 0.0;

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Calculate new particle temperature
    T_ = calcHeatTransfer(td, dt, Reynolds, Pr, kappas, NCpW, Sh, celli, dhsTrans, Sph);

    // Calculate new particle velocity
    U_ = calcVelocity(td, dt, celli, Reynolds, muc_, mass0, Su, dUTrans, dworksTrans, Spu);

    U_ += UCorrect_;

    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Update momentum transfer
        td.cloud().momentumTrans()[celli] += np0*dUTrans;

        // Update momentum transfer coefficient
        td.cloud().momentumCoeff()[celli] += np0*Spu;

        // Update sensible enthalpy transfer
        td.cloud().qTrans()[celli] += np0*dhsTrans;

        // Update sensible enthalpy coefficient
        td.cloud().qCoeff()[celli] += np0*Sph;

        td.cloud().wTrans()[celli] += np0*dworksTrans;

}

Foam::scalar Foam::solidParcel::calcHeatTransfer
(
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW,
    const scalar Sh,
    const label celli,
    scalar& dhsTrans,
    scalar& Sph
)
{
    if (!td.cloud().heatTransfer().active())
    {
        return T_;
    }

    // Calc heat transfer coefficient
    scalar htc = td.cloud().heatTransfer().htc(d_, Re, Pr, kappa, NCpW, Mac());

    const scalar As = areaS(d_);

    // Calculate the integration coefficients
    const scalar bcp = htc*As/(mass()*Cp_);
    const scalar acp = bcp*Tc_;
    scalar ancp = Sh;
    
    ancp /= mass()*Cp_;
    // Integrate to find the new parcel temperature
    const scalar deltaT = td.cloud().TIntegrator().delta(T_, dt, acp + ancp, bcp);//Tc-T
    const scalar deltaTncp = ancp*dt;
    const scalar deltaTcp = deltaT - deltaTncp;

    // Calculate the new temperature and the enthalpy transfer terms
    scalar Tnew = T_ + deltaT;

    Tnew =  min
            (
                max
                (
                    Tnew,
                    td.cloud().constProps(typeId_).TMin()
                ),
                td.cloud().constProps(typeId_).TMax()
            );

    Sph = dt*htc*As;

    dhsTrans -= mass()*Cp_*deltaTcp;

    return Tnew;
}

const Foam::vector Foam::solidParcel::calcVelocity
(
    trackingData& td,
    const scalar dt,           // timestep
    const label celli,         // owner cell
    const scalar Re,           // Reynolds number
    const scalar mu,           // local carrier viscosity
    const scalar mass,         // mass
    const vector& Su,          // explicit particle momentum source
    vector& dUTrans,           // momentum transfer to carrier
    scalar& dworksTrans,
    scalar& Spu                // linearised drag coefficient
) //const
{
    // Momentum transfer coefficient
    const forceSuSp Fcp = td.cloud().forces().calcCoupled(*this, dt, mass, Re, mu);
    const forceSuSp Fncp = td.cloud().forces().calcNonCoupled(*this, dt, mass, Re, mu);
    const scalar massEff = td.cloud().forces().massEff(*this, mass);
    
    // New particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~

    // Shortcut splitting assuming no implicit non-coupled force ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*Uc_ + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;
    
    // Integrate to find the new parcel velocity
    const vector deltaU = td.cloud().UIntegrator().delta(U_, dt, acp + ancp, bcp);
    const vector deltaUncp = ancp*dt;
    const vector deltaUcp = deltaU - deltaUncp;

    // Calculate the new velocity and the momentum transfer terms
    vector Unew = U_ + deltaU;

    dUTrans -= massEff*deltaUcp;

    Spu = dt*Fcp.Sp();
    //- W=FS
    //- F= massEff*(Unew-p.U())/dt
    dworksTrans += (massEff*deltaUcp)&Unew;

    return Unew;
}

// * MOVE * //
bool Foam::solidParcel::move
(
    trackingData& td,
    const scalar trackTime
)
{
    solidParcel& p = static_cast<solidParcel&>(*this);
    td.switchProcessor = false;

    switch (td.part())
    {
        case trackingData::LinearTrack:
        {
            LinearMove(p, td, trackTime);
            break;
        }
        case trackingData::DampingNoTrack:
        {
            p.UCorrect() = Zero;
            p.UCorrect() = td.cloud().dampingModel().velocityCorrection(p, trackTime);
            td.keepParticle = true;
            td.switchProcessor = false;

            break;
        }
        case trackingData::PackingNoTrack:
        {
            p.UCorrect() = Zero;
            p.UCorrect() = td.cloud().packingModel().velocityCorrection(p, trackTime);
            td.keepParticle = true;
            td.switchProcessor = false;

            break;
        }
        case trackingData::CorrectTrack:
        {
            vector U = p.U();

            scalar f = p.stepFraction();

            p.U() = (1.0 - f)*p.UCorrect();

            LinearMove(p, td, trackTime);

            p.U() = U + (p.stepFraction() - f)*p.UCorrect();

            break;
        }
    }
    
    return td.keepParticle;
}

void Foam::solidParcel::LinearMove
(
    solidParcel& p,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    if(newParcel()!=-1)
    {
        // note: this justifies that freshly inserted parcels should be
        // tracked as if they had passed the boundary face on which they
        // have been inserted in the time step in which they are inserted.
        stepFraction() = td.cloud().rndGen().sample01<scalar>();
        newParcel() = -1;
    }

    const cloudSolution& solution = td.cloud().solution();
    const scalarField& cellLengthScale = td.cloud().cellLengthScale();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = solution.deltaTMax(trackTime);

    bool tracking = true;
    label nTrackingStalled = 0;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        if(td.cloud().PreduceD())
        {
            meshTools::constrainToMeshCentre(mesh, p.position());
        }
        
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Cache the parcel current cell as this will change if a face is hit
        const label celli = p.cell();

        const scalar magU = mag(U_);
        if (p.active() && tracking && (magU > ROOTVSMALL))
        {
            const scalar dl = dt*magU;
            const scalar deltaLMax = solution.deltaLMax(cellLengthScale[celli]);
            const scalar dCorr = min(dl, deltaLMax);
            dt *= dCorr/dl*p.trackToFace(p.position() + dCorr*(U_)/magU, td);
        }

        tEnd -= dt;

        const scalar newStepFraction = 1.0 - tEnd/trackTime;

        if (tracking)
        {
            if
            (
                mag(p.stepFraction() - newStepFraction)
            < particle::minStepFractionTol
            )
            {
                nTrackingStalled++;

                if (nTrackingStalled > 1)
                {
                    tracking = false;
                }
            }
            else
            {
                nTrackingStalled = 0;
            }
        }

        p.stepFraction() = newStepFraction;
        
        // Avoid problems with extremely small timesteps
        if ((dt > ROOTVSMALL))
        {
            // Update cell based properties
            p.setCellValues(td, dt, celli);

            if (solution.cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(td, dt, celli);
            }

            p.calc(td, dt, celli);
        }

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
            forAll(td.cloud().boundaries().solidCyclicBoundaryModels(), c)
            {
                const labelList& faces = td.cloud().boundaries().solidCyclicBoundaryModels()[c]->allFaces();

                if (findIndex(faces, face()) != -1)
                {
                    td.cloud().boundaries().solidCyclicBoundaryModels()[c]->controlParticle(*this, td);
                }
            }
        }

        
    }
}

bool Foam::solidParcel::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}

bool Foam::solidParcel::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label
)
{
    return false;
}

void Foam::solidParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //-find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().solidPatchToModelIds()[patchIndex];

    td.cloud().boundaries().solidPatchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::solidParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //- find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().solidPatchToModelIds()[patchIndex];

    //- apply a boundary model when a particle collides with this poly patch
    td.cloud().boundaries().solidPatchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::solidParcel::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::solidParcel::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::solidParcel::wallImpactDistance(const vector&) const
{
    return 0.5*d_;
}

// ************************************************************************* //
#include "solidParcelIO.C"
