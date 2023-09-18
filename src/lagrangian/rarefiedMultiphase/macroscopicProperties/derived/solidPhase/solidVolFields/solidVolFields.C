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

Measures overall temperature, including vibrational temperature, for a single 
species gas or a gas mixture and writes the results to a volume scalar field 
that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be 
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements
class and are also written.

\*---------------------------------------------------------------------------*/

#include "solidVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidVolFields, 0);

addToRunTimeSelectionTable(solidField, solidVolFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidVolFields::solidVolFields
(
    const Time& t,
    const polyMesh& mesh,
    solidParticleCouplingCloud& cloud,
    const dictionary& dict
)
:
    solidField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    n_(),
    t1_(),
    t2_(),
    sampleInterval_(1),
    sampleCounter_(0),
//     mfpReferenceTemperature_(273.0),
    fieldName_(propsDict_.get<word>("fieldName")),
    solidRhoN_
    (
        IOobject
        (
            "solidRhoN_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    solidRhoNMean_
    (
        IOobject
        (
            "solidRhoNMean_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    loadBalanceWeight_
    (
        IOobject
        (
            "loadBalanceWeight_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    overallTsolid_
    (
        IOobject
        (
            "overallTsolid_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    arithmeticAveT_
    (
        IOobject
        (
            "arithmeticAveT_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    surfaceNumberFlux_
    (
        IOobject
        (
            "surfaceNumberFlux_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime/dimLength/dimLength, Zero)
    ),
    surfaceMassFlux_
    (
        IOobject
        (
            "surfaceMassFlux_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime/dimLength/dimLength, Zero)
    ),
    surfaceConductiveHeatFlux_
    (
        IOobject
        (
            "surfaceConductiveHeatFlux_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime/dimLength/dimLength, Zero)
    ),
    //- 2021/3/21
    liquidMassFraction_
    (
        IOobject
        (
            "liquidMassFraction_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
//     volFraction_
//     (
//         IOobject
//         (
//             "volFraction_"+ fieldName_,
//             mesh_.time().timeName(),
//             mesh_,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh_,
//         dimensionedScalar(dimless/dimVolume, Zero)
//     ),
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/dimTime, Zero)
    ),
    nTimeSteps_(0),
    typeIdSolids_(),
//     nTotColls_
//     (
//         IOobject
//         (
//             "nTotInterphsColls_"+fieldName_,
//             mesh_.time().timeName(),
//             mesh_,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         Pstream::nProcs()
//     ),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNInstantaneous_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    liquidMassXnParticle_(mesh_.nCells(), 0.0),
//     volumeXnParticle_(mesh_.nCells(), 0.0),
    cXm_(mesh_.nCells(), 0.0),
    cXmXT_(mesh_.nCells(), 0.0),
    TXnParticle_(mesh_.nCells(), 0.0),
    nColls_(mesh_.nCells(), 0.0),
    momentumMean_(mesh_.nCells(), Zero),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    nTotColls_
    (
        IOobject
        (
            "nTotColls",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    boundaryCells_(),
    nParcels_(),
    nParcelsXnParticle_(),
    rhoNBF_(),
    rhoMBF_(),
//     volFractionBF_(),
    numberFluxBF_(),
    massFluxBF_(),
    conductiveHeatFluxBF_(),
    momentumBF_(),
    averagingAcrossManyRuns_(false),
    densityOnly_(false)
{

    // standard to reading typeIds ------------ 
    typeIdSolids_ = cloud_.getTypeIDs(propsDict_);
    
    nParcels_.setSize(typeIdSolids_.size());

    for (auto& n : nParcels_)
    {
        n.setSize(mesh_.nCells());
    }
    
    nParcelsXnParticle_.setSize(typeIdSolids_.size());
    
    for (auto& n : nParcelsXnParticle_)
    {
        n.setSize(mesh_.nCells());
    }
    
    nTotColls_.setSize(Pstream::nProcs(),0.0);
    
    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }
    
    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
//     volFractionBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    numberFluxBF_.setSize(mesh_.boundaryMesh().size());
    massFluxBF_.setSize(mesh_.boundaryMesh().size());
    conductiveHeatFluxBF_.setSize(mesh_.boundaryMesh().size());
    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());
        
    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        
        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
//         volFractionBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        numberFluxBF_[j].setSize(patch.size(), 0.0);
        massFluxBF_[j].setSize(patch.size(), 0.0);
        conductiveHeatFluxBF_[j].setSize(patch.size(), 0.0);
        
        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }
    
    calculateWallUnitVectors();
    
    
    sampleInterval_ =
        propsDict_.readIfPresent<label>("sampleInterval", sampleInterval_);
    
    if (propsDict_.readIfPresent<bool>("densityOnly", densityOnly_))
    {
        if (densityOnly_)
        {
            Info<< "densityOnly initiated" << endl;
        }
    }
    
    if
    (
        propsDict_.readIfPresent<bool>
        (
            "averagingAcrossManyRuns",
            averagingAcrossManyRuns_
        )
    )
    {
        // read in stored data from dictionary
        if (averagingAcrossManyRuns_)
        {
            Info<< "averagingAcrossManyRuns initiated." << nl << endl;
            readIn();
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidVolFields::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "volFieldsMethod_"+fieldName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );
    
    dict.readIfPresent("nTimeSteps", nTimeSteps_);
    dict.readIfPresent("rhoNMean", rhoNMean_);
    dict.readIfPresent("rhoNInstantaneous", rhoNInstantaneous_);
    dict.readIfPresent("rhoNMeanXnParticle", rhoNMeanXnParticle_);
    dict.readIfPresent("rhoMMean", rhoMMean_);
    dict.readIfPresent("rhoMMeanXnParticle", rhoMMeanXnParticle_);
    dict.readIfPresent("liquidMassXnParticle", liquidMassXnParticle_);
//     dict.readIfPresent("volumeXnParticle", volumeXnParticle_);
    dict.readIfPresent("cXm", cXm_);
    dict.readIfPresent("cXmXT", cXmXT_);
    dict.readIfPresent("TXnParticle", TXnParticle_);
    dict.readIfPresent("nColls", nColls_);
    dict.readIfPresent("momentumMean", momentumMean_);
    dict.readIfPresent("momentumMeanXnParticle", momentumMeanXnParticle_);
    dict.readIfPresent("nParcels", nParcels_);
    dict.readIfPresent("nParcelsXnParticle", nParcelsXnParticle_);
    dict.readIfPresent("rhoNBF", rhoNBF_);
    dict.readIfPresent("rhoMBF", rhoMBF_);
//     dict.readIfPresent("volFractionBF", volFractionBF_);
    dict.readIfPresent("numberFluxBF", numberFluxBF_); 
    dict.readIfPresent("massFluxBF_", massFluxBF_); 
    dict.readIfPresent("conductiveHeatFluxBF_", conductiveHeatFluxBF_); 
    dict.readIfPresent("momentumBF", momentumBF_);
}

void solidVolFields::writeOut()
{
    if (mesh_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_"+fieldName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);
        dict.add("rhoNMean", rhoNMean_);
        dict.add("rhoNInstantaneous", rhoNInstantaneous_);
        dict.add("rhoNMeanXnParticle", rhoNMeanXnParticle_);
        dict.add("rhoMMean", rhoMMean_);
        dict.add("rhoMMeanXnParticle", rhoMMeanXnParticle_);
        dict.add("liquidMassXnParticle", liquidMassXnParticle_);
//         dict.add("volumeXnParticle", volumeXnParticle_);
        dict.add("cXm", cXm_);
        dict.add("cXmXT", cXmXT_);
        dict.add("TXnParticle", TXnParticle_);
        dict.add("nColls", nColls_);
        dict.add("momentumMean", momentumMean_);
        dict.add("momentumMeanXnParticle", momentumMeanXnParticle_);
        dict.add("nParcels", nParcels_);
        dict.add("nParcelsXnParticle", nParcelsXnParticle_);
        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
//         dict.add("volFractionBF", volFractionBF_);
        dict.add("numberFluxBF", numberFluxBF_); 
        dict.add("massFluxBF_", massFluxBF_); 
        dict.add("conductiveHeatFluxBF_", conductiveHeatFluxBF_); 
        dict.add("momentumBF", momentumBF_);
        
        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}

void solidVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];
       
        if(isA<wallPolyPatch>(pPatch))
        {           
            const vectorField& fC = pPatch.faceCentres();
           
            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei];
                n_[patchi][facei].normalise();
               
                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] - 
                            mesh_.points()[mesh_.faces()[pPatch.start() 
                            + facei][0]];
                t1_[patchi][facei].normalise();
               
                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei].normalise();
            }
        }
    }
}

//- initial conditions
void solidVolFields::createField()
{
    Info << "Initialising solidVolFields field" << endl;
}


void solidVolFields::calculateField()
{ 
    ++sampleCounter_;
    
    const scalar nParticle = cloud_.nSolidParticles();
    
    rhoNInstantaneous_ = 0.0;
    loadBalanceWeight_ = 0.0;
    
    if(sampleInterval_ <= sampleCounter_)
    {
        ++nTimeSteps_;
        
        if(densityOnly_)
        {
            forAllConstIters(cloud_, iter)
            {
                const solidParticleCoupling& p = iter();

                if(typeIdSolids_.find(p.typeID()) != -1)
                {
                    const label cell = p.cell();
                    //cell not cellI in parcel. It can be found in particle.H
                    //- calculate individual solid particle mass
                    const scalar mass = cloud_.constSolidProps(p.typeID()).massSphere();

                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    
                    scalar RWF =
                        cloud_.dsmcCloudReference()->axiRWF(cloud_.mesh().cellCentres()[cell]);

                        rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                }
            }
        }
        else
        {
            forAllConstIters(cloud_, iter)
            {
                const solidParticleCoupling& p = iter();
                const label& iD = typeIdSolids_.find(p.typeID());

                if(iD != -1)
                {
                    const label cell = p.cell();
                    const scalar mass = cloud_.constSolidProps(p.typeID()).massSphere();     
                    const vector& U = p.U();
                    const scalar CzRatio = p.CzRatio();
                                    
                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    rhoMMean_[cell] += mass;
                    momentumMean_[cell] += mass*U;
                    nParcels_[iD][cell] += 1.0;
                    
                    scalar RWF =
                        cloud_.dsmcCloudReference()->axiRWF(cloud_.mesh().cellCentres()[cell]);

                        nParcelsXnParticle_[iD][cell] += RWF*nParticle;//- for specific species
                        rhoNMeanXnParticle_[cell] += RWF*nParticle;
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                        
//                         scalar volume = pi*pow(p.D(),3.0)/6.0;
//                         volumeXnParticle_[cell] += (volume*RWF*nParticle);
                        
                        //- calculate particle mass of liquid phase
                        if(cloud_.solidPhaseChange().active())
                        {
                            if(p.phaseState()==1)
                            {
                                //- phaseState = 1, the particle temperature is lower than
                                //- the melting temperature  so that the 
                                //- density is not corrected through liquid phase density equation
                                liquidMassXnParticle_[cell]+=(pow(CzRatio,3.0)*mass*RWF*nParticle);
                            }
                            else if(p.phaseState()==2)
                            {
                                //- the case of melting
                                liquidMassXnParticle_[cell]+= mass-(1.0/6.0)*pi*pow(CzRatio*cloud_.constSolidProps(p.typeID()).d(),3)*cloud_.constSolidProps(p.typeID()).rho();
                            }
                            else if(p.phaseState()==3)
                            {
                                //- pure liquid 
                                liquidMassXnParticle_[cell]+= mass*RWF*nParticle;
                            }
                            
                        }
                        
                        cXm_[cell] += cloud_.constSolidProps(p.typeID()).Cp()*mass*RWF*nParticle;
                        cXmXT_[cell] += cloud_.constSolidProps(p.typeID()).Cp()*mass*p.T()*RWF*nParticle;
                        momentumMeanXnParticle_[cell] += (mass*(U)*RWF*nParticle);
                        TXnParticle_[cell] += p.T()*RWF*nParticle;
                }
            }
            
            
            // Obtain collision quality measurements

            forAll(cloud_.solidCellPropMeasurements().nColls(), cell)
            {
//                 collisionSeparation_[cell] +=
//                     cloud_.cellPropMeasurements().collisionSeparation()[cell];
                nColls_[cell] += cloud_.solidCellPropMeasurements().nColls()[cell];
            }
            
            forAll(cloud_.solidCellPropMeasurements().nCollsTotal(),procNo)
            {
                nTotColls_[procNo] += cloud_.solidCellPropMeasurements().nCollsTotal()[procNo];
            }
            
            // obtain boundary measurements
            auto& bm = cloud_.boundaryFluxMeasurements();
            forAll(bm.rhoNBF(), i)
            {
                const label iD = typeIdSolids_.find(i);
                
                if(iD != -1)
                { 
                    forAll(bm.rhoNBF()[i], j)
                    {                
                        forAll(bm.rhoNBF()[i][j], k)
                        {
                            rhoNBF_[j][k] += bm.rhoNBF()[i][j][k];
                            rhoMBF_[j][k] += bm.rhoMBF()[i][j][k];
//                             volFractionBF_[j][k] += bm.volumeFractionBF()[i][j][k];
                            momentumBF_[j][k] += bm.momentumBF()[i][j][k];
                            numberFluxBF_[j][k] += bm.numberFluxBF()[i][j][k];
                            massFluxBF_[j][k] += bm.massFluxBF()[i][j][k];
                            conductiveHeatFluxBF_[j][k] += bm.conductiveHeatFluxBF()[i][j][k];
                        }
                    }
                }
            }
        }
        sampleCounter_ = 0;
    }
    
    if(mesh_.time().writeTime())
    {
        const scalar nAvTimeSteps = nTimeSteps_;
        
        if(densityOnly_)
        {
            forAll(rhoNMean_, cell)
            {
                if(rhoNMean_[cell] > VSMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    
                    rhoN_[cell] = 
                        (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = 
                        (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);

                }
                else
                {
                     // not zero so that weighted decomposition still works
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }
                
                if (rhoNInstantaneous_[cell] > VSMALL)
                {
                    solidRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    solidRhoN_[cell] = 0.001;
                }
                
                if(cloud_.dsmcCloudReference()->fields().dsmcWeight()[cell] > 0.001
                    ||solidRhoN_[cell] > 0.001)
                {
                    loadBalanceWeight_[cell] = 
                        cloud_.dsmcCloudReference()->fields().dsmcWeight()[cell]
                        + cloud_.solidWeightingAmplificationFactor()*nColls_[cell];
                }
                else
                {
                    loadBalanceWeight_[cell] = 0.001;
                }
            }
        }
        else
        {                  
            
            forAll(rhoNMean_, cell)
            {                
                if(rhoNMean_[cell] > VSMALL)
                {                  
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    
                    rhoN_[cell] = 
                        (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = 
                        (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    scalar rhoMMean = 
                        rhoMMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);
                    UMean_[cell] = momentumMeanXnParticle_[cell] / 
                                    (rhoMMean*cellVolume*nAvTimeSteps);
                    
                    if(cloud_.solidPhaseChange().active())
                    {
                        liquidMassFraction_[cell] = liquidMassXnParticle_[cell]/rhoMMeanXnParticle_[cell];
                    }
                    
//                     volFraction_[cell] = volumeXnParticle_[cell]/(cellVolume*nAvTimeSteps);
                                    
                    arithmeticAveT_[cell] = TXnParticle_[cell]/rhoNMeanXnParticle_[cell];

                    overallTsolid_[cell] = cXmXT_[cell]/cXm_[cell];

                }
                else
                {
                    // not zero so that weighted decomposition still works
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    overallTsolid_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
//                     volFraction_[cell] = 0.0;
                    if(cloud_.solidPhaseChange().active())
                    {
                        liquidMassFraction_[cell] = 0.0;
                    }
                    arithmeticAveT_[cell] = 0.0;
                }
                
                if (rhoNInstantaneous_[cell] > VSMALL)
                {
                    solidRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    solidRhoN_[cell] = 0.001;
                }
                
                if(cloud_.dsmcCloudReference()->fields().dsmcWeight()[cell] > 0.001
                    ||solidRhoN_[cell] > 0.001)
                {
                    loadBalanceWeight_[cell] = 
                        cloud_.dsmcCloudReference()->fields().dsmcWeight()[cell]
                        + cloud_.solidWeightingAmplificationFactor()*nColls_[cell];
                }
                else
                {
                    loadBalanceWeight_[cell] = 0.001;
                }
                
            }
            
            // computing boundary measurements
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                if(isA<wallPolyPatch>(patch))
                {                               
                    forAll(rhoN_.boundaryFieldRef()[j], k)
                    {                        
                        rhoN_.boundaryFieldRef()[j][k] = 
                                rhoNBF_[j][k]*nParticle/nAvTimeSteps;
                        rhoM_.boundaryFieldRef()[j][k] = 
                                rhoMBF_[j][k]*nParticle/nAvTimeSteps;
                                
//                         volFraction_.boundaryFieldRef()[j][k] = 
//                                 volFractionBF_[j][k]*nParticle/nAvTimeSteps;
                        
                        if(rhoM_.boundaryFieldRef()[j][k] > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = 
                                    momentumBF_[j][k]*nParticle/
                                    (rhoM_.boundaryFieldRef()[j][k]*nAvTimeSteps);
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                        }
                        
                        const scalar& deltaT = mesh_.time().deltaTValue();

                        surfaceNumberFlux_.boundaryFieldRef()[j][k] = numberFluxBF_[j][k]/(nAvTimeSteps*deltaT);
                        surfaceMassFlux_.boundaryFieldRef()[j][k] = massFluxBF_[j][k]/(nAvTimeSteps*deltaT);
                        
                        surfaceConductiveHeatFlux_.boundaryFieldRef()[j][k] = conductiveHeatFluxBF_[j][k]/(nAvTimeSteps*deltaT);

                    }


                }
            }
            
            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                
                const labelList& bCs = boundaryCells_[j];
                
                forAll(bCs, k)
                {
                    if
                    (
                        isA<polyPatch>(patch)
                     && !isA<emptyPolyPatch>(patch)
                     && !isA<cyclicPolyPatch>(patch)
                    )
                    {
                        rhoN_.boundaryFieldRef()[j][k] = 
                                            rhoN_[bCs[k]];
                        rhoM_.boundaryFieldRef()[j][k] = 
                                            rhoM_[bCs[k]];
//                         volFraction_.boundaryFieldRef()[j][k] = 
//                                             volFraction_[bCs[k]];
                        
                                    
                        if(!isA<wallPolyPatch>(patch))
                        {    

                            UMean_.boundaryFieldRef()[j][k] = 
                                            UMean_[bCs[k]]; 
                            
                            overallTsolid_.boundaryFieldRef()[j][k] = 
                                            overallTsolid_[bCs[k]]; 
                                            
                            arithmeticAveT_.boundaryFieldRef()[j][k] = 
                                            arithmeticAveT_[bCs[k]];
                        }
                        
                        
                    }
                }
            }
            
            overallTsolid_.write();
            arithmeticAveT_.write();
//             volFraction_.write();
            UMean_.write();
            surfaceMassFlux_.write();
            surfaceNumberFlux_.write();
            surfaceConductiveHeatFlux_.write();
            if(cloud_.solidPhaseChange().active())
            {
                liquidMassFraction_.write();
            }
            
            nTotColls_.write();
            
            
        }
        
        //- reset
        if(timeVel_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0;
            
            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = scalar(0.0);
                rhoMMean_[c] = scalar(0.0);
                momentumMean_[c] = vector::zero;
                rhoNMeanXnParticle_[c] = scalar(0.0);
                rhoMMeanXnParticle_[c] = scalar(0.0);
//                 volumeXnParticle_[c] = scalar(0.0);
                cXm_[c] = scalar(0.0);
                cXmXT_[c] = scalar(0.0);
                momentumMeanXnParticle_[c] = vector::zero;
                TXnParticle_[c] = scalar(0.0);
                liquidMassXnParticle_[c] = scalar(0.0);
            }
            
            
            // reset boundary information
            
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
//                 volFractionBF_[j] = 0.0;
                momentumBF_[j] = vector::zero;
                numberFluxBF_[j] = 0.0;
                massFluxBF_[j] = 0.0;
                conductiveHeatFluxBF_[j] = 0.0;
                
            }
            
            forAll(cloud_.solidCellPropMeasurements().nCollsTotal(),procNo)
            {
                nTotColls_[procNo] = 0.0;
            }
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
        
    }
}


//- write field
void solidVolFields::writeField()
{}


// scalarField& solidVolFields::rhoNInstantaneous()
// {
//     return rhoNInstantaneous_;
// }

// volVectorField& solidVolFields::UMean()
// {
//     return UMean_;
// }
// 
// volScalarField& solidVolFields::rhoM()
// {
//     return rhoM_;
// }

// word& solidVolFields::fieldName()
// {
//     return fieldName_;
// }

void solidVolFields::updateProperties(const dictionary& dict)
{
    //- the main properties should be updated first
    solidField::updateProperties(dict);
    
    propsDict_ = dict.subDict(typeName + "Properties");
    
    Info << "Update Field Properties !"<< endl;
    
    propsDict_.readIfPresent("densityOnly", densityOnly_);
    if (densityOnly_)
    {
        Info<< nl << "densityOnly initiated." << nl << endl;
    }
    
    propsDict_.readIfPresent
    (
        "averagingAcrossManyRuns",
        averagingAcrossManyRuns_
    );
    if (averagingAcrossManyRuns_)
    {
        Info << "averagingAcrossManyRuns initiated." << endl;
    }

}


} // End namespace Foam

// ************************************************************************* //

