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

#include "dsmcMorrisNozzleExitPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcMorrisNozzleExitPatch, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcMorrisNozzleExitPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMorrisNozzleExitPatch::dsmcMorrisNozzleExitPatch
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    velocities_(),
    transTemp_(),
    rotatTemp_(),
    vibrTemp_(),
    elecTemp_(),
    nozzleExitRadius_(),
    numberDensities_(),
    rhoX_(17,0.0),
    rhoy_(17,0.0),
    TX_(17,0.0),
    TY_(17,0.0),
    axialUx_(8,0.0),
    axialUy_(8,0.0),
    radialVx_(12,0.0),
    radialVy_(12,0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;
    
    //- Do not put "setProperties()" here. It will cause an error.
    //- Because the constProps in dsmcCloud has not been built during constructing.
//     setProperties();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcMorrisNozzleExitPatch::initialConfiguration()
{
    setProperties();
}

void dsmcMorrisNozzleExitPatch::calculateProperties()
{}

void dsmcMorrisNozzleExitPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        transTemp_,
        velocities_,
        numberDensities_
    );
    
    insertParcels
    (
        transTemp_,
        rotatTemp_,
        vibrTemp_,
        elecTemp_,
        velocities_
    );
    
}

void dsmcMorrisNozzleExitPatch::controlParcelsBeforeCollisions()
{}

void dsmcMorrisNozzleExitPatch::controlParcelsAfterCollisions()
{}

void dsmcMorrisNozzleExitPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void dsmcMorrisNozzleExitPatch::updateProperties(const dictionary& dict)
{
    //- the main properties should be updated first
    dsmcGeneralBoundary::updateProperties(dict);
}



void dsmcMorrisNozzleExitPatch::setProperties()
{
    initialiseRhoX();
    initialiseRhoY();
    initialiseTX();
    initialiseTY();
    initialiseUx();
    initialiseUy();
    initialiseVx();
    initialiseVy();
    
    nozzleExitRadius_ = propsDict_.get<scalar>("nozzleExitRadius");

    typeIds_ = cloud_.getTypeIDs(propsDict_);
    
    // Read in the mass density per specie

//     const dictionary& numberDensitiesDict
//     (
//         propsDict_.subDict("numberDensities")
//     );
    
    numberDensities_.clear();
    numberDensities_.setSize(typeIds_.size());
    forAll(numberDensities_, i)
    {
        numberDensities_[i].setSize(faces_.size(), 0.0);
    }
    
    velocities_.setSize(faces_.size(), Zero);
    transTemp_.setSize(faces_.size(), 0.0);
    rotatTemp_.setSize(faces_.size(), 0.0);
    vibrTemp_.setSize(faces_.size(), 0.0);
    elecTemp_.setSize(faces_.size(), 0.0);

    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const vector fC = mesh_.faceCentres()[faceI];
        
        scalar RatioRR0 = fC.y()/nozzleExitRadius_;
        velocities_[f].x() = linearInterpolation(axialUx_,axialUy_,RatioRR0);
        velocities_[f].y() = velocities_[f].x()*linearInterpolation(radialVx_,radialVy_,RatioRR0)*(-0.01);
        velocities_[f].z() = 0.0;
        
        transTemp_[f] = linearInterpolation(TX_,TY_,RatioRR0);
        rotatTemp_[f] = transTemp_[f];
        vibrTemp_[f] = transTemp_[f];
        elecTemp_[f] = transTemp_[f];
        
        forAll(numberDensities_, i)
        {
            const label typeId = typeIds_[i];
            const scalar& mass = cloud_.constProps(typeId).mass();
            numberDensities_[i][f] = linearInterpolation(rhoX_,rhoy_,RatioRR0)/mass;
        }
    }
    
//     Info<<numberDensities_<<endl;

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }
}


scalar dsmcMorrisNozzleExitPatch::linearInterpolation(scalarList& listX,scalarList& listY, scalar Rratio)
{
    scalar Yinterplotated = 0.0;
    
    for(label i=0; i<(listX.size()-1);i++)
    {
        if(Rratio > listX[i] && Rratio < listX[i+1])
        {
            Yinterplotated = listY[i+1]-(listX[i+1]-Rratio)*(listY[i+1]-listY[i])/(listX[i+1]-listX[i]);
            break;
        }
        else if( Rratio == listX[i])
        {
            Yinterplotated = listY[i];
            break;
        }
    }
    
    
    return Yinterplotated;
}



void dsmcMorrisNozzleExitPatch::initialiseRhoX()
{
    //-initialise rhoX_
    rhoX_[0] = 0.0; rhoX_[1] = 0.19826087; rhoX_[2] = 0.3426087; rhoX_[3] = 0.37391304;
    rhoX_[4] = 0.39826087; rhoX_[5] = 0.43130435; rhoX_[6] = 0.45217391; rhoX_[7] = 0.49217391;
    rhoX_[8] = 0.57043478; rhoX_[9] = 0.63478261; rhoX_[10] = 0.88173913; rhoX_[11] = 0.91304348;
    rhoX_[12] = 0.9426087; rhoX_[13] = 0.95478261; rhoX_[14] = 0.96521739; rhoX_[15] = 0.97391304;
    rhoX_[16] = 1.0;
}


void dsmcMorrisNozzleExitPatch::initialiseRhoY()
{
    rhoy_[0] = 0.00037343; rhoy_[1] = 0.00036503; rhoy_[2] = 0.00035105; rhoy_[3] = 0.00034266;
    rhoy_[4] = 0.00036503; rhoy_[5] = 0.00052448; rhoy_[6] = 0.0005972; rhoy_[7] = 0.00065594;
    rhoy_[8] = 0.00078182; rhoy_[9] = 0.00091888; rhoy_[10] = 0.00157343; rhoy_[11] = 0.00163776;
    rhoy_[12] = 0.00167972; rhoy_[13] = 0.00168811; rhoy_[14] = 0.00168252; rhoy_[15] = 0.00160979;
    rhoy_[16] = 0.00089371;
}

void dsmcMorrisNozzleExitPatch::initialiseTX()
{
    TX_[0] = 0.0; TX_[1] = 0.25018851; TX_[2] = 0.37020067; TX_[3] = 0.38932198;
    TX_[4] = 0.40492247; TX_[5] = 0.43089389; TX_[6] = 0.45868957; TX_[7] = 0.60465491;
    TX_[8] = 0.75233506; TX_[9] = 0.85311645; TX_[10] = 0.92783825; TX_[11] = 0.95910307;
    TX_[12] = 0.96949833; TX_[13] = 0.97464275; TX_[14] = 0.98479173; TX_[15] = 1.0;
}

void dsmcMorrisNozzleExitPatch::initialiseTY()
{
    TY_[0] = 646.85314685; TY_[1] = 646.60839161; TY_[2] = 634.61538462; TY_[3] = 639.86013986;
    TY_[4] = 669.58041958; TY_[5] = 736.01398601; TY_[6] = 753.4965035; TY_[7] = 823.42657343;
    TY_[8] = 907.34265734; TY_[9] = 958.04195804; TY_[10] = 993.00699301; TY_[11] = 1015.73426573;
    TY_[12] = 1038.46153846; TY_[13] = 1080.41958042; TY_[14] = 1244.75524476; TY_[15] = 1283.21678322;
}


void dsmcMorrisNozzleExitPatch::initialiseUx()
{
    axialUx_[0] = 0.0; axialUx_[1] = 0.10086957; axialUx_[2] = 0.36521739; axialUx_[3] = 0.47304348;
    axialUx_[4] = 0.94956522; axialUx_[5] = 0.96695652; axialUx_[6] = 0.97565217; axialUx_[7] = 1.0;
}

void dsmcMorrisNozzleExitPatch::initialiseUy()
{
    axialUy_[0] = -3291.95804196; axialUy_[1] = -3291.95804196; axialUy_[2] = -3243.00699301; axialUy_[3] = -3194.05594406;
    axialUy_[4] = -2979.8951049; axialUy_[5] = -2961.53846154; axialUy_[6] = -2881.99300699; axialUy_[7] = -1321.67832168;
}

void dsmcMorrisNozzleExitPatch::initialiseVx()
{
    radialVx_[0] = 0.0; radialVx_[1] = 0.331404; radialVx_[2] = 0.374819; radialVx_[3] = 0.428365;
    radialVx_[4] = 0.505065; radialVx_[5] = 0.610709; radialVx_[6] = 0.694645; radialVx_[7] = 0.775687;
    radialVx_[8] = 0.859624; radialVx_[9] = 0.936324; radialVx_[10] = 0.975398; radialVx_[11] = 1.0;
}

void dsmcMorrisNozzleExitPatch::initialiseVy()
{
    radialVy_[0] = 1.007194; radialVy_[1] = 17.42446; radialVy_[2] = 16.719424; radialVy_[3] = 13.395683;
    radialVy_[4] = 14.705036; radialVy_[5] = 16.820144; radialVy_[6] = 16.215827; radialVy_[7] = 17.323741;
    radialVy_[8] = 18.733813; radialVy_[9] = 22.964029; radialVy_[10] = 31.928058; radialVy_[11] = 69.798561;
}

} // End namespace Foam

// ************************************************************************* //
