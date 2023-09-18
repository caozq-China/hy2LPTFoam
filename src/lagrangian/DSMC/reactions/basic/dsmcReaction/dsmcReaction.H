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
    dsmcReaction

Description
    Basic/abstract class of a state controller

SourceFiles
    dsmcReactionI.H
    dsmcReaction.C
    dsmcReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcReaction_H
#define dsmcReaction_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "Random.H"
#include "dsmcParcel.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class dsmcReaction Declaration
\*---------------------------------------------------------------------------*/

class dsmcReaction
{
    // Private Member Functions

        //- Dissociate particle P|Q
        void dissociatePQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const;

        //- Ionise particle P|Q
        void ionisePQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the cloud database
        dsmcCloud& cloud_;

        //- Total number of reactions
        label nTotReactions_;

        //- Number of reactions per time step
        label nReactionsPerTimeStep_;

        //- Properties dictionary
        dictionary propsDict_;

        //- Reaction name
        word reactionName_;

        //- Reactant names
        wordPair reactants_;

        //- Reactant IDs
        labelPair reactantIds_;

        label rDof1_;
        label rDof2_;
        label vDof1_;
        label vDof2_;
        label charge1_;
        label charge2_;

        //- Relax flag
        bool relax_;

        // remove if not used for heat bath
        bool allowSplitting_;

        bool writeRatesToTerminal_;

        scalar volume_;

        List<scalar> numberDensities_;


    // Protected Member Functions

        void setCommonReactionProperties();


        //- Dissociate particle P
        void dissociateP
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const
        {
            dissociatePQ(heatOfReaction, productIDs, p, q);
        }

        //- Dissociate particle Q
        void dissociateQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const
        {
            dissociatePQ(heatOfReaction, reverse(productIDs), q, p);
        }

        //- Ionise particle P
        void ioniseP
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const
        {
            ionisePQ(heatOfReaction, productIDs, p, q);
        }

        //- Ionise particle Q
        void ioniseQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            dsmcParcel& p,
            dsmcParcel& q
        ) const
        {
            ionisePQ(heatOfReaction, reverse(productIDs), q, p);
        }

        //- Particle vibrational energy
        scalar EVib(const dsmcParcel& p) const;

        //- Particle rotational energy
        scalar ERot(const dsmcParcel& p) const;

        //- Particle electrical energy level
        scalar EEle(const dsmcParcel& p) const;

        //- Particle mass
        scalar m(const dsmcParcel& p) const;

        //- Reduced mass
        scalar mR(const dsmcParcel& p, const dsmcParcel& q) const;

        //- Tranlational energy
        scalar translationalEnergy
        (
            const dsmcParcel& p,
            const dsmcParcel& q
        ) const;

        //- Omega
        scalar omega(const dsmcParcel& p, const dsmcParcel& q) const;

        //- Particle thetaD
        scalar thetaD(const dsmcParcel& p) const;

        //- Particle thetaV
        scalar thetaV(const dsmcParcel& p) const;

        //- Particle Zref
        scalar Zref(const dsmcParcel& p) const;

        //- Particle refTempZv
        scalar refTempZv(const dsmcParcel& p) const;

        //- Species dissociation level
        label charDissLevel(const dsmcParcel& p) const;

        //- Species number of electronic levels
        label jMax(const dsmcParcel& p) const;

        //- Species rotational models
        label rotationalDof(const dsmcParcel& p) const;

        //- Centre of mass velocity for particle pair
        vector Ucm(const dsmcParcel& p, const dsmcParcel& q) const;

        //- Species electronic energies
        scalarList EEList(const dsmcParcel& p) const;

        //- Species degeneracies
        labelList gList(const dsmcParcel& p) const;

        //- Perform associative ionisation
        void associativeIonisation
        (
            const scalar& heatOfReactionIntermediateIonisation,
            const scalar& heatOfReactionRecombination,
            const labelPair& assIonProductIds,
            dsmcParcel& p,
            dsmcParcel& q
        );
        
        //- Perform exchange reaction
        void exchangePQ
        (
            const scalar& heatOfReactionExchJoules,
            const labelPair& exchangeProductIds,
            dsmcParcel& p,
            dsmcParcel& q
        );
        
        //- Perform charge exchange reaction
        void chargeExchangePQ
        (
            const scalar& heatOfReactionExchJoules,
            const labelPair& chargeExchangeProductIds,
            dsmcParcel& p,
            dsmcParcel& q
        );


public:

    //- Runtime type information
    TypeName("dsmcReaction");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcReaction,
            dictionary,
            (
                const Time& t,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );

    // Constructors

        //- Construct from components
        dsmcReaction
        (
            const Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcReaction> New
        (
            const Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dsmcReaction() = default;

    // Member Functions

        virtual void initialConfiguration() = 0;

        //- Apply a reaction between parcel p and q
        virtual void reaction(dsmcParcel& p, dsmcParcel& q) = 0;

        //- Return bool for relax
        //  true : run normal collision model
        //  false : apply reaction to particles
        virtual bool relax() const;

        //- Return bool if particles belong to this reaction model
        virtual bool tryReactMolecules
        (
            const label typeIdP,
            const label typeIdQ
        ) const = 0;

        virtual bool outputResults(const label counterIndex);

        //- Return const access to the cloud
        const dsmcCloud& cloud() const;

        //- Return const access to the total number of reactions
        inline label nTotReactions() const;

        //- Return const access to the number of reactions per time step
        inline label nReactionsPerTimeStep() const;

        //- Return access to the number of reactions per time step
        inline label& nReactionsPerTimeStep();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "dsmcReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
