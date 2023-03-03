/*---------------------------------------------------------------------------*\
Copyright (C) 2023-2023 University of Latvia
-------------------------------------------------------------------------------
License
    This file is part of ERapid-Library.
    
    ERapid-Library is based on EOF-Library (https://github.com/jvencels/eof-library)
    and adapted for use with RapidCFD (https://github.com/Atizar/RapidCFD-dev)
    instead of OpenFOAM.

    ERapid-Library is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ERapid-Library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ERapid-Library.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Elmer

Description
    Coupling with FEM open source multiphysical simulation software Elmer

SourceFiles
    Elmer.C

Original author:   Juris Vencels (EOF Consulting, Latvia)

Adaptation for RapidCFD by:       Didzis Berenis

Adaptation date: 03.03.2023.

\*---------------------------------------------------------------------------*/

#ifndef Elmer_H
#define Elmer_H

// Stuff for compatibility with OF2.4
#ifndef FatalErrorInFunction
#ifdef __GNUC__
   #define FUNCTION_NAME __PRETTY_FUNCTION__
#else
   #define FUNCTION_NAME __func__
#endif
 //- Report an error message using Foam::FatalError
 //  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)
#endif

#include "fvCFD.H"
#include "PstreamGlobals.H"
#include "findCellFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
typedef gpuField<vector> vectorgpuField;
typedef gpuList<vector> vectorgpuList;
typedef UList<vector> vectorUList;
typedef List<vector> vectorList;
typedef gpuList<label> labelgpuList;
typedef Field<label> labelField;

/*---------------------------------------------------------------------------*\
                          Class Elmer Declaration
\*---------------------------------------------------------------------------*/

template <class meshT>
class Elmer
{

public:

    // Public data

        //- Construct and initialize if needed
        Elmer(const meshT& mesh, int mode, bool init=true);

        //- Construct and initialize if needed
        Elmer(const meshT& mesh, int mode, point minPoint, point maxPoint, bool init=true);

        //- Initialize
        void initialize();

        //- Send interpolated scalar field at element centres to Elmer
        void sendScalar(volScalarField& field);

        //- Send interpolated vector field at element centres to Elmer
        void sendVector(volVectorField& field);

        //- Get interpolated scalar field at cell centres from Elmer
        void recvScalar(volScalarField& field);

        //- Get interpolated vector field at cell centres from Elmer
        void recvVector(volVectorField& field);//Foam::List<vector>& list_v);//

        //- Simulation status (1=ok, 0=lastIter, -1=error)
        void sendStatus(int status);

        //- Get bound box for local mesh partition
        void getboundBox(const meshT& mesh);

        //- Compare bounding box of this OF partition with all EL ranks
        void findOverlappingBoxes();

        //- Test for incoming message that consumes less CPU
        void MPI_Test_Sleep(MPI_Request& req);

        //- Data structure for each Elmer MPI rank
        class ElmerProc_t
        {
            public:
                // MPI stuff
		        MPI_Request reqSend, reqRecv;
		        int globalRank;

                // Receiving fields from Elmer
		        int *foundCellsIndx;
		        int nFoundCells;
				double *recvBuffer0;

                // For interpolation
		int **interpCellIndx;
        	point *interpCoeff;
		labelgpuList foundElmerGPUIndx;

                // Sending fields to Elmer
                int nElem;
                double *sendBuffer0, *sendBuffer1, *sendBuffer2;
		
                label *foundElement;
                int *foundElementIndx;
                int *foundElementCellIndx;
		        int nFoundElements;
                point *positions;

                // Bounding box
                bool boxOverlap;
        };


private:

    // Private data

        //- Find cells
        void findCells();//new

        //- Prepare for interpolation
        void interpPrepare();//new

        //- Number of cells for this OpenFOAM MPI rank
        const meshT& mesh_;

        //- Coupling mode: 1=send, -1=receive, 0=both
        const int mode_;

        //- Number of cells for this OpenFOAM MPI rank
        const int nCells;
        
        //- Number of cells for this OpenFOAM MPI rank as label
	const label Lsize;

        //- My rank's cell centres
	const vectorgpuField  cellCentres_;
	
        //- My rank's cell centres
        double *cellCentres_x, *cellCentres_y, *cellCentres_z;

	int nCommElem;

	double localTime;

        //- list for keeping scalar values of field
	scalarList list_s;
	
        //- list for keeping vector values of field
	vectorList list_v;
	
        //- reference list used for setting scalar values to gpuField
	const scalarUList& Ulist_s;
	
        //- reference list used for setting vector values to gpuField
	const vectorUList& Ulist_v;
	
        //- Total global MPI ranks
        int totGlobalRanks;

        //- My rank at MPI_COMM_WORLD
        int myGlobalRank;

        //- Total OpenFOAM ranks
        int totLocalRanks;

        //- My rank at MPI_OF_WORLD
        int myLocalRank;

        //- Total ELmer ranks
        int totElmerRanks;
    
        //- Elmer ranks start in MPI_COMM_WORLD
        int ElmerRanksStart;

        //- OpenFOAM ranks start in MPI_COMM_WORLD
        int OFRanksStart;

        //- Data structure for each Elmer MPI rank
        ElmerProc_t *ELp;

        // Interpolation types used for interpolating Elmer's elemental field
        dictionary interpolationDict;

        // Bounding box for this OF rank
        boundBox myBoundBox;

        //- Bounding boxes for all OF ranks [rank][min,max][x,y,z]
        double *OFboundBoxes;

        //- Bounding boxes for all EL ranks [rank][min,max][x,y,z]
        double *ELboundBoxes;

        //- Boolean matrix for storing overlappings [OFranks][ELranks]
        int *OF_EL_overlap;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
