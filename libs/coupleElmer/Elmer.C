/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include "Elmer.H"
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class meshT>
Foam::Elmer<meshT>::Elmer(const meshT& mesh, int mode, bool init)
:
mesh_(mesh),
mode_(mode),
Ulist_s(list_s),
Ulist_v(list_v),
cellCentres_(mesh.C().internalField()),
Lsize(mesh.C().size()),
nCells(mesh.C().size()),
myBoundBox(mesh.points(),false)
{
	list_s.setSize(Lsize, 0 );
	list_v.setSize(Lsize, vector(0,0,0) );
	if(init) initialize();
}

template <class meshT>
Foam::Elmer<meshT>::Elmer(const meshT& mesh, int mode, point minPoint, point maxPoint, bool init)
:
mesh_(mesh),
mode_(mode),
Ulist_s(list_s),
Ulist_v(list_v),
cellCentres_(mesh.C().internalField()),
Lsize(mesh.C().size()),
nCells(mesh.C().size()),
myBoundBox(minPoint,maxPoint)
{
	list_s.setSize(Lsize, 0 );
	list_v.setSize(Lsize, vector(0,0,0) );
	if(init) initialize();
}

template <class meshT>
void Foam::Elmer<meshT>::initialize()
{
    int i, j, k;
    double commTime = MPI_Wtime();
    localTime = 0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totGlobalRanks);

    MPI_Comm_size(PstreamGlobals::MPI_COMM_FOAM, &totLocalRanks);
    MPI_Comm_rank(PstreamGlobals::MPI_COMM_FOAM, &myLocalRank);

    totElmerRanks = totGlobalRanks-totLocalRanks;

    if (myLocalRank==0) OFRanksStart = myGlobalRank;

    MPI_Bcast(&OFRanksStart, 1, MPI_INT, 0, PstreamGlobals::MPI_COMM_FOAM);

    if (OFRanksStart==0) {
        ElmerRanksStart = totLocalRanks;
    } else {
        ElmerRanksStart = 0;
    }

    ELp = new ElmerProc_t[totElmerRanks];

    if (ELp == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        ELp[i].globalRank = i+ElmerRanksStart;
    }

    OFboundBoxes = new double[totLocalRanks*2*3];
    ELboundBoxes = new double[totElmerRanks*2*3];
    OF_EL_overlap = new int[totLocalRanks*totElmerRanks];
    findOverlappingBoxes();

		
    scalargpuField cx_ = cellCentres_ & vector(1,0,0);
    scalargpuField cy_ = cellCentres_ & vector(0,1,0);
    scalargpuField cz_ = cellCentres_ & vector(0,0,1);
    
    cellCentres_x = new double[nCells];
    cellCentres_y = new double[nCells];
    cellCentres_z = new double[nCells];

    if (cellCentres_x == nullptr || cellCentres_y == nullptr || cellCentres_z == nullptr) {
		FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
    }
    else
    {
		for ( k=0; k<nCells; k++ )
		{ 
			cellCentres_x[k] = cx_.get(k);
			cellCentres_y[k] = cy_.get(k);
			cellCentres_z[k] = cz_.get(k);
		}
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Receiving from Elmer
    if (mode_==-1) {

        Info<< "Sending data to Elmer.." << endl;

        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Isend(&nCells, 1, MPI_INT, ELp[i].globalRank, 999, MPI_COMM_WORLD, &ELp[i].reqSend);
        }
	
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);

            MPI_Isend(cellCentres_x, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);
            MPI_Request_free(&ELp[i].reqSend);
            MPI_Isend(cellCentres_y, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);
            MPI_Request_free(&ELp[i].reqSend);
            MPI_Isend(cellCentres_z, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);

        }

        for ( i=0; i<totElmerRanks; i++ ) {
            ELp[i].nFoundCells = 0; // keep this
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
            MPI_Irecv(&ELp[i].nFoundCells, 1, MPI_INT, ELp[i].globalRank, 995,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        int totCellsFound = 0;
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            totCellsFound += ELp[i].nFoundCells;
        }

        if (totCellsFound < nCells) {
            Pout << "OpenFOAM2Elmer Warning: "<< "OpenFOAM #" << myLocalRank << " has " << nCells
                                 << " cells, Elmer found " << totCellsFound << endl;
        }

        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            ELp[i].foundCellsIndx = new int[ELp[i].nFoundCells];
            ELp[i].recvBuffer0 = new double[ELp[i].nFoundCells];

            if (ELp[i].foundCellsIndx == nullptr || ELp[i].recvBuffer0 == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }

            MPI_Irecv(ELp[i].foundCellsIndx, ELp[i].nFoundCells, MPI_INT, ELp[i].globalRank, 994,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        for ( i=0; i<totElmerRanks; i++ ) { 
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
        }
    }

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Sending fields to Elmer
    if (mode_==1) {

        // Extract the dictionary from the database
        const dictionary& fvSchemes = mesh_.template lookupObject<IOdictionary>
        (
           "fvSchemes"
        );

        // Exctract subdictionary from the main dictionary
        interpolationDict = fvSchemes.subDict("interpolationSchemes");
    
        nCommElem = 0;

        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Irecv(&ELp[i].nElem, 1, MPI_INT, ELp[i].globalRank, 899, MPI_COMM_WORLD, &ELp[i].reqRecv);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);

            nCommElem += ELp[i].nElem;

            ELp[i].sendBuffer0 = new double[ELp[i].nElem];
            ELp[i].sendBuffer1 = new double[ELp[i].nElem];
            ELp[i].sendBuffer2 = new double[ELp[i].nElem];
            ELp[i].foundElement = new label[ELp[i].nElem];

            if (ELp[i].sendBuffer0 == nullptr || ELp[i].sendBuffer1 == nullptr || 
                ELp[i].sendBuffer2 == nullptr || ELp[i].foundElement == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }

            MPI_Irecv(ELp[i].sendBuffer0, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
            MPI_Request_free(&ELp[i].reqRecv);
            MPI_Irecv(ELp[i].sendBuffer1, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
            MPI_Request_free(&ELp[i].reqRecv);
            MPI_Irecv(ELp[i].sendBuffer2, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        if (nCommElem ==  0 ) {
            FatalErrorInFunction 
            //Info << "OpenFOAM2Elmer Warning: "
            << "OpenFOAM2Elmer: Mesh of size 0 received for interpolation!\n"
			<< "OpenFOAM2Elmer: No defined mesh counterpart for OF mesh in Elmer!" 
			//<< endl;
			<< Foam::abort(FatalError);
        } 

		Info<< "OpenFOAM2Elmer: Searching for cells.." << endl;

	for ( i=0; i<totElmerRanks; i++ )
	{
		ELp[i].nFoundElements = 0; // keep this
		if (!ELp[i].boxOverlap) continue;
		
		 MPI_Test_Sleep(ELp[i].reqRecv);
		/**Write or read from file**/

            string fname = "O2E-" + std::to_string(myLocalRank) + "-" + std::to_string(i) + ".out";
            std::ifstream f(fname);

            if( f.good() ) {
                Pout<< "O2E pair file exist.." << endl;
                for ( j=0; j<ELp[i].nElem; j++ ) {
                    f >> ELp[i].foundElement[j];
                    if (ELp[i].foundElement[j] > -1) ELp[i].nFoundElements++;
                }
            }
            else {
                	Pout<< "O2E pair file does not exist, creating.." << endl;

                	f.close();
                	std::ofstream f(fname);
					findCells();
		for ( j=0; j<ELp[i].nElem; j++ ) {
			f << ELp[i].foundElement[j] << std::endl;
		}
		}
            	f.close();
		
       	Pout<<  "Found " << ELp[i].nFoundElements << " out of " << ELp[i].nElem << " elements from Elmer #" << i << endl;
		MPI_Isend(&ELp[i].nFoundElements, 1, MPI_INT, ELp[i].globalRank, 895,
		MPI_COMM_WORLD, &ELp[i].reqSend);
	}
		
        Pout << "Search complete!" << endl;

        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
        

            if (ELp[i].nFoundElements == 0) continue;
            ELp[i].foundElementIndx = new int[ELp[i].nFoundElements];
            ELp[i].foundElementCellIndx = new int[ELp[i].nFoundElements];
            ELp[i].positions = new point[ELp[i].nFoundElements];
        

            if (ELp[i].foundElementIndx == nullptr ||
                ELp[i].foundElementCellIndx == nullptr ||
                ELp[i].positions == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }
            
            k = 0;
            for ( j=0; j<ELp[i].nElem; j++ ) {
                if (ELp[i].foundElement[j] > -1) {
                    ELp[i].foundElementIndx[k] = j;
                    ELp[i].foundElementCellIndx[k] = ELp[i].foundElement[j];
                    ELp[i].positions[k].x() = ELp[i].sendBuffer0[j];
                    ELp[i].positions[k].y() = ELp[i].sendBuffer1[j];
                    ELp[i].positions[k].z() = ELp[i].sendBuffer2[j];
                    k++;
                }
            }
        
            MPI_Isend(ELp[i].foundElementIndx, ELp[i].nFoundElements, MPI_INT, ELp[i].globalRank, 894,
                      MPI_COMM_WORLD, &ELp[i].reqSend);
        
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (ELp[i].nFoundElements == 0) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
        }
        
        Pout << "OpenFOAM2Elmer: Calculating field interpolation prerequisites.." << endl;
		interpPrepare();
        Pout << "Calculation of prerequisites complete!" << endl;

        MPI_Barrier(PstreamGlobals::MPI_COMM_FOAM);
        Info<< "OpenFOAM2Elmer Init = " << MPI_Wtime()-commTime << " s" << nl << endl;
    }
}


template <class meshT>
void Foam::Elmer<meshT>::recvScalar(volScalarField& field)
{
    int i, j;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells == 0 ) continue;
        MPI_Irecv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 1000,
                  MPI_COMM_WORLD, &ELp[i].reqRecv);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqRecv);
        for (j=0; j<ELp[i].nFoundCells; j++ ) {
	    list_s[ELp[i].foundCellsIndx[j]] = ELp[i].recvBuffer0[j];
        }
    }	
	field.internalField() = scalargpuList(Ulist_s);
}


template <class meshT>
void Foam::Elmer<meshT>::recvVector(volVectorField& field)
{
    int i, j, dim;

    for (dim=0; dim<3; dim++) { 
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Irecv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 1000,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            for (j=0; j<ELp[i].nFoundCells; j++ ) {
		list_v[ELp[i].foundCellsIndx[j]][dim] = ELp[i].recvBuffer0[j];
            }
        }
    }
    
	field.internalField() = vectorgpuList(Ulist_v);
}


template <class meshT>
void Foam::Elmer<meshT>::sendScalar(volScalarField& field)
{
	int i,j;
	double commTime = MPI_Wtime();
	scalarField iFs( field.internalField() );
	
    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;

        for (j=0; j<ELp[i].nFoundElements; j++) {
	double yc = iFs[ ELp[i].foundElementCellIndx[j] ];
	
	if( word(interpolationDict.lookup(field.name())) == "cellPoint") { 
	
		//Info<< "OpenFOAM2Elmer: interpolation = ON" << endl;
		
			double y0 = iFs[ ELp[i].interpCellIndx[j][0] ];
			double y1 = iFs[ ELp[i].interpCellIndx[j][1] ];
			double y2 = iFs[ ELp[i].interpCellIndx[j][2] ];
			
			yc = yc + ( (y0 + ( ( y1 + (y2-y1)*ELp[i].interpCoeff[j][0] ) - y0)*ELp[i].interpCoeff[j][1] ) - yc)*ELp[i].interpCoeff[j][2];
	}
            	ELp[i].sendBuffer0[j] = yc;//interpolateScalar(field, i, j);
	}
        MPI_Isend(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                  MPI_COMM_WORLD, &ELp[i].reqSend);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqSend);
    }

    Info<< "OpenFOAM2Elmer: send scalar = " << MPI_Wtime()-commTime << " s" << endl;
}

template <class meshT>
void Foam::Elmer<meshT>::sendVector(volVectorField& field)
{
	int i,j;
	double commTime = MPI_Wtime();
	vectorField iFv( field.internalField() );

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;

        for (j=0; j<ELp[i].nFoundElements; j++) {
			point yc_v( iFv[ ELp[i].foundElementCellIndx[j] ] );

			if( word(interpolationDict.lookup(field.name())) == "cellPoint") {
		
				//Info<< "OpenFOAM2Elmer: interpolation = ON" << endl;
				point 		y0_v( iFv[ ELp[i].interpCellIndx[j][0] ] );
				point 		y1_v( iFv[ ELp[i].interpCellIndx[j][1] ] );
				point 		y2_v( iFv[ ELp[i].interpCellIndx[j][2] ] );
		
				yc_v = yc_v + ( (y0_v + ( ( y1_v + (y2_v-y1_v)*ELp[i].interpCoeff[j][0] ) - y0_v)*ELp[i].interpCoeff[j][1] ) - yc_v)*ELp[i].interpCoeff[j][2]; 
			}
			ELp[i].sendBuffer0[j] = yc_v.x();
			ELp[i].sendBuffer1[j] = yc_v.y();
			ELp[i].sendBuffer2[j] = yc_v.z();
		}

        MPI_Isend(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
        MPI_Request_free(&ELp[i].reqSend);
        MPI_Isend(ELp[i].sendBuffer1, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
        MPI_Request_free(&ELp[i].reqSend);
        MPI_Isend(ELp[i].sendBuffer2, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqSend);
    }

    Info<< "OpenFOAM2Elmer: send vector = " << MPI_Wtime()-commTime << " s" << endl;
}

template <class meshT>
void Foam::Elmer<meshT>::findCells()
{
	Info<< "OpenFOAM2Elmer: Searching for cells.." << endl;
	int i,j;
		
	/*******************************Find cell bounds for each cell in OF mesh*******************************/
	int k,n;
	double minx, miny, minz, maxx, maxy, maxz;
	double x_, y_, z_;
	
	int cellIndxSize = 0;

	labelListList cellNeighbours = mesh_.cellCells();
	double **cellBounds = new double*[nCells];

	//Info<< "OpenFOAM2Elmer: Preparing mesh for interpolation" << endl;
	for (k=0; k<nCells; k++ )
	{
		cellBounds[k] = new double[6];
		if (cellBounds[k] == nullptr) {
			FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
		} 
		cellIndxSize = cellNeighbours[k].size();
		if(cellIndxSize < 1) continue;

		minx = cellCentres_x[k];
		miny = cellCentres_y[k];
		minz = cellCentres_z[k];

		maxx = minx;
		maxy = miny;
		maxz = minz;

		for (n=0; n<cellIndxSize; n++ )
		{
			x_ = cellCentres_x[ cellNeighbours[k][n] ];
			y_ = cellCentres_y[ cellNeighbours[k][n] ];
			z_ = cellCentres_z[ cellNeighbours[k][n] ];

			minx = min(minx,x_);
			miny = min(miny,y_);
			minz = min(minz,z_);

			maxx = max(maxx,x_);
			maxy = max(maxy,y_);
			maxz = max(maxz,z_);
		}
		
		if (minx == cellCentres_x[k]) minx = 2*cellCentres_x[k] - maxx;
		if (miny == cellCentres_y[k]) miny = 2*cellCentres_y[k] - maxy;
		if (minz == cellCentres_z[k]) minz = 2*cellCentres_z[k] - maxz;

		if (maxx == cellCentres_x[k]) maxx = 2*cellCentres_x[k] - minx;
		if (maxy == cellCentres_y[k]) maxy = 2*cellCentres_y[k] - miny;
		if (maxz == cellCentres_z[k]) maxz = 2*cellCentres_z[k] - minz;

		maxx = max(maxx,2*cellCentres_x[k] - minx);
		maxy = max(maxy,2*cellCentres_y[k] - miny);
		maxz = max(maxz,2*cellCentres_z[k] - minz);

		minx = min(minx,2*cellCentres_x[k] - maxx);
		miny = min(miny,2*cellCentres_y[k] - maxy);
		minz = min(minz,2*cellCentres_z[k] - maxz);

		cellBounds[k][0] = minx;
		cellBounds[k][1] = maxx;

		cellBounds[k][2] = miny;
		cellBounds[k][3] = maxy;
		
		cellBounds[k][4] = minz;
		cellBounds[k][5] = maxz;
		
	}
	
	/*******************************Find closest OF cell for each point in elmer mesh*******************************/

	double oldDist;
	double newDist;

	int nElemDone = 0;
	int nElemDonePrev = 0;
	for ( i=0; i<totElmerRanks; i++ )
	{
        if (!ELp[i].boxOverlap) continue;
		const label ELsize = ELp[i].nElem;
		vectorList rawElmerPoints;
		rawElmerPoints.setSize(ELsize);
		for (j=0; j<ELp[i].nElem; j++ )
		{
			rawElmerPoints[j] = vector(	ELp[i].sendBuffer0[j],
							ELp[i].sendBuffer1[j],
							ELp[i].sendBuffer2[j]);
		}

		const vectorUList& constPtr = rawElmerPoints;
		vectorgpuList rawElmerGPUPoints(constPtr);
		ELp[i].foundElmerGPUIndx.setSize(ELsize,-1);
	
        localTime = MPI_Wtime();
		for (k=0; k<nCells; k++ )
		{
         	nElemDone++;
         	if(MPI_Wtime()-localTime>30 || nElemDone==(nCells*totElmerRanks))
			{
             	Pout 	<< 100.0*nElemDone/(nCells*totElmerRanks) << "% done, search speed "
            			<< 2*(nElemDone-nElemDonePrev) << " points/min, remaining "
            			<< (nCells*totElmerRanks)-nElemDone << " points" << endl;
            	localTime = MPI_Wtime();
            	nElemDonePrev = nElemDone;			
			}

			const double minX = cellBounds[k][0];
			const double maxX = cellBounds[k][1];

			const double minY = cellBounds[k][2];
			const double maxY = cellBounds[k][3];

			const double minZ = cellBounds[k][4];
			const double maxZ = cellBounds[k][5];

			const int cellI = k;

			thrust::transform
			(
				//iterator and first input variable
				rawElmerGPUPoints.begin(),
				rawElmerGPUPoints.end(),
				//second input variable
				ELp[i].foundElmerGPUIndx.begin(),
				//return variable
				ELp[i].foundElmerGPUIndx.begin(),
				//functor call with input parameters
				findCellFunctor(minX,maxX,minY,maxY,minZ,maxZ,cellI)
			);
		}

		labelField foundElmerIndx(ELp[i].foundElmerGPUIndx);

		for ( j=0; j<ELp[i].nElem; j++ )
		{
			ELp[i].foundElement[j] = foundElmerIndx[j];
			//use the values found with findCellFunctor instead of very slow built in algorithm
			//ELp[i].foundElement[j] = mesh_.findCell(tmpPoint,polyMesh::FACEDIAGTETS);
			if (ELp[i].foundElement[j] > -1)
			{
				ELp[i].nFoundElements++;
                point tmpPoint(ELp[i].sendBuffer0[j],ELp[i].sendBuffer1[j],ELp[i].sendBuffer2[j]);
				k = ELp[i].foundElement[j];

				//Distance to cell center
				point tmpDist(	tmpPoint[0]-cellCentres_x[k],
								tmpPoint[1]-cellCentres_y[k],
								tmpPoint[2]-cellCentres_z[k]);
				oldDist = (tmpDist & tmpDist);

				for (n=0; n<cellNeighbours[k].size(); n++ )
				{
					//Prioritize neighbour cell if distance  is smaller
                	point tmpDist(	tmpPoint[0]-cellCentres_x[cellNeighbours[k][n]],
									tmpPoint[1]-cellCentres_y[cellNeighbours[k][n]],
									tmpPoint[2]-cellCentres_z[cellNeighbours[k][n]]);
			
					newDist = (tmpDist & tmpDist);
					if(newDist < oldDist )
					{
						oldDist = newDist;
						ELp[i].foundElement[j] = cellNeighbours[k][n];
					}
				}
			}
			else ELp[i].foundElement[j] = -1;
		}
	}
}


template <class meshT>
void Foam::Elmer<meshT>::interpPrepare()
{
	/******Finds neighbour cells and coefficients for scalar and vector field interpolation*******/
	
	int i,j,n,n1,n2;
	labelListList cellNeighbours = mesh_.cellCells();

	for ( i=0; i<totElmerRanks; i++ )
	{
        if ( ELp[i].nFoundElements == 0 ) continue;
        
		int elmerElements = ELp[i].nFoundElements;
		int interpCounter = 0;
		
		ELp[i].interpCellIndx = new int *[elmerElements];
		ELp[i].interpCoeff = new point[elmerElements];
        if (ELp[i].interpCellIndx == nullptr || ELp[i].interpCoeff == nullptr)
		{
			FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
		} 

        for (j=0; j<elmerElements; j++)
		{
			int thisCellIndx = ELp[i].foundElementCellIndx[j];
			int cellIndxSize = cellNeighbours[thisCellIndx].size();

			ELp[i].interpCellIndx[j] = new int[3];
        	if (ELp[i].interpCellIndx[j] == nullptr)
			{
        		FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
			} 
			ELp[i].interpCellIndx[j][0] = thisCellIndx;
			ELp[i].interpCellIndx[j][1] = thisCellIndx;
			ELp[i].interpCellIndx[j][2] = thisCellIndx;
			ELp[i].interpCoeff[j] = point(0,0,0);

			if(cellIndxSize < 3) continue;//Only for 2D. Not neccessary???

			int* cellNeighboursIndx = new int[cellIndxSize];
			double* cellNeighboursDist = new double[cellIndxSize];
        	if (cellNeighboursIndx == nullptr || cellNeighboursDist == nullptr)
			{
        		FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
        	}

			point tmpPoint = ELp[i].positions[j];
			//Prioritize neighbour cell if distance  is smaller
			for (n=0; n<cellIndxSize; n++ )
			{
				//Prioritize neighbour cell if distance  is smaller
				point tmpDist(	tmpPoint[0]-cellCentres_x[cellNeighbours[thisCellIndx][n]],
								tmpPoint[1]-cellCentres_y[cellNeighbours[thisCellIndx][n]],
								tmpPoint[2]-cellCentres_z[cellNeighbours[thisCellIndx][n]]);

				cellNeighboursDist[n] = mag(tmpDist);
				cellNeighboursIndx[n] = n;
			}

			for (n=0; n<cellIndxSize; n++ )
			{
				for (n1=n+1; n1<cellIndxSize; n1++ )
				{
					if(cellNeighboursDist[n1] < cellNeighboursDist[n] )
					{
						//if smaller distance switch places
						double oldDist = cellNeighboursDist[n];
						int oldIndx = cellNeighboursIndx[n];
						cellNeighboursDist[n] = cellNeighboursDist[n1];
						cellNeighboursIndx[n] = cellNeighboursIndx[n1];

						cellNeighboursDist[n1] = oldDist;
						cellNeighboursIndx[n1] = oldIndx;
					}
				}
			}

			/*********************************************************************************
			To minimize the likelyhood of having multiple cell centers on the same line we can 
			try to find three neighbour cells that are themselves neighbours.
			
			To do so we start by keeping the closest neighbour cell 
			and try to find two more cells with common neighbours
			*********************************************************************************/
			
			int firstNeighbour = cellNeighbours[thisCellIndx][cellNeighboursIndx[0]];
			int secondNeighbour = cellNeighbours[thisCellIndx][cellNeighboursIndx[1]];
			int thirdNeighbour = cellNeighbours[thisCellIndx][cellNeighboursIndx[2]];
			bool secondNeighbourGood = false;
			bool thirdNeighbourGood = false;

			for (n=0; n<cellIndxSize-2; n++ )
			{
				if(!secondNeighbourGood)
				{
					secondNeighbour = cellNeighbours[thisCellIndx][cellNeighboursIndx[n+1]];
					for (n1=0; n1<cellNeighbours[secondNeighbour].size(); n1++ )
					{
						if(cellNeighbours[secondNeighbour][n1] == thisCellIndx) continue;

						for (n2=0; n2<cellNeighbours[firstNeighbour].size(); n2++ )
						{
							if(cellNeighbours[secondNeighbour][n1]==cellNeighbours[firstNeighbour][n2] )
							{
								secondNeighbourGood = true;
								break;
							}
						}
						if(secondNeighbourGood) break;
					}
				}

				if(secondNeighbourGood)
				{
					bool commonFirstNeighhbour = false;
					bool commonSecondNeighbour = false;
					thirdNeighbour = cellNeighbours[thisCellIndx][cellNeighboursIndx[n+2]];
					for (n1=0; n1<cellNeighbours[thirdNeighbour].size(); n1++ )
					{
						if( cellNeighbours[thirdNeighbour][n1] == thisCellIndx) continue;

						if(!commonFirstNeighhbour) for (n2=0; n1<cellNeighbours[firstNeighbour].size(); n2++ )
						{
							if(cellNeighbours[thirdNeighbour][n1]==cellNeighbours[firstNeighbour][n2])
							{
								commonFirstNeighhbour = true;
								break;
							}
						}
						//third connection does not need to be through the same neighbour
						if(!commonSecondNeighbour) for (n2=0; n1<cellNeighbours[secondNeighbour].size(); n2++ )
						{
							if(cellNeighbours[thirdNeighbour][n1]==cellNeighbours[secondNeighbour][n2])
							{
								commonSecondNeighbour = true;
								break;
							}
						}
						//if this cell has neighbours with both previous cells then break
						if(commonFirstNeighhbour&&commonSecondNeighbour)
						{
							thirdNeighbourGood = true;
							break;
						}
					}

					if(thirdNeighbourGood)
					{
						bool done = false;
						int tempCells[3] = {firstNeighbour,secondNeighbour,thirdNeighbour};

        				vector la(	cellCentres_x[thisCellIndx],
									cellCentres_y[thisCellIndx],
									cellCentres_z[thisCellIndx]);
                		vector lab = tmpPoint-la;

						for (n1=0; n1<3; n1++ )
						{
                			vector p0(	cellCentres_x[tempCells[(0+n1)%3]],
										cellCentres_y[tempCells[(0+n1)%3]],
										cellCentres_z[tempCells[(0+n1)%3]]);

                			vector p1(	cellCentres_x[tempCells[(1+n1)%3]],
										cellCentres_y[tempCells[(1+n1)%3]],
										cellCentres_z[tempCells[(1+n1)%3]]);

                			vector p2(	cellCentres_x[tempCells[(2+n1)%3]],
										cellCentres_y[tempCells[(2+n1)%3]],
										cellCentres_z[tempCells[(2+n1)%3]]);

							// if vectors (p1-p0) and (p2-p0) happen to be on the same line then skip interpolation
							if(abs(((p1-p0) & (p2-p0))/mag(p1-p0)/mag(p2-p0))==1) break;
							vector p12 = (p1-p0) ^ (p2-p0);
							double t = -( p12 & (la-p0) )/( lab & p12 );
							vector px = la + t*lab;

							// if vectors (la-p0) and (px-p0) happen to be on the same line then shuffle points
							if(abs(abs(((la-p0) & (px-p0))/mag(la-p0)/mag(px-p0)) - 1.0) < 1e-6 ) continue;

							//if succesfully found cells for interpolation then prepare interpolation coeffs
							vector p12_ = (la-p0) ^ (px-p0);
							double t_ = -( p12_ & (p1-p0) )/( (p2-p1) & p12_ );
							vector px_ = p1 + t_*(p2-p1);

							point tmpCoeffs(sign((px_-p1)& (p2-p1) )*mag(px_-p1)/mag(p2-p1),
											sign((px-p0) & (px_-p0) )*mag(px-p0)/mag(px_-p0),
											sign((lab) & (px-la) )*mag(lab)/mag(px-la) );
											
							//next part tests the found interpolation coefficients
							//with some values: la=1, p0=2, p1=3, p2=4
							double testPoint = 1.0 + ( (2.0 + ( ( 3.0 + (4.0-3.0)*tmpCoeffs[0] ) - 2.0)*tmpCoeffs[1] ) - 1.0)*tmpCoeffs[2];
							//Arbitrary filter but if returned interpolated value many times exceeds the test value then probably something went wrong
							if(abs(testPoint) < 40.0 )
							{
								//if the interpolated value is not much larger than the input values then 
								//save coefficients and cell numbers for interpolation
								ELp[i].interpCoeff[j] = tmpCoeffs;

								ELp[i].interpCellIndx[j][0] = tempCells[(0+n1)%3];
								ELp[i].interpCellIndx[j][1] = tempCells[(1+n1)%3];
								ELp[i].interpCellIndx[j][2] = tempCells[(2+n1)%3];

								interpCounter++;
								done = true;
								break;
							}
						}
						//when done return to elmer element loop
						if(done) break;
						else
						{
							//if failed to find neighbour cells with good interpolation conditions then continue search
							thirdNeighbourGood = false;
						}
					}
					else
					{
						//if this cell has neighbours only with one or none of previous cells then continue search
						thirdNeighbourGood = false;
					}
				}
			}
			
		}
		if(interpCounter != elmerElements)
		{
		Pout 	<< "Couldn't find favorable condition for interpolation and will use closest cell values for "  << elmerElements-interpCounter 
			<< " or " << (elmerElements-interpCounter)*100.0/elmerElements 
			<< "% out of " << elmerElements  << " elmer elements of rank " << i 
			<< endl;
		}
	}
	
}


template <class meshT>
void Foam::Elmer<meshT>::sendStatus(int status)
{
    int i;

    if (myLocalRank==0) {
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Isend(&status, 1, MPI_INT, ELp[i].globalRank,799, MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Test_Sleep(ELp[i].reqSend);
        }
    }
}


template <class meshT>
void Foam::Elmer<meshT>::findOverlappingBoxes()
{
    int OFoffset = myLocalRank*totElmerRanks;

    if ( myLocalRank==0 ) {
        MPI_Recv(ELboundBoxes, totElmerRanks*2*3, MPI_DOUBLE, ELp[0].globalRank, 1001,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(ELboundBoxes, totElmerRanks*2*3, MPI_DOUBLE, 0, PstreamGlobals::MPI_COMM_FOAM);

    for (int i=0; i<totElmerRanks; i++ ) {
        point tmpPointMin;
        point tmpPointMax;

        int ELoffset = i*2*3;
        tmpPointMin.x() = ELboundBoxes[ELoffset];
        tmpPointMin.y() = ELboundBoxes[ELoffset+1];
        tmpPointMin.z() = ELboundBoxes[ELoffset+2];
        tmpPointMax.x() = ELboundBoxes[ELoffset+3];
        tmpPointMax.y() = ELboundBoxes[ELoffset+4];
        tmpPointMax.z() = ELboundBoxes[ELoffset+5];
        boundBox tmpBoundBox(tmpPointMin, tmpPointMax);

        ELp[i].boxOverlap = myBoundBox.overlaps(tmpBoundBox);
              
        OF_EL_overlap[OFoffset+i] = int(ELp[i].boxOverlap);
    }

    MPI_Allgather(&OF_EL_overlap[OFoffset], totElmerRanks, MPI_INT,
        OF_EL_overlap, totElmerRanks, MPI_INT, PstreamGlobals::MPI_COMM_FOAM);

    if ( myLocalRank==0 ) {
        MPI_Send(OF_EL_overlap, totLocalRanks*totElmerRanks, MPI_INT, ELp[0].globalRank, 1002,
                 MPI_COMM_WORLD);
    }
}


template <class meshT>
void Foam::Elmer<meshT>::MPI_Test_Sleep(MPI_Request& req)
{
    int flag;

    while ( true ) {
        MPI_Test( &req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        nanosleep((const struct timespec[]){{0, 1000000L}}, NULL);
    }
}

// explicit instantiations
template class Foam::Elmer<fvMesh>;

// ************************************************************************* //
