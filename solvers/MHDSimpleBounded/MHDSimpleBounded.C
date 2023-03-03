/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Foam::Elmer<fvMesh> sending(mesh,1, minPoint, maxPoint); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    // Receive fields from Elmer
    Foam::Elmer<fvMesh> receiving(mesh,-1, minPoint, maxPoint); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JxB);

    double OFstep = 0;
    double Elmerstep = 0;
    while (simple.loop())
    {
	OFstep = runTime.elapsedClockTime();
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();
	
        runTime.write();
	OFstep = runTime.elapsedClockTime()-OFstep;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Elmerstep = runTime.elapsedClockTime();
	
        dimensionedScalar smallU
        (
            "smallU",
            dimensionSet(0, 1, -1, 0, 0, 0 ,0),
            1e-6
        );

        // Check whether we need to update electromagnetic stuff with Elmer
        scalar maxRelDiff_local = (max(mag(U_old-U)/(average(mag(U))+smallU))).value();

        bool doElmer = false;
        if(maxRelDiff_local>maxRelDiff && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0)) {
            doElmer = true;
        }
        if(doElmer && runTime.run()) {
            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(1);
            sending.sendVector(U);

            // Receive fields form Elmer
            receiving.sendStatus(1);
            receiving.recvVector(JxB);
        }
	Elmerstep=runTime.elapsedClockTime()-Elmerstep;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  OF step = " << OFstep << " s"
            << "  Elmer step = " << Elmerstep << " s" 
            //<< "  err = ( " << maxErr << " )" 
		<< nl << endl;

    }
        //Final iter for Elmer
        U_old = U;
        // Send fields to Elmer
        sending.sendStatus(0);
        sending.sendVector(U);
        // Receive fields form Elmer
        receiving.sendStatus(0);
        receiving.recvVector(JxB);

        Info<< "Final iteration: "
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s = " 
	<< runTime.elapsedClockTime()/60.0 << " min = " 
	<< runTime.elapsedClockTime()/3600.0 << " h = " 
	<< runTime.elapsedClockTime()/3600.0/24.0 << " days"
	<< nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
