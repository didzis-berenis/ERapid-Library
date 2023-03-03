/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    pimpleFoam

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable finite volume options, e.g. MRF, explicit porosity
  
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H" 
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
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

    pimpleControl pimple(mesh); 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    double Elmerstep = runTime.elapsedClockTime(); 

    // Send fields to Elmer
    //Foam::Elmer<fvMesh> sending(mesh,1, minPoint, maxPoint); // 1=send, -1=receive
    Foam::Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    sending.sendVector(U);

    // Receive fields from Elmer
    //Foam::Elmer<fvMesh> receiving(mesh,-1, minPoint, maxPoint); // 1=send, -1=receive
    Foam::Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(J);
    receiving.recvVector(B);
    JxB = J ^ B;

	Elmerstep=runTime.elapsedClockTime()-Elmerstep;
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< "  Elmer init = " << Elmerstep << " s" 
		<< endl;
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
	double OFstep = runTime.elapsedClockTime();
    if (runTime.run()) 
    {
		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"
		scalar writeLastSteps( readScalar(runTime.controlDict().lookup("writeLastSteps")) );

    	Info<< "\nStarting time loop\n" << endl;         

    	runTime++;

    	Info<< "\nTime = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())    
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
		if(!runTime.outputTime() && runTime.time().value() > runTime.endTime().value()-runTime.deltaT().value()*(floor(writeLastSteps) - 0.5) )
		{	runTime.writeNow();	}
	}
	OFstep = runTime.elapsedClockTime()-OFstep;
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    while (runTime.run()) 
    {
		Elmerstep = runTime.elapsedClockTime(); 
        // Check whether we need to update electromagnetic stuff with Elmer
        scalar maxRelDiff_local = (max(mag(U_old-U)/(average(mag(U))+smallU))).value();
        if(maxRelDiff_local>maxRelDiff && (maxRelDiff<SMALL || maxRelDiff+SMALL<=1.0))
		{
            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(runTime.run());
            sending.sendVector(U);

            // Receive fields form Elmer
            receiving.sendStatus(runTime.run());
    	    receiving.recvVector(J);
    	    receiving.recvVector(B);
	    JxB = J ^ B;
        }
		Elmerstep=runTime.elapsedClockTime()-Elmerstep;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  OF step = " << OFstep << " s"
            << "  Elmer step = " << Elmerstep << " s" 
			<< endl;
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
		OFstep = runTime.elapsedClockTime();

		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"
		scalar writeLastSteps( readScalar(runTime.controlDict().lookup("writeLastSteps")) );
		p.storePrevIter();
		U.storePrevIter();       

        runTime++;

        Info<< "\nTime = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())    
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
		if(!runTime.outputTime() && runTime.time().value() > runTime.endTime().value()-runTime.deltaT().value()*(floor(writeLastSteps) - 0.5) )
		{	runTime.writeNow();	}

		OFstep = runTime.elapsedClockTime()-OFstep;
    }
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
	//Info<< "\nHanding final iteration to elmer\n" << endl; 
	Elmerstep = runTime.elapsedClockTime();  
	// Send fields to Elmer
	sending.sendStatus(runTime.run());
	sending.sendVector(U);

	// Receive fields form Elmer
	receiving.sendStatus(runTime.run());
	receiving.recvVector(J);
	receiving.recvVector(B);
	//JxB = J ^ B;
	Elmerstep=runTime.elapsedClockTime()-Elmerstep;

	Info<< "Final iteration: "
		<< "  OF step = " << OFstep << " s"
		<< "  Elmer step = " << Elmerstep << " s" 
		<< nl << endl;

	Info<< "Calculation complete: "
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
