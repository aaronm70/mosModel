#include<iostream>
#include"lModH.h"
#include <iostream>
#include <map>
#include <random>
#include <chrono>
#include <fstream>



int main()
{


	modParms initParms;
	obsDatX garkiDatGam;
	// read in garki village data anoph gamb values for each time point (time points are in continuous days and linked to corresponding time points in the rainfall data)
	garkiDatGam.garki154 = { {14,1},{28,9},{42,8},{56,22},{70,27},{84,22},{98,104},{112,52},{140,26},{154,18},{186,4},{215,0},{242,0},{270,0},{298,1} };

	garkiDatGam.garki408 = { {4,2},{18,45},{32,382},{46,272},{60,486},{74,645},{88,261},{102,288},{116,38},{130,20},{159,13},{187,3},{215,1},{243,1},{271,1} };

	garkiDatGam.garki553 = { {11,12},{25,29},{39,59},{53,58},{67,111},{81,152},{95,60},{109,9},{123,4},{150,6},{178,0},{206,0} };//, { 220,0 }, { 234,0 }, { 262,0 }, { 290,0 }, { 318,1 }, { 346,8 }, { 361,17 }, { 375,34 }, { 389,23 }, { 403,7 }, { 417,7 }, { 433,50 }, { 447,25 }, { 461,14 }, { 475,2 }, { 489,0 }, { 510,0 }, { 538,0 }, { 566,0 }, { 594,0 }, { 622,0 }, { 650,0 }, { 678,0 }, { 706,0 }, { 727,3 }, { 741,0 }, { 755,8 }, { 769,57 }, { 783,226 }, { 797,71 }, { 811,58 }, { 825,15 }};
	garkiDatGam.garki553_2 = { { 290,0 }, { 318,1 }, { 346,8 }, { 361,17 }, { 375,34 }, { 389,23 }, { 403,7 }, { 417,7 }, { 433,50 }, { 447,25 }, { 461,14 }, { 475,2 }, { 489,0 } };

	//MW vill 55
	garkiDatGam.garki55 = { {8,0},{22,1},{36,17},{50,39},{64,57},{78,137},{92,406},{106,126},{120,55},{134,41},{148,65},{172,56},{203,5},{231,1},{259,0},{287,0} };

	garkiDatGam.garki802 = { {3,6},{17,27},{31,66},{45,145},{59,9},{73,48},{87,23},{101,14},{115,13},{131,8},{163,0},{191,0} };//, { 247,0 }, { 275,0 }, { 303,1 }, { 331,0 }, { 367,43 }, { 381,6 }, { 395,4 }, { 409,7 }, { 423,12 }, { 439,13 }, { 453,8 }, { 467,2 }, { 481,0 }, { 495,0 }, { 523,0 }, { 551,0 }, { 579,0 }, { 607,0 } };
	garkiDatGam.garki802_2 = { { 303,1 }, { 331,0 }, { 367,43 }, { 381,6 }, { 395,4 }, { 409,7 }, { 423,12 }, { 439,13 }, { 453,8 }, { 467,2 }, { 481,0 } };


	//MW vill 202
	garkiDatGam.garki202 = { {7,1},{21,11},{35,7},{49,24},{63,97},{77,85},{91,54},{105,20},{119,36},{138,13},{169,8},{197,0},{225,0},{253,0},{280,0} };
	//MW vill 304
	garkiDatGam.garki801 = { {1,0},{15,0},{29,4},{43,4},{57,21},{71,69},{85,244},{99,197},{113,106},{141,14},{155,11},{189,3},{217,0},{273,0},{301,0} };
	//MW vill 218
	garkiDatGam.garki801_2 = { {4,3},{18,29},{32,38},{46,86},{60,63},{74,87},{88,103},{102,46},{116,40},{137,13},{169,2},{197,0},{225,0},{253,0} };


	//starting proposal distribution SD's for each parameter
		vector<double> sdProps = {
			0.001, 0.001, 0.01,0.01,1,0.0001,1, 
			1,1,1,1,1, 1,1,1,1,1,
			1,1,1,1,1, 1,1,1,1,1
			,0.01,0.05,0.05,5,1,0.5,0.001,1,1,1};

		//Maximum proposal distribution SD's
		vector<double> maxSdProps = {
			0.05, 0.05, 0.8, 0.5,6,0.01,
			5, //n
			5, 5, 5,5,5,5,5,5,5,5,//z1:10
			5,5,5, 5,5,5,5,5,5,5,//sf1:10
			0.1,0.1,0.1,5,1,1,0.01,0,0,0};

		//Acceptance ratios
		vector<double> acptRs = {
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25 ,0.25 };



		//Read in pMCMC options from text file to enable multiple instances to run on the cluster with differing setups
		pmcmcOptions pmcmcOpt;
		pmcmcOpt = optionsReader("C:\\lmodCPP\\x64\\powerClumped\\paramOptions.txt", pmcmcOpt);//read in pMCMC options
		string outputFolder = pmcmcOpt.outputFolder; //output folder for results
		string dFunc = pmcmcOpt.dFunc; //Which density/egg laying functions to use: "expClumped", "linearClumped","powerClumped","expNoClumped", "linearNoClumped" or "powerNoClumped"

		string initParamsLoc = pmcmcOpt.initParamsLoc; //location of initial parameter values - taken from previous pMCMC run and found in pMCMC options text file
		initParms = initParamsReader(initParamsLoc, initParms);//initial parameters, taken from previous runs

		vector<tuple<string, double>> fitPrms = { { "uoE", initParms.uoE },{ "uoL", initParms.uoL },{ "uP", initParms.uP },{ "uM", initParms.uM },{ "Y", initParms.Y },
		{ "w", initParms.w },{ "n", initParms.n },{ "z1", initParms.z1 },{ "z2", initParms.z2 },{ "z3", initParms.z3 },{ "z4", initParms.z4 },{ "z5", initParms.z5 },{ "z6", initParms.z6 },{ "z7", initParms.z7 },{ "z8", initParms.z8 },{ "z9", initParms.z9 },{ "z10", initParms.z10 },
		{ "sf1", initParms.sf1 } ,{ "sf2", initParms.sf2 } ,{ "sf3", initParms.sf3 }   ,{ "sf4", initParms.sf4 } ,{ "sf5", initParms.sf5 } ,{ "sf6", initParms.sf6 } ,{ "sf7", initParms.sf7 } ,{ "sf8", initParms.sf8 } ,{ "sf9", initParms.sf9 } ,{ "sf10", initParms.sf10 }
		,{ "dE", initParms.dE },{ "dL", initParms.dL } ,{ "dP", initParms.dP } ,{ "o", initParms.o },{ "tau", initParms.tau },{ "Mg", initParms.Mg },{"p",initParms.p },{ "lK",initParms.lK } ,{ "lKs",initParms.lKs } ,{ "lKm",initParms.lKm } };


		//Main function for running model fitting
		pMMHres results = pMMHSampler(
			initParms,//initial parameters
			dFunc,//density function to use: "expClumped", "linearClumped","powerClumped","expNoClumped", "linearNoClumped" or "powerNoClumped"
			sdProps,//initial sd for param proposals
			acptRs,//acceptance ratios
			fitPrms,//tuple of initial parm values plus names - needed as no reflection...
			maxSdProps,//max sd for each parameter proposal in tuner
			50000,//iterations
			25,//particles
			pmcmcOpt.nburn,//nburn 
			pmcmcOpt.monitoring,//monitoring
			pmcmcOpt.startAdapt,//start adapt
			pmcmcOpt.tell,//tell
			garkiDatGam//observed data
		);

		//write results to csv
		resultsWriter("C:\\ImperialOld\\Imperial\\lModCpp\\Data\\meds",outputFolder, results);
		string fileNameFit = "C:\\ImperialOld\\Imperial\\lModCpp\\Data\\meds";
		fileNameFit.append(outputFolder);
		pFitFunc(250, results, garkiDatGam, initParms, fileNameFit, dFunc);
		cout << "End";

	cin.get();
}







