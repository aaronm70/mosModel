#include"lModH.h"
#include <algorithm>
boost::mt19937 rng(std::time(0));
double inf = std::numeric_limits<double>::infinity();


vector<double> rainfall_05 = txtReader("C:\\lmodCPP\\data\\rf05.txt", 0.25);
vector<double> rainfall_07 = txtReader("C:\\lmodCPP\\data\\rf07.txt", 0.25);
vector<double> rainfall_08 = txtReader("C:\\lmodCPP\\data\\rf08.txt", 0.25);
vector<double> rainfall_04 = txtReader("C:\\lmodCPP\\data\\rf04.txt", 0.25);
vector<double> rainfall_02 = txtReader("C:\\lmodCPP\\data\\rf02.txt", 0.25);
vector<double> rainfall_01 = txtReader("C:\\lmodCPP\\data\\rf01.txt", 0.25);
vector<double> rainfall_03 = txtReader("C:\\lmodCPP\\data\\rf03.txt", 0.25);
vector<double> rainfall_06 = txtReader("C:\\lmodCPP\\data\\rf06.txt", 0.25);


/*proposal sd tuning function
@param current s.d.
@param target acceptance ratio
@param current acceptance ratio
@param maximum proposal s.d.
@return tuned s.d.*/
double tuner(double curSd, double acptR, double curAcptR, double maxSd){
	boost::math::normal dist(0.0, 1.0);
	if (curAcptR == 1)
		curAcptR = 0.99;
	if (curAcptR == 0)
		curAcptR = 0.01;
	double res = (curSd *quantile(dist, (acptR / 2))) / quantile(dist, (curAcptR / 2));
	if (res > maxSd) res = maxSd;
	if (res <= 0) res = 0;
	return res;
}


//random num function 0-1
double rn01(void)
{
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}


//structure updating function - could be improved?
modParms parmUpdt(modParms prms, string prmName, double propPrm) {
	if (prmName == "uoE")
		prms.uoE = propPrm;
	if (prmName == "uoL")
		prms.uoL = propPrm;
	if (prmName == "uP")
		prms.uP = propPrm;
	if (prmName == "Y")
		prms.Y = propPrm;
	if (prmName == "z1")
		prms.z1 = propPrm;
	if (prmName == "z2")
		prms.z2 = propPrm;
	if (prmName == "z3")
		prms.z3 = propPrm;
	if (prmName == "z4")
		prms.z4 = propPrm;
	if (prmName == "z5")
		prms.z5 = propPrm;
	if (prmName == "z6")
		prms.z6 = propPrm;
	if (prmName == "z7")
		prms.z7 = propPrm;
	if (prmName == "z8")
		prms.z8 = propPrm;
	if (prmName == "z9")
		prms.z9 = propPrm;
	if (prmName == "z10")
		prms.z10 = propPrm;
	if (prmName == "w")
		prms.w = propPrm;
	if (prmName == "sf1")
		prms.sf1 = propPrm;
	if (prmName == "sf2")
		prms.sf2 = propPrm;
	if (prmName == "sf3")
		prms.sf3 = propPrm;
	if (prmName == "sf4")
		prms.sf4 = propPrm;
	if (prmName == "sf5")
		prms.sf5 = propPrm;
	if (prmName == "sf6")
		prms.sf6 = propPrm;
	if (prmName == "sf7")
		prms.sf7 = propPrm;
	if (prmName == "sf8")
		prms.sf8 = propPrm;
	if (prmName == "sf9")
		prms.sf9 = propPrm;
	if (prmName == "sf10")
		prms.sf10 = propPrm;
	if (prmName == "n")
		prms.n = propPrm;
	if (prmName == "dE")
		prms.dE = propPrm;
	if (prmName == "dL")
		prms.dL = propPrm;
	if (prmName == "dP")
		prms.dP = propPrm;
	if (prmName == "o")
		prms.o = propPrm;
	if (prmName == "uM")
		prms.uM = propPrm;
	if (prmName == "Mg")
		prms.Mg = propPrm;
	if (prmName == "p")
		prms.p = propPrm;
	if (prmName == "tau")
		prms.tau = propPrm;
	if (prmName == "lK")
		prms.lK = propPrm;
	if (prmName == "lKs")
		prms.lKs = propPrm;
	if (prmName == "lKm")
		prms.lKm = propPrm;
	return prms;
}

/*log likelihood function
@param particles number of particles in pMMH
@param prms current parameters
@param obsDatX observed data
@param fixedParam fixed parameters
@return log likelihood value*/
double llFunc(int particles, modParms prms, obsDatX obsDat, string dFunc) {
	vector<double> pfiltRes;

	for (auto j = 0; j != 10; ++j) {

		vector<tuple<int, int>> oDat;
		if (j == 0) {
			oDat = obsDat.garki408;
			prms.sf = prms.sf1;
			prms.z = prms.z1;
			prms.rF = rainfall_03;
		}
		else if (j == 1) {
			oDat = obsDat.garki154;
			prms.sf = prms.sf2;
			prms.z = prms.z2;
			prms.rF = rainfall_05;
		}
		else if (j == 2) {
			oDat = obsDat.garki801;//MW304
			prms.sf = prms.sf3;
			prms.z = prms.z3;
			prms.rF = rainfall_08;
		}
		else if (j == 3) {
			oDat = obsDat.garki802;
			prms.sf = prms.sf4;
			prms.z = prms.z4;
			prms.rF = rainfall_01;
		}
		else if (j == 4) {
			oDat = obsDat.garki553;
			prms.sf = prms.sf5;
			prms.z = prms.z5;
			prms.rF = rainfall_02;
		}
		else if (j == 5) {
			oDat = obsDat.garki801_2;//MW 218
			prms.sf = prms.sf6;
			prms.z = prms.z6;
			prms.rF = rainfall_07;
		}
		else if (j == 6) {
			oDat = obsDat.garki802_2;
			prms.sf = prms.sf7;
			prms.z = prms.z7;
			prms.rF = rainfall_01;
		}
		else if (j == 7) {
			oDat = obsDat.garki553_2;
			prms.sf = prms.sf8;
			prms.z = prms.z8;
			prms.rF = rainfall_02;
		}

		else if (j == 8) {
			oDat = obsDat.garki202;//MW 202
			prms.sf =  prms.sf9;
			prms.z =  prms.z9;
			prms.rF = rainfall_04;
		}
		else if (j == 9) {
			oDat = obsDat.garki55;//MW 55
			prms.sf =  prms.sf10;
			prms.z =  prms.z10;
			prms.rF = rainfall_06;
		}

		//run particle filter
		prms.sf = pow(10,prms.sf);
		prms.z = pow(10,prms.z);
		pfiltRes.emplace_back(pFilt(particles,
			oDat,//garki data
			prms,//parameters
			false,//full output or just likelihood
			"C:\\lmodCPP\\data\\test.txt",
			false,
			dFunc//density function 
		));
	}

	return boost::accumulate(pfiltRes, 0.0);
}

/*Parameter proposal function
@param sd parameter standard deviation
@param parameter
@return proposed parameter*/
double propPrmFunc(double sd, double parm) {
	boost::normal_distribution<> nd(0.0, sd);
	double ran = nd(rng);
	double prop = ran + parm;
	if (prop <= 0) prop = 0;
	return prop;
}

/*log prior function
@param current parameters
@return sum loglikelihood for each parameter*/
double lprior(modParms prms, string dFunc) {
	double res = 0;
	boost::math::normal_distribution<double> d1(0.035, 0.0066);//uoE
	res = res + (log(pdf(d1, prms.uoE)));


	boost::math::normal_distribution<double> d2(0.035, 0.0066);//uoL
	res = res + (log(pdf(d2, prms.uoL)));


	boost::math::normal_distribution<double> d3(0.10, 0.035);//uP
	res = res + (log(pdf(d3, prms.uP)));


	boost::math::normal_distribution<double> uM1(0.091,0.005);//uM 
	res = res + (log(pdf(uM1, prms.uM)));


	boost::math::normal_distribution<double> d4(13.06, 2.5);//Y
	res = res + (log(pdf(d4, prms.Y)));

	//boost::math::uniform_distribution<double> d4x(0.1, 70);//Y
	//res = res + (log(pdf(d4x, prms.Y)));

	//boost::math::uniform_distribution<double> u5(0.000001, 1);//w unif
	//res = res + (log(pdf(u5, prms.w)));
	
	boost::math::uniform_distribution<double> u61(0.01, 9);//z1 unif
	boost::math::uniform_distribution<double> u62(0.01, 9);//z2 unif
	boost::math::uniform_distribution<double> u63(0.01, 9);//z3 unif
	boost::math::uniform_distribution<double> u64(0.01, 9);//z4 unif
	boost::math::uniform_distribution<double> u65(0.01, 9);//z5 unif
	boost::math::uniform_distribution<double> u66(0.01, 9);//z6 unif
	boost::math::uniform_distribution<double> u67(0.01, 9);//z7 unif
	boost::math::uniform_distribution<double> u68(0.01, 9);//z8 unif
	boost::math::uniform_distribution<double> u69(0.01, 9);//z7 unif
	boost::math::uniform_distribution<double> u610(0.01, 9);//z8 unif

	res = res + (log(pdf(u61, prms.z1)));
	res = res + (log(pdf(u62, prms.z2)));
	res = res + (log(pdf(u63, prms.z3)));
	res = res + (log(pdf(u64, prms.z4)));
	res = res + (log(pdf(u65, prms.z5)));
	res = res + (log(pdf(u66, prms.z6)));
	res = res + (log(pdf(u67, prms.z7)));
	res = res + (log(pdf(u68, prms.z8)));
	res = res + (log(pdf(u69, prms.z9)));
	res = res + (log(pdf(u610, prms.z10)));

	if (dFunc == "powerClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> d8(0.01, 20);//o
		if ((log(pdf(d8, prms.o))) == -inf) res = -inf; else res = res;
	}
	if (dFunc == "powerNoClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> d8(0.01, 20);//o
		if ((log(pdf(d8, prms.o))) == -inf) res = -inf; else res = res;
	}
	if (dFunc == "expClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;
	}
	if (dFunc == "expNoClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;
	}
	if (dFunc == "linearClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;

	}
	if (dFunc == "linearNoClumped") {
		boost::math::uniform_distribution<double> u71(0.01, 9);//sf1 unif
		boost::math::uniform_distribution<double> u72(0.01, 9);//sf2 unif
		boost::math::uniform_distribution<double> u73(0.01, 9);//sf3 unif
		boost::math::uniform_distribution<double> u74(0.01, 9);//sf4 unif
		boost::math::uniform_distribution<double> u75(0.01, 9);//sf5 unif
		boost::math::uniform_distribution<double> u76(0.01, 9);//sf6 unif
		boost::math::uniform_distribution<double> u77(0.01, 9);//sf7 unif
		boost::math::uniform_distribution<double> u78(0.01, 9);//sf8 unif
		boost::math::uniform_distribution<double> u79(0.01, 9);//sf9 unif
		boost::math::uniform_distribution<double> u710(0.01, 9);//sf10 unif
		if ((log(pdf(u71, prms.sf1))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u72, prms.sf2))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u73, prms.sf3))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u74, prms.sf4))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u75, prms.sf5))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u76, prms.sf6))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u77, prms.sf7))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u78, prms.sf8))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u79, prms.sf9))) == -inf) res = -inf; else res = res;
		if ((log(pdf(u710, prms.sf10))) == -inf) res = -inf; else res = res;
		boost::math::uniform_distribution<double> u8(1, 30);//n unif
		if ((log(pdf(u8, prms.n))) == -inf) res = -inf; else res = res;
	}

	//dE = 0.15, dL = 0.269, dP = 1.563

	

	boost::math::normal_distribution<double> d5(0.149925, 0.03);//dE
	res = res + (log(pdf(d5, prms.dE)));

	boost::math::normal_distribution<double> d6(0.2398082, 0.1);//dL
	res = res + (log(pdf(d6, prms.dL)));

	boost::math::uniform_distribution<double> d6Uni(0.1, 0.5);//dL
	res = res + (log(pdf(d6Uni, prms.dL)));


	boost::math::normal_distribution<double> d7(1, 0.23);//dP
	res = res + (log(pdf(d7, prms.dP)));

	boost::math::uniform_distribution<double> d7p(0.2, 1.5);//dP
	res = res + (log(pdf(d7p, prms.dP)));

	boost::math::normal_distribution<double> dtau(7, 2.5);//tau
	res = res + (log(pdf(dtau, round(prms.tau))));




	boost::math::uniform_distribution<double> mm(0.1, 10);//mg
	res = res + (log(pdf(mm, prms.Mg)));


	//boost::math::uniform_distribution<double> lk(0.0000001, 1);//mg
	//res = res + (log(pdf(lk, prms.lK)));

	return(res);
}



//std::vector<int> gen_seeds(int num) {
//	std::set<int> s;
//	std::random_device g;
//	while (s.size() < num)
//		s.insert(g());
//	std::vector<int> seeds(s.begin(), s.end());
//	return seeds;
//}



/*pMMH sampler function
@param initParams initial parameter values
@param fixedParam fixed parameter(s)
@param sdProps starting s.d. for proposals for each parameter
@param acptRs target acceptance ratios for each parameter
@param fitPrms tuple containing name and values for each parameter 
@param maxSdProps maximum s.d. for each parameter proposal to stop it going too high during tuning
@param niter number of iterations in MMH algorithm
@param nburn number of iterations to automatically burn
@param monitoring if true print output every x number of iterations according to @param tell
@param startAdapt when to start adapting proposal s.d.'s
@param tell how often to print outputs
@param oDat struct containing observed data@return struct containing results of pMMH*/
pMMHres pMMHSampler(
	modParms initParams,
	string dFunc,
	vector<double> sdProps,
	vector<double> acptRs,
	vector<tuple<string, double>> fitPrms,
	vector<double> maxSdProps,
	int niter,
	int particles,
	int nburn,
	bool monitoring,
	int startAdapt,
	int	tell,
	obsDatX	oDat){
	int iter = 0;
	unsigned int parmNum = 0;
	double llProp;//proposed ll
	double llCur;//current ll
	vector<double> llRes;
	double llRatio; //ratio between proposed and current ll
	modParms prms = initParams;
	double propPrm; //proposed new parameter
	//int seedIter;

	if (dFunc != "expNoClumped" && dFunc != "linearNoClumped" && dFunc != "powerNoClumped" && dFunc != "expClumped" && dFunc != "linearClumped" && dFunc != "powerClumped" && dFunc != "logisticClumped" && dFunc != "logisticNoClumped") {
		cerr << "dFunc must equal correct value: linearNoClumped, powerNoClumped, expNoClumped, linearClumped, powerClumped or expClumped";
		cin.get();
	}

	//std::vector<int> seedStore = gen_seeds(particles*niter);
	//vector<int> seeds(seedStore.begin(), seedStore.begin() + particles);
	if (dFunc == "linearNoClumped" || dFunc == "linearClumped" || dFunc == "expClumped" || dFunc == "expNoClumped") {
		prms.o = 1;
	}
		llCur = llFunc(particles, prms, oDat, dFunc) + lprior(prms, dFunc); //get value for initial ll
	vector<double> acptRcur = acptRs;//current acceptance ratio
	vector<double> acpts(acptRs.size(), 0.0);//number of acceptances (use acptRs to get correct vector length)
	vector<int> parmIter(acptRs.size(), 0);//iteration number for specific parameters
	pMMHres results;



	for (auto iter = 0; iter != niter; ++iter) {

		if (maxSdProps[parmNum] != 0) {

			propPrm = propPrmFunc(sdProps[parmNum], get<1>(fitPrms[parmNum]));//propose new parameter

			prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), propPrm);//update parameter struct

			if ((monitoring = true && iter % tell == 0)) {

				cout << endl << "||-----------------------||" << dFunc << "||-----------------------||" << endl;
				cout << "iteration " << iter << " of " << niter << endl;
				cout << " uoE = " << prms.uoE << " uoL = " << prms.uoL << " uP = " << prms.uP << " uM = " << prms.uM << " Y = " << prms.Y << " w = " << prms.w << " n = " << prms.n << " z1 = " << prms.z1 << endl
					<< " z2 = " << prms.z2 << " z3 = " << prms.z3 << " z4 = " << prms.z4 << " z5 = " << prms.z5 << " z6 = " << prms.z6 << " z7 = " << prms.z7 << " z8 = " << prms.z8 << " z9 = " << prms.z9 << " z10 = " << prms.z10
					<< " sf1 = " << prms.sf1 << " sf2 = " << prms.sf2 << " sf3 = " << prms.sf3 << " sf4 = " << prms.sf4 << " sf5 = " << prms.sf5 << " sf6 = " << prms.sf6 << " sf7 = " << prms.sf7 << " sf8 = " << prms.sf8 << " sf9 = " << prms.sf9 << " sf10 = " << prms.sf10
					<< "dE = " << prms.dE << " dL = " << prms.dL << " dP = " << prms.dP << " o = " << prms.o << " Mg = " << prms.Mg << " p = " << prms.p << " tau = " << prms.tau << " lK = " << prms.lK << " lKs = " << prms.lKs << " lKm = " << prms.lKm << endl;
				cout << "||---------aratio--------||" << endl;


				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << acptRcur[iter] << " ";
				}
				cout << endl << "||---------aNum--------||" << endl;
				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << acpts[iter] << " ";
				}
				cout << endl << "||---------sd--------||" << endl;
				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << sdProps[iter] << " ";
				}
				cout << endl << "||-----------------------||" << endl;
			}

			/*ofstream myfile;
			myfile.open("Q:\\test.csv");
				myfile << " Mg = " << prms.Mg << endl;*/

			//vector<int> seeds(seedStore.begin()+ seedIter, seedStore.begin()+ seedIter + particles);
			if (dFunc == "linearNoClumped" || dFunc == "linearClumped" || dFunc == "expClumped" || dFunc == "expNoClumped") {
				prms.o = 1;
			}
			llProp = llFunc(particles, prms, oDat, dFunc) + lprior(prms, dFunc);//+prms.sf+prms.z;//find log likelihood from particle filter
			//seedIter = seedIter + particles;

			//boost::math::normal_distribution<double> curr(get<1>(fitPrms[parmNum]), sdProps[parmNum]);//normal dist for metHast 
			//boost::math::normal_distribution<double> prp(propPrm, sdProps[parmNum]);//normal dist for metHast 

			//llRatio = (llProp -log(pdf(curr, propPrm))) - (llCur- log(pdf(prp, get<1>(fitPrms[parmNum]))));//ratio between previous ll and proposed ll
			double llRatio = llProp - llCur;//ratio between previous ll and proposed ll


			//print outputs
			if ((monitoring = true && iter % tell == 0)) {
				cout << "current = " << llCur << endl;
				cout << "proposed = " << llProp << endl;
				cout << "llratio = " << llRatio << endl;
				cout << "current parameter = " << get<0>(fitPrms[parmNum]) << endl;

			}


			if (llRatio >= 0 || rn01() <= exp(llRatio)) {
				llCur = llProp; //update current ll if ll is accepted
				acpts[parmNum] = acpts[parmNum] + 1;
				get<1>(fitPrms[parmNum]) = propPrm;
			}
			else {
				prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), get<1>(fitPrms[parmNum])); //change parameter back if not accepted
			}

			if (iter > startAdapt) {
				//start adapting with tuner 
				acptRcur[parmNum] = (acpts[parmNum] / parmIter[parmNum]); //calc current acceptance ratio for parameter
				sdProps[parmNum] = tuner(sdProps[parmNum], acptRs[parmNum], acptRcur[parmNum], maxSdProps[parmNum]);
			}

			//if passed burn-in, start adding to results
			if (iter > nburn) {
				results.uoE.emplace_back(prms.uoE);
				results.uoL.emplace_back(prms.uoL);
				results.uP.emplace_back(prms.uP);
				results.uM.emplace_back(prms.uM);

				results.Y.emplace_back(prms.Y);
				results.w.emplace_back(prms.w);
				results.n.emplace_back(prms.n);
				results.z1.emplace_back(prms.z1);
				results.z2.emplace_back(prms.z2);
				results.z3.emplace_back(prms.z3);
				results.z4.emplace_back(prms.z4);
				results.z5.emplace_back(prms.z5);
				results.z6.emplace_back(prms.z6);
				results.z7.emplace_back(prms.z7);
				results.z8.emplace_back(prms.z8);
				results.z9.emplace_back(prms.z9);
				results.z10.emplace_back(prms.z10);

				results.sf1.emplace_back(prms.sf1);
				results.sf2.emplace_back(prms.sf2);
				results.sf3.emplace_back(prms.sf3);
				results.sf4.emplace_back(prms.sf4);
				results.sf5.emplace_back(prms.sf5);
				results.sf6.emplace_back(prms.sf6);
				results.sf7.emplace_back(prms.sf7);
				results.sf8.emplace_back(prms.sf8);
				results.sf9.emplace_back(prms.sf9);
				results.sf10.emplace_back(prms.sf10);

				results.dE.emplace_back(prms.dE);
				results.dL.emplace_back(prms.dL);
				results.dP.emplace_back(prms.dP);
				results.o.emplace_back(prms.o);
				results.Mg.emplace_back(prms.Mg);
				results.p.emplace_back(prms.p);
				results.tau.emplace_back(prms.tau);
				results.lK.emplace_back(prms.lK);
				results.lKs.emplace_back(prms.lKs);
				results.lKm.emplace_back(prms.lKm);


				results.ll.emplace_back(llCur);

				if ((monitoring = true && iter % (tell+5000) == 0)) {

					try {
						string medFile = "C:\\ImperialOld\\Imperial\\lModCpp\\Data\\meds";
						medFile.append(dFunc);
						medFile.append("resUpdates.txt");
						double o_min_value = *std::min_element(results.o.begin(), results.o.end());

						ofstream myfile;
						myfile.open(medFile);
						myfile << "current dP median = " << medianFnc(results.dP) << endl << "current n median = " << medianFnc(results.n) << endl << "current uoE median = " << medianFnc(results.uoE) << endl << "current Y median = " << medianFnc(results.Y) << endl
							<< "current dL median = " << medianFnc(results.dL) << endl << "current dE median = " << medianFnc(results.dE) << endl << "current o minimum = " << o_min_value << endl
							<< "current o median = " << medianFnc(results.o) << endl << "current uM median = " << medianFnc(results.uM) << endl << "current tau median = " << medianFnc(results.tau) << endl
							<< "current ll = " << medianFnc(results.ll) << endl << " iteration = " << iter;
					}
					catch (...) {
						cerr << "something wrong with writing median vals...better luck next time?";
					}
					
					}
			}

			parmNum++;

			; //parameter number, if greater than number of parms, revert back to 0
			if (parmNum >= boost::size(sdProps)) {
				parmNum = 0;
			}

			parmIter[parmNum]++;
		}
		else {
			parmNum++;
			if (parmNum >= boost::size(sdProps)) {
				parmNum = 0;
			}
			parmIter[parmNum]++;
		}

	}

	return results;
}