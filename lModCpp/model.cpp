#include "lModH.h"


vector<tuple<double, double, double, double,double>> mPmod(modParms parmsx, boost::mt19937_64 rd, string dFunc) {

	int t = parmsx.startTime;
	int time = parmsx.endTime;
	double Be = 0;
	double Bl = 0;
	double Bp = 0;
	double Bm = 0;
	double nt = 0;
	double E = parmsx.E0;
	double L = parmsx.L0;
	double P = parmsx.P0;
	double M = parmsx.M0;
	double K;
	double trx = round(parmsx.tau / parmsx.dt);
	double dt = parmsx.dt;
	double  uoE = parmsx.uoE;
	double  uoL = parmsx.uoL;
	double  Y = parmsx.Y;
	double  n = parmsx.n;
	double sf = parmsx.sf;
	double o = parmsx.o;
	double uP = parmsx.uP;
	double S = parmsx.S;
	double dE = parmsx.dE;// parmsx.dE;
	double dL = parmsx.dL;
	double uM = parmsx.uM;
	double uE;
	double uL;
	double uN;
	double rEff;
	vector<double> sortrF;
	double mRan;
	double dP = parmsx.dP;
	double Mg = parmsx.Mg;
	double uEn;
	double uLn;
	double lK = pow(10,parmsx.lK);
	double lKs = pow(10,parmsx.lKs);
	double lKm = parmsx.lKm;
	double rFsum = 0.0;

	vector<double> rF = parmsx.rF;
	vector<tuple<double, double, double, double,double>> r;
	vector<double> rFsum2;

	for (int r = 0; r <= time; r++) {

		rFsum2.emplace_back(exp((-(time - r)) / trx));

	}

	while (t < time) {

	//	for (int r = 0; r <= t; r++) {
		//	rFsum = rFsum + (rFsum2[r] * rF[r]);
		//}

		rFsum = 0;

		for (int r = 0; r <= t; r++) {
			rFsum = rFsum + rFsum2[time - t + r] * rF[r];
			}

		K = sf * (1.0 / (trx * (1 - exp(-(t + 1) / trx)))) * rFsum;

		

	
	
		if (dFunc == "powerNoClumped" || dFunc == "linearNoClumped" || dFunc == "powerClumped" || dFunc == "linearClumped" || dFunc == "logisticClumped"|| dFunc == "logisticNoClumped") {
			uE = uoE * (1 + pow(((E + L) / (K)), o));
			uL = uoL * (1 + (Y*pow(((E + L) / (K)), o)));

		}
		else {

			uE = uoE * exp(((E + L) / (K)));
			uL = uoL * exp(Y*(((E + L) / (K))));

		}


		if (uL < 0)
			uL = 0;
		if (uE < 0)
			uE = 0;

			if ((dE + uE)*dt < 1) {
				boost::binomial_distribution<int> distributionBe(round(E), (dE + uE)*dt);
				Be = distributionBe(rd);
			}
			else {
				boost::binomial_distribution<int> distributionBe(round(E), 1);
				Be = distributionBe(rd);
			}
		
		
	
			if ((dL + uL)*dt < 1) {
				boost::binomial_distribution<int> distributionBl(round(L), (dL + uL)*dt);
				Bl = distributionBl(rd);
			}
			else {
				boost::binomial_distribution<int> distributionBl(round(L),1);
				Bl = distributionBl(rd);
			}


			if ((dP + uP)*dt < 1) {
				boost::binomial_distribution<int> distributionBp(round(P), (dP + uP)*dt);
				Bp = distributionBp(rd);
			}
			else {
				boost::binomial_distribution<int> distributionBp(round(P), 1);
				Bp = distributionBp(rd);
			}
		
	
			if (M >= 1) {
				boost::binomial_distribution<int> distributionBm(round(M), uM*dt);
				Bm = distributionBm(rd);
			}
			else Bm = 0;


			if (dFunc == "expClumped"|| dFunc == "linearClumped"|| dFunc == "powerClumped" || dFunc == "logisticClumped") {
				if (M >= 1) {
					boost::binomial_distribution<int> distributionNt(round(M), (dt / S));
					nt = distributionNt(rd);
				
				}
				else nt = 0;

				if (n*nt > 0) {
					boost::poisson_distribution<long unsigned int> distributionRp(n*nt);
					E = (E - Be + distributionRp(rd));
				}
				else E = (E - Be);
			}

			else {
				if (M * (n*dt) > 0) {
				//	boost::poisson_distribution<long unsigned int> distributionRp2(M * (n*dt));
					E = (E - Be + (M * (n*dt)));
				}
				else E = (E - Be);
				
			}


		if (Be >= 1) {
			boost::binomial_distribution<int> distributionL(Be, (dE / (uE + dE)));
			L = (L - Bl + distributionL(rd));
		}
		else 
			L = (L - Bl);

		if (Bl >= 1) {
			boost::binomial_distribution<int> distributionP(Bl, (dL / (uL + dL)));
			P = (P - Bp + distributionP(rd));
		}
		else 
			P = (P - Bp);

		if (Bp >= 1) {
			boost::binomial_distribution<int> distributionM(Bp, (dP / (uP + dP)));
			mRan = distributionM(rd);
		}
		else
			mRan = 0; 

		if (Mg > 0) {
			boost::poisson_distribution<long unsigned int> distributionRpMG(Mg);
			M = (M + (0.5*mRan) - Bm + distributionRpMG(rd));
		}
		
		else 
			M = (M + (0.5*mRan) - Bm);

	//	if (M < 1)
		//	M = 1;

		//round M up or down 50% of the time
		if (round(M) - M == 0.5){
			boost::binomial_distribution<int> distributionMtss(1, 0.5);
			int Mtss = distributionMtss(rd);
			if (Mtss == 1) M = M - 0.5; else M = round(M);
		}

		if (dFunc == "expClumped" || dFunc == "linearClumped" || dFunc == "powerClumped" || dFunc == "logisticClumped") {
			rEff = 0.5*((n) / (exp(uM*S) - 1))*(1 / (1 + uE / dE))*(1 / (1 + uL / dL))*(1 / (1 + (uP) / dP));

		} else { 
			double Es = n * (exp(S*uM) - 1) / uM;//calculate value for oviposition cycle from daily rate
			rEff = 0.5*((Es) / (exp(uM*S) - 1))*(1 / (1 + uE / dE))*(1 / (1 + uL / dL))*(1 / (1 + (uP) / dP));}

		t++;
		r.emplace_back(make_tuple(E, L, P, M, rEff));

	}


	return r;
}

