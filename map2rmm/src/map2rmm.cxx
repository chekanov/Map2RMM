/**************************************************
*       Project event into RMM	
**************************************************/

using namespace std;
#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "SystemOfUnits.h"
#include "LParticle.h"
#include "CParticle.h"
#include <string>
#include <vector>

/*

 We build matrix like this assuming maxN=2:
 j -jet
 m- muons
 e- electron
 g -photon

 All variables are ratio to CMS 

 
     |   MET    j1       j2       mu1     mu2     e1      e2     g1      g2
 ------------------------------------------------------------------------- 
 MET |  MT^2   MT(j1)  MT(j2)   MT(m1)   MT(m2)  MT(e1)  MT(e2) MT(g1) MT(g2)  
     |      
 j1  |  h_j1     et(j1) m(j1,j2) 
     |  
 j2  |  h_j2               et(j1)
     |
 m1  |                         et(m1) 
     | 
 m2  |                                 et(m2) 
     |
 e1  |   
     |
 e2  |   
     |
 g1  |   
     |
 g2  |  


*/



// return angle product
float getAngle(const float CMS, const TLorentzVector p1, const TLorentzVector p2){
	//float ang=p1.Angle(p2.Vect());
	//ang=(1-TMath::Cos(ang));
	//double  ene=sqrt(p1.E()*p2.E())/CMS;
	double y1=p1.Rapidity();
	double y2=p2.Rapidity();
	double HL=TMath::CosH( 0.5*(y2-y1) )-1;
	return float(HL);
}

// get masses
float getMass(const float CMS, const TLorentzVector p1, const TLorentzVector p2){
	TLorentzVector pp=p1+p2;
	float xmass=pp.M()/CMS;
	//if (xmass>1) cout << "xmass=" << xmass << " CMS=" << CMS << endl;
	return xmass;
}

// get cosh(y)-1
// y=0 correspongs 0
// y=4 correspondsi 2.7/100
float getHL(const TLorentzVector p1){
	double y=p1.Rapidity();
	double HL=TMath::CosH(y)-1;
	//cout << "y=" << y << " HL=" << HL << endl;
	return (float)(HL);
}


// create matrix with event projection
// maxN max number of particles of a given type
// maxNumberTypes total number of types
float**  map2rmm(const float CMS, const int maxN, const int maxNumberTypes,
                      const vector<LParticle> missing,
                      const vector<LParticle> jets,
                      const vector<LParticle> muons,
                      const vector<LParticle> electrons,
                      const vector<LParticle> photons) {


	//just to make sure, sort
	/*
	if (jets.size()>1)  std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
	if (muons.size()>1) std::sort(muons.begin(), muons.end(), greater<LParticle>() ) ;
	if (electrons.size()>1) std::sort(electrons.begin(), electrons.end(), greater<LParticle>() ) ;
	if (photons.size()>1) std::sort(photons.begin(), photons.end(), greater<LParticle>() ) ;
	*/

	const int maxNumber=maxN; // max number for each object
	const int maxTypes=maxNumberTypes;  // max numbers of types
	const int maxSize=maxNumber*maxTypes+1; // + MET

	const int height=maxSize;
	const int width=maxSize;

	float** outMatrix = 0;
	outMatrix = new float*[height];
	for (int h=0; h<height; h++) outMatrix[h] = new float[width];

	unsigned  int INCR=0;

	for(int i=0; i<maxSize; i++)
		for(int j=0; j<maxSize; j++)  outMatrix[i][j]=0;


	// invariant masses of same objects
	INCR=0;
	int ii=0;
	LParticle MET=(LParticle)missing.at(ii);
	TLorentzVector LMET=MET.GetP();
	outMatrix[ii+INCR*maxNumber][ii+INCR*maxNumber]=LMET.Et()/CMS;


	INCR=0;
	unsigned  int mJets=maxNumber;
	if (jets.size()<mJets) mJets=jets.size();
	if (mJets>0) {
		for (unsigned int k1 = 0; k1<mJets; k1++) {
			LParticle LPP1=(LParticle)jets.at(k1);
			TLorentzVector LP1=LPP1.GetP();
			TLorentzVector LM=LP1+LMET;
			if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=LM.Mt()/CMS;
			outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

			if (k1 == 0) {
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
			} else if (k1>0)  {
				LParticle LPP3=(LParticle)jets.at(k1-1);
				TLorentzVector LP3=LPP3.GetP();
				float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
			}

			//  non-diagonal
			for (unsigned int k2 = 0; k2<mJets; k2++) {
				LParticle LPP2=(LParticle)jets.at(k2);
				TLorentzVector LP2=LPP2.GetP();
				if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
				if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
			}
		}
	}


	INCR=1;
	unsigned  int mMuons=maxNumber;
	if (muons.size()<mMuons) mMuons=muons.size();
	if (mMuons>0) {
		for (unsigned int k1 = 0; k1<mMuons; k1++) {
			LParticle LPP1=(LParticle)muons.at(k1);
			TLorentzVector LP1=LPP1.GetP();
			TLorentzVector LM=LP1+LMET;
			if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=LM.Mt()/CMS;
			outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

			if (k1 == 0) {
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
			} else if (k1>0)  {
				LParticle LPP3=(LParticle)muons.at(k1-1);
				TLorentzVector LP3=LPP3.GetP();
				float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
			}


			for (unsigned int k2 = 0; k2<mMuons; k2++) {
				LParticle LPP2=(LParticle)muons.at(k2);
				TLorentzVector LP2=LPP2.GetP();
				if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
				if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
			}}
	}


	INCR=2;
	unsigned  int mEle=maxNumber;
	if (electrons.size()<mEle) mEle=electrons.size();
	if (mEle>0) {
		for (unsigned int k1 = 0; k1<mEle; k1++) {
			LParticle LPP1=(LParticle)electrons.at(k1);
			TLorentzVector LP1=LPP1.GetP();
			TLorentzVector LM=LP1+LMET;
			if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=LM.Mt()/CMS;
			outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

			if (k1 == 0) {
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
			} else if (k1>0)  {
				LParticle LPP3=(LParticle)electrons.at(k1-1);
				TLorentzVector LP3=LPP3.GetP();
				float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
			}

			// non-diagonal
			for (unsigned int k2 = 0; k2<mEle; k2++) {
				LParticle LPP2=(LParticle)electrons.at(k2);
				TLorentzVector LP2=LPP2.GetP();
				if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
				if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
			}}
	}


	INCR=3;
	unsigned  int mPho=maxNumber;
	if (photons.size()<mPho) mPho=photons.size();
	if (mPho>0) {
		for (unsigned int k1 = 0; k1<mPho; k1++) {
			LParticle LPP1=(LParticle)photons.at(k1);
			TLorentzVector LP1=LPP1.GetP();
			TLorentzVector LM=LP1+LMET;
			if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=LM.Mt()/CMS;
			outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

			if (k1 == 0) {
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
			} else if (k1>0)  {
				LParticle LPP3=(LParticle)photons.at(k1-1);
				TLorentzVector LP3=LPP3.GetP();
				float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
				outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
			}


			// non-diagonal
			for (unsigned int k2 = 0; k2<mPho; k2++) {
				LParticle LPP2=(LParticle)photons.at(k2);
				TLorentzVector LP2=LPP2.GetP();
				if (k1<k2) outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
				if (k1>k2) outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
			}}
	}




	// diagonal ellements, such as jet-lep, etc. First, do jets:
	unsigned  int INCR_SHIFT=0; // for jets =1
	for (unsigned int k1 = 0; k1<mJets; k1++) {
		LParticle LPP1=(LParticle)jets.at(k1);
		TLorentzVector LP1=LPP1.GetP();

		INCR=1;
		for (unsigned int k2 = 0; k2<mMuons; k2++) {
			LParticle LPP2=(LParticle)muons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

		INCR=2;
		for (unsigned int k2 = 0; k2<mEle; k2++) {
			LParticle LPP2=(LParticle)electrons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

		INCR=3;
		for (unsigned int k2 = 0; k2<mPho; k2++) {
			LParticle LPP2=(LParticle)photons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

	}



	// muons
	INCR_SHIFT=1; // for muons =2
	for (unsigned int k1 = 0; k1<mMuons; k1++) {
		LParticle LPP1=(LParticle)muons.at(k1);
		TLorentzVector LP1=LPP1.GetP();

		INCR=2;
		for (unsigned int k2 = 0; k2<mEle; k2++) {
			LParticle LPP2=(LParticle)electrons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

		INCR=3;
		for (unsigned int k2 = 0; k2<mPho; k2++) {
			LParticle LPP2=(LParticle)photons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

	}


	// electrons
	INCR_SHIFT=2; // for electrons =3
	for (unsigned int k1 = 0; k1<mEle; k1++) {
		LParticle LPP1=(LParticle)electrons.at(k1);
		TLorentzVector LP1=LPP1.GetP();
		INCR=3;
		for (unsigned int k2 = 0; k2<mPho; k2++) {
			LParticle LPP2=(LParticle)photons.at(k2);
			TLorentzVector LP2=LPP2.GetP();
			TLorentzVector LPP=LP1+LP2;
			outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
			outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
		}

	}




	return outMatrix;


}
