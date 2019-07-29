#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include <TLorentzVector.h>
#include "AddXsecBranch/MT2_ROOT.h"
#include "AddXsecBranch/MyAMT2.h"

using namespace std;

void Execute( string treeN, string inpath, string fileN, string sample, string outpath, string varnames, string camp ) {

	inpath += "/";
	outpath += "/";

	cout << "Processing: " << fileN << endl;

	// opening file/tree within a RDataFrame object
	ROOT::RDataFrame df( treeN.c_str(), (inpath+fileN).c_str() );
	cout << "RDataFrame successfully loaded!!!!" << endl;

	// defining lamda functions of the new variables to add

	// DeltaEta
	//auto cDeltaEta = []( vector<float> LepEta ) { return ( LepEta.size()>1 ? TMath::Abs(LepEta[0]-LepEta[1]) : -1. ); };

	// claculate mT
	auto calculateMT = []( TLorentzVector tlv, float eT_miss, float EtMissPhi ) {
		TVector2 v2;
		v2.SetMagPhi( eT_miss, EtMissPhi );
		TLorentzVector met_vec;
		met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
		float mm = (tlv.Mt() + met_vec.Mt())*(tlv.Mt() + met_vec.Mt()) - (tlv+met_vec).Perp2();
		return ( mm>=0. ? TMath::Sqrt(mm) : TMath::Sqrt(-mm) );
	};

	// mT
	auto cmT = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>0 ) {
			TLorentzVector lep_vec;
                        lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TVector2 v2; 
			v2.SetMagPhi( eT_miss, EtMissPhi );
			TLorentzVector met_vec;
			met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
			float mm = (lep_vec.Mt() + met_vec.Mt())*(lep_vec.Mt() + met_vec.Mt()) - (lep_vec+met_vec).Perp2();
			return ( mm>=0. ? TMath::Sqrt(mm) : TMath::Sqrt(-mm) );
		}
		//else 
			return (double)(-1);
	};

	// meff
	auto cmeff = []( vector<float> LepPt, vector<float> JetPt, float eT_miss ) {
		if ( LepPt.size()>0 ) {
			float mm = eT_miss + LepPt[0] + LepPt[1];
			for ( int i=0 ; i<JetPt.size() ; i++ ) mm += JetPt[i];
			return mm;
		}
		//else
			return  (float)(-1);
	};

	// mt2
	auto cmt2 = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>0 ) {
			TLorentzVector Lep1, Lep2;
			Lep1.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			Lep2.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			MyAMT2 * mt2tool = new MyAMT2();
			TLorentzVector vds = TLorentzVector(0, 0, 0, 0);
			TVector2 ptm;
			ptm.SetMagPhi( eT_miss, EtMissPhi );
			float m_t2 = mt2tool->MT2calc( Lep1, Lep2, vds, ptm );
			return m_t2;
		}
			return (float)(-1);
	};

	// MET/meff
	auto cMETmeff = []( vector<float> LepPt, vector<float> JetPt, float eT_miss ) {
		if ( LepPt.size()>0 ) {
			float mm = eT_miss + LepPt[0] + LepPt[1];
			for ( int i=0 ; i<JetPt.size() ; i++ ) mm += JetPt[i];
			return ( eT_miss/mm );
		}
		//else
			return (float)(-1);
	};

	// lepton closest to jet(s)
	auto closestLep = []( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		if ( LepPt.size()>0 ) {
			if ( JetPt.size() == 1 ) {	
				TLorentzVector jet_vec; 
				jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
				TLorentzVector lep_vec;
				lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
				float dr = jet_vec.DeltaR( lep_vec ); // distance from leading lepton
				TLorentzVector lep_vec1;
				lep_vec1.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] ); 
				float dr1 = jet_vec.DeltaR( lep_vec1 ); // distance from sub-leading lepton
				if ( dr <= dr1 )
					return 0;
				else
					return 1;
			}
			if ( JetPt.size() >= 2 ) {	
				TLorentzVector jet_vec; 
				jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
				TLorentzVector jet_vec1; 
				jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
				TLorentzVector dijet_vec = jet_vec + jet_vec1;
				TLorentzVector lep_vec;
				lep_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] ); 
				float dr = dijet_vec.DeltaR( lep_vec ); // distance from leading lepton
				TLorentzVector lep_vec1;
				lep_vec1.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] ); 
				float dr1 = dijet_vec.DeltaR( lep_vec1 ); // distance from sub-leading lepton
				if ( dr <= dr1 )
					return 0;
				else
					return 1;
			}
		}
		return (-1);
	};

	// mlj
	auto cmlj = [=]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		int icl = closestLep( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
		if ( icl>-1 ) { 
			TLorentzVector lep_vec;
			lep_vec.SetPtEtaPhiM( LepPt[icl], LepEta[icl], LepPhi[icl], LepM[icl] ); 
			TLorentzVector jet_vec; 
			if ( JetPt.size() == 1 )	
				jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
			else if ( JetPt.size() >= 2 ) {	
				TLorentzVector jet_vec0; 
				jet_vec0.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
				TLorentzVector jet_vec1; 
				jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
				jet_vec = jet_vec0 + jet_vec1;
			}
			return ( ( jet_vec + lep_vec ).M() );
		}
		return (double)(-1);
	};

	auto cLumi = [=]( float TotalWeight ) {
		if ( TotalWeight == 0. )	return 0.;
		if ( camp == "a" )		return 36100.;
		else if ( camp == "d" )		return 44100.;
		else if ( camp == "e" )		return 60000.;
		else 													return 0.;
	};

	// mTW
	auto cmTW = [=]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		int icl = closestLep( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
		if ( icl==-1) return (float)(-1);
		int ifl = (icl==0) ? 1 : 0;

		TLorentzVector lep_vec;
		lep_vec.SetPtEtaPhiM( LepPt[ifl], LepEta[ifl], LepPhi[ifl], LepM[ifl] ); 
	
		float mt = calculateMT( lep_vec, eT_miss, EtMissPhi );
		return mt;
	};

	// mTH
	auto cmTH = [=]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		int icl = closestLep( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
		if ( icl==-1) return (float)(-1);

		TLorentzVector lep_vec;
		lep_vec.SetPtEtaPhiM( LepPt[icl], LepEta[icl], LepPhi[icl], LepM[icl] ); 
	
		TLorentzVector jet_vec; 
		if ( JetPt.size() == 1 )	
			jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
		else if ( JetPt.size() >= 2 ) {	
			TLorentzVector jet_vec0; 
			jet_vec0.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
			TLorentzVector jet_vec1; 
			jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
			jet_vec = jet_vec0 + jet_vec1;
		}

		float mt = calculateMT( lep_vec+jet_vec, eT_miss, EtMissPhi );
		return mt;
	};

	// mt2
	auto cmt2tag = [=]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {

		int icl = closestLep( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
		if ( icl==-1) return (float)(-1);
		int ifl = (icl==0) ? 1 : 0;
		TLorentzVector lepc_vec, lepf_vec;
		lepc_vec.SetPtEtaPhiM( LepPt[icl], LepEta[icl], LepPhi[icl], LepM[icl] );
		lepf_vec.SetPtEtaPhiM( LepPt[ifl], LepEta[ifl], LepPhi[ifl], LepM[ifl] );

		TLorentzVector jet_vec; 
		if ( JetPt.size() == 1 )	
			jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
		else if ( JetPt.size() >= 2 ) {	
			TLorentzVector jet_vec0; 
			jet_vec0.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
			TLorentzVector jet_vec1; 
			jet_vec1.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
			jet_vec = jet_vec0 + jet_vec1;
		}

		MyAMT2 * mt2tool = new MyAMT2();
		TLorentzVector vds = TLorentzVector(0, 0, 0, 0);
		TVector2 ptm;
		ptm.SetMagPhi( eT_miss, EtMissPhi );
		float m_t2 = mt2tool->MT2calc( lepc_vec+jet_vec, lepc_vec, vds, ptm );
		return m_t2;
	};

	// mTl2
	auto cmTl2 = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>0 ) {
			TLorentzVector lep_vec;
                        lep_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			TVector2 v2; 
			v2.SetMagPhi( eT_miss, EtMissPhi );
			TLorentzVector met_vec;
			met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
			float mm = (lep_vec.Mt() + met_vec.Mt())*(lep_vec.Mt() + met_vec.Mt()) - (lep_vec+met_vec).Perp2();
			return ( mm>=0. ? TMath::Sqrt(mm) : TMath::Sqrt(-mm) );
		}
		//else 
			return (double)(-1);
	};

	// mTlmin
	auto cmTlmin = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>0 ) {
			TLorentzVector lep1_vec;
                        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec;
                        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			TVector2 v2; 
			v2.SetMagPhi( eT_miss, EtMissPhi );
			TLorentzVector met_vec;
			met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
			float mm1 = (lep1_vec.Mt() + met_vec.Mt())*(lep1_vec.Mt() + met_vec.Mt()) - (lep1_vec+met_vec).Perp2();
			float mt1 = mm1>=0. ? TMath::Sqrt(mm1) : TMath::Sqrt(-mm1);
			float mm2 = (lep2_vec.Mt() + met_vec.Mt())*(lep2_vec.Mt() + met_vec.Mt()) - (lep2_vec+met_vec).Perp2();
			float mt2 = mm2>=0. ? TMath::Sqrt(mm2) : TMath::Sqrt(-mm2);
			if ( mt1 <= mt2 ) return mt1;
			else							return mt2;
		}
		//else 
			return (float)(-1);
	};

	// mTlmax
	auto cmTlmax = [=]( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>0 ) {
			TLorentzVector lep1_vec;
                        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec;
                        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			TVector2 v2; 
			v2.SetMagPhi( eT_miss, EtMissPhi );
			TLorentzVector met_vec;
			met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
			float mm1 = (lep1_vec.Mt() + met_vec.Mt())*(lep1_vec.Mt() + met_vec.Mt()) - (lep1_vec+met_vec).Perp2();
			float mt1 = mm1>=0. ? TMath::Sqrt(mm1) : TMath::Sqrt(-mm1);
			float mm2 = (lep2_vec.Mt() + met_vec.Mt())*(lep2_vec.Mt() + met_vec.Mt()) - (lep2_vec+met_vec).Perp2();
			float mt2 = mm2>=0. ? TMath::Sqrt(mm2) : TMath::Sqrt(-mm2);
			if ( mt1 >= mt2 ) return mt1;
			else							return mt2;
		}
		//else 
			return (float)(-1);
	};

	// METrel
	auto cMETrel = []( vector<float> LepPhi, vector<float> JetPhi, float eT_miss, float EtMissPhi ) {
		TVector2 v1; v1.SetMagPhi( eT_miss, EtMissPhi );
		TVector2 v2; v2.SetMagPhi( 1., LepPhi[0] );
		float dphi = fabs( v1.DeltaPhi(v2) );
		for ( size_t i=0; i<LepPhi.size(); i++ ) {
			TVector2 vt; vt.SetMagPhi( 1., LepPhi[i] );
			float dphi_t = fabs( v1.DeltaPhi(vt) );
			if ( dphi_t<dphi ) dphi = dphi_t;
		}
		for ( size_t i=0; i<JetPhi.size(); i++ ) {
			TVector2 vt; vt.SetMagPhi( 1., JetPhi[i] );
			float dphi_t = fabs( v1.DeltaPhi(vt) );
			if ( dphi_t<dphi ) dphi = dphi_t;
		}
		float METrel = eT_miss;
		if ( dphi < TMath::Pi()/2 ) METrel = eT_miss*TMath::Sin(dphi);
		return METrel;
	};

	// mll
	auto cmll = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM ) {
    if ( LepPt.size()>0 ) {
			TLorentzVector lep1_vec;
                        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec;
                        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			return ( (lep1_vec+lep2_vec).M() );
		}
		return (double)(-1);
	};

	// |mll-mZ|
	auto cdistmZ = []( vector<float> LepPt, vector<float> LepEta, vector<float> LepPhi, vector<float> LepM ) {
    if ( LepPt.size()>0 ) {
			TLorentzVector lep1_vec;
                        lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec;
                        lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			double mll = (lep1_vec+lep2_vec).M();
			return ( fabs(mll-91.1876) );
		}
		return (double)(-1);
	};

	// HT
	auto cHT = []( vector<float> JetPt ) {
		float HT = 0.;
		for ( size_t i=0 ; i< JetPt.size() ; i++ )
			HT += JetPt[i];
		return HT;
	};

	// HT1
	auto cHT1 = []( vector<float> JetPt ) {
		if ( JetPt.size()<1 ) return (float)(-1.);
		float HT = JetPt[0];
		if ( JetPt.size()>=2 ) HT += JetPt[1];
		if ( JetPt.size()>=3 ) HT += JetPt[2];
		if ( JetPt.size()>=4 ) HT += JetPt[3];
		return HT;
	};

	// HT2
	auto cHT2 = []( vector<float> JetPt ) {
		if ( JetPt.size()<2 ) return (float)(-1.);
		float HT = JetPt[1];
		if ( JetPt.size()>=3 ) HT += JetPt[2];
		if ( JetPt.size()>=4 ) HT += JetPt[3];
		return HT;
	};

	// HT3
	auto cHT3 = []( vector<float> JetPt ) {
		if ( JetPt.size()<3 ) return (float)(-1.);
		float HT = JetPt[2];
		if ( JetPt.size()>=4 ) HT += JetPt[3];
		return HT;
	};

	// LT
	auto cLT = []( vector<float> LepPt ) {
		if ( LepPt.size()<1 ) return (float)(-1.);
		float LT = LepPt[0] + LepPt[1];
		return LT;
	};

	// mjj 
	auto cmjj = []( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		if ( JetPt.size()<1 ) return (float)(-1.);
		TLorentzVector jet1_vec; jet1_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
		float mjj = JetM[0];
		if ( JetPt.size()>=2 ) {
			TLorentzVector jet2_vec; jet2_vec.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
			mjj = (float)(jet1_vec+jet2_vec).M();
		}
		return mjj;
	};

	// angular variables

	int wDist = 0; // 0->DR, 1->DPhi, 2->DEta
	int wLep = 0; // 0->lep1, 1->lep2, 2->dilep, 3->lep-closest-to-jet
	int wJet = 0; // 0->jet1, 1->jet2, 2->jet3, 3->dijet-12, 4->dijet-13, 5->dijet-23 

	auto cDistll = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM ) {
		if ( LepPt.size()>1 ) {
			TLorentzVector lep1_vec; lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec; lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			float dr = lep1_vec.DeltaR( lep2_vec );
			float dphi = fabs(lep1_vec.DeltaPhi( lep2_vec ));
			float deta = fabs(LepEta[0]-LepEta[1]);
			if			( wDist==0 ) return dr;
			else if	( wDist==1 ) return dphi;
			else if	( wDist==2 ) return deta;
			else return (float)(-1);
		}
		return (float)(-1);
	};
	auto cdRll = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM ) {
		wDist = 0; return cDistll( LepPt, LepPhi, LepEta, LepM );
        };
	auto cdPhill = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM ) {
		wDist = 1; return cDistll( LepPt, LepPhi, LepEta, LepM );
        };
	auto cdEtall = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM ) {
		wDist = 2; return cDistll( LepPt, LepPhi, LepEta, LepM );
        };

	auto cDistjMET = [&]( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		/// jet setup
		TLorentzVector jet_vec;
		if ( JetPt.size()==1 ) 
			jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
		else if ( JetPt.size()>=2 ) {
			TLorentzVector jet1_vec; jet1_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
			TLorentzVector jet2_vec; jet2_vec.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
			jet_vec = jet1_vec + jet2_vec;
		} /// finished jet setup
		else 			return (float)(-1);
		TVector2 v2; v2.SetMagPhi( eT_miss, EtMissPhi );
		TLorentzVector met_vec; met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
		float dr = jet_vec.DeltaR( met_vec );
		float dphi = fabs(jet_vec.DeltaPhi( met_vec ));
		float deta = fabs( jet_vec.Eta() - met_vec.Eta() );
		if			( wDist==0 ) return dr;
		else if	( wDist==1 ) return dphi;
		else if	( wDist==2 ) return deta;
		return (float)(-1);
	};
	// DeltaAngjetMET
	auto cdRjMET = [&]( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		wDist=0; return cDistjMET( JetPt, JetPhi, JetEta, JetM, eT_miss, EtMissPhi );
        };
	auto cdPhijMET = [&]( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		wDist=1; return cDistjMET( JetPt, JetPhi, JetEta, JetM, eT_miss, EtMissPhi );
        };
	auto cdEtajMET = [&]( vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM, float eT_miss, float EtMissPhi ) {
		wDist=2; return cDistjMET( JetPt, JetPhi, JetEta, JetM, eT_miss, EtMissPhi );
        };

	auto cDistlMET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		if ( LepPt.size()>1 ) {
			TLorentzVector lep1_vec; lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec; lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			TLorentzVector lep_vec;
			/// lepton setup
			if 			( wLep==0 ) lep_vec = lep1_vec;
			else if ( wLep==1 ) lep_vec = lep2_vec;
			else if ( wLep==2 ) lep_vec = lep1_vec + lep2_vec;
			else								return (float)(-1); // finished lepton setup
			TVector2 v2; v2.SetMagPhi( eT_miss, EtMissPhi );
			TLorentzVector met_vec; met_vec.SetPxPyPzE( v2.Px(), v2.Py(), 0., eT_miss );
			float dr = lep_vec.DeltaR( met_vec );
			float dphi = fabs(lep_vec.DeltaPhi( met_vec ));
			float deta = fabs( lep_vec.Eta() - met_vec.Eta() );
			if			( wDist==0 ) return dr;
			else if	( wDist==1 ) return dphi;
			else if	( wDist==2 ) return deta;
			else return (float)(-1);
		}
		return (float)(-1);
	};
	auto cdRl1MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=0; wLep=0; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdPhil1MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=1; wLep=0; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdEtal1MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=2; wLep=0; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdRl2MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=0; wLep=1; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdPhil2MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=1; wLep=1; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdEtal2MET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=2; wLep=1; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdRdlMET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=0; wLep=2; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdPhidlMET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=1; wLep=2; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };
	auto cdEtadlMET = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, float eT_miss, float EtMissPhi ) {
		wDist=2; wLep=2; return cDistlMET( LepPt, LepPhi, LepEta, LepM, eT_miss, EtMissPhi );
        };

	auto cDistlj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		if ( LepPt.size()>1 ) {
			TLorentzVector lep1_vec; lep1_vec.SetPtEtaPhiM( LepPt[0], LepEta[0], LepPhi[0], LepM[0] );
			TLorentzVector lep2_vec; lep2_vec.SetPtEtaPhiM( LepPt[1], LepEta[1], LepPhi[1], LepM[1] );
			TLorentzVector dilep_vec = lep1_vec + lep2_vec;
			TLorentzVector lep_vec;
			/// jet setup
			TLorentzVector jet_vec;
			if ( JetPt.size()==1 ) 
				jet_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
			else if ( JetPt.size()>=2 ) {
				TLorentzVector jet1_vec; jet1_vec.SetPtEtaPhiM( JetPt[0], JetEta[0], JetPhi[0], JetM[0] );
				TLorentzVector jet2_vec; jet2_vec.SetPtEtaPhiM( JetPt[1], JetEta[1], JetPhi[1], JetM[1] );
				jet_vec = jet1_vec + jet2_vec;
			} /// finished jet setup
			else 			return (float)(-1);
			/// closest lepton
			float drl1 = lep1_vec.DeltaR( jet_vec );
			float drl2 = lep2_vec.DeltaR( jet_vec );
			TLorentzVector lepc_vec = drl2>=drl1 ? lep1_vec : lep2_vec;
			/// lepton setup
			if 			( wLep==0 ) lep_vec = lep1_vec;
			else if ( wLep==1 ) lep_vec = lep2_vec;
			else if ( wLep==2 ) lep_vec = dilep_vec;
			else if ( wLep==3 ) lep_vec = lepc_vec;
			else 									return (float)(-1); // finished lepton setup
			/// return based on wDist variable
			float dr = lep_vec.DeltaR(jet_vec);
			float dphi = fabs(lep_vec.DeltaPhi(jet_vec));
			float deta = fabs(lep_vec.Eta()-jet_vec.Eta());
			if			( wDist==0 ) return dr;
			else if ( wDist==1 ) return dphi;
			else if ( wDist==2 ) return deta;
			else								 return (float)(-1);
		}
            return (float)(-1);
        };

	// deltaR functions - lep1
	auto cdRl1j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=0; wLep=0; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdRl2j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=0; wLep=1; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdRdlj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=0; wLep=2; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdRlcj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=0; wLep=3; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	// deltaPhi functions - lep1
	auto cdPhil1j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=1; wLep=0; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdPhil2j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=1; wLep=1; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdPhidlj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=1; wLep=2; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdPhilcj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=1; wLep=3; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	// deltaEta functions - lep1
	auto cdEtal1j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=2; wLep=0; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdEtal2j = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=2; wLep=1; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdEtadlj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=2; wLep=2; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};
	auto cdEtalcj = [&]( vector<float> LepPt, vector<float> LepPhi, vector<float> LepEta, vector<float> LepM, vector<float> JetPt, vector<float> JetPhi, vector<float> JetEta, vector<float> JetM ) {
		wDist=2; wLep=3; return cDistlj( LepPt, LepPhi, LepEta, LepM, JetPt, JetPhi, JetEta, JetM );
	};


	cout << " variables functions defined!" << endl;

	// creating new RDataFrame Object from the previous one and adding the new branches
	auto dfPP = df.Define( "mTlMET",   cmT,       { "LepPt", "LepEta", "LepPhi", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "mt2",      cmt2,      { "LepPt", "LepEta", "LepPhi", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "meff",     cmeff,     { "LepPt", "JetPt", "eT_miss" } )
	              .Define( "METmeff",  cMETmeff,  { "LepPt", "JetPt", "eT_miss" } )
	              .Define( "mlj",      cmlj,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "lumi", 		 cLumi,			{ "TotalWeight" } )
                      .Define( "mTH", 	 	 cmTH,		      { "LepPt", "LepEta", "LepPhi", "LepM", "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } )
                      .Define( "mTW",	  	 cmTW,		      { "LepPt", "LepEta", "LepPhi", "LepM", "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } )
                      .Define( "mt2tag", 	 cmt2tag,	      { "LepPt", "LepEta", "LepPhi", "LepM", "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } )
                      .Define( "mTl2",  	 cmTl2,		      { "LepPt", "LepEta", "LepPhi", "LepM", "eT_miss", "EtMissPhi" } )
                      .Define( "mTlmax",   cmTlmax,       { "LepPt", "LepEta", "LepPhi", "LepM", "eT_miss", "EtMissPhi" } )
                      .Define( "mTlmin",   cmTlmin,       { "LepPt", "LepEta", "LepPhi", "LepM", "eT_miss", "EtMissPhi" } )
                      .Define( "METrel",   cMETrel,       { "LepPhi", "JetPhi", "eT_miss", "EtMissPhi" } )
                      .Define( "mll",   	 cmll,	        { "LepPt", "LepEta", "LepPhi", "LepM" } )
                      .Define( "distmZ",   cdistmZ,	      { "LepPt", "LepEta", "LepPhi", "LepM" } )
	              .Define( "HT",       cHT,       { "JetPt" } )
//	              .Define( "HT1",      cHT1,      { "JetPt" } )
//	              .Define( "HT2",      cHT2,      { "JetPt" } )
//	              .Define( "HT3",      cHT3,      { "JetPt" } )
	              .Define( "LT",       cLT,       { "LepPt" } )
	              .Define( "mjj",      cmjj,      { "JetPt", "JetPhi", "JetEta", "JetM" } )
                      // DeltaAngll variables
                      .Define( "dRll",     cdRll,     { "LepPt", "LepPhi", "LepEta", "LepM" } )
	              .Define( "dPhill",   cdPhill,   { "LepPt", "LepPhi", "LepEta", "LepM" } )
	              .Define( "dEtall",   cdEtall,   { "LepPt", "LepPhi", "LepEta", "LepM" } )
                      // DeltaAngl1j variables
                      .Define( "dRl1j",  			cdRl1j,        { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dPhil1j",  		cdPhil1j,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dEtal1j",  		cdEtal1j,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      // DeltaAngl2j variables
                      .Define( "dRl2j",  			cdRl2j,        { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dPhil2j",  		cdPhil2j,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dEtal2j",  		cdEtal2j,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      // DeltaAngdlj variables
                      .Define( "dRdlj",  			cdRdlj,        { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dPhidlj",  		cdPhidlj,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dEtadlj",  		cdEtadlj,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      // DeltaAnglcj variables
                      .Define( "dRlcj",  			cdRlcj,        { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dPhilcj",  		cdPhilcj,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      .Define( "dEtalcj",  		cdEtalcj,      { "LepPt", "LepPhi", "LepEta", "LepM", "JetPt", "JetPhi", "JetEta", "JetM" } )
                      // DeltaAngl1MET
	              .Define( "dRl1MET",			cdRl1MET,   	   { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dPhil1MET",		cdPhil1MET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dEtal1MET",		cdEtal1MET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
								// DeltaAngl2MET
	              .Define( "dRl2MET",			cdRl2MET,   	   { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dPhil2MET",		cdPhil2MET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dEtal2MET",		cdEtal2MET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
								// DeltaAngdlMET
	              .Define( "dRdlMET",			cdRdlMET,   	   { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dPhidlMET",		cdPhidlMET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
	              .Define( "dEtadlMET",		cdEtadlMET,      { "LepPt", "LepPhi", "LepEta", "LepM", "eT_miss", "EtMissPhi" } )
								// DeltaAngjMET
	              .Define( "dRjMET",			cdRjMET,			   	 { "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } )
	              .Define( "dPhijMET",		cdPhijMET,		   	 { "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } )
	              .Define( "dEtajMET",		cdEtajMET,		   	 { "JetPt", "JetPhi", "JetEta", "JetM", "eT_miss", "EtMissPhi" } );

	cout << "   new branches added" << endl;

	/// getting all the column names and adding them to a vector
	auto colNames = df.GetColumnNames();
	vector< std::string > vars;
	for (auto &&colName : colNames ) {
		if (colName.find(".") == std::string::npos)
			vars.push_back(colName);
	}

	/// adding new variables names to vector
	vars.push_back( "mTlMET" );
	vars.push_back( "mt2" );
	vars.push_back( "meff" );
	vars.push_back( "METmeff" );
	vars.push_back( "mlj" );
	vars.push_back( "lumi" );
	vars.push_back( "mTl2" );
	vars.push_back( "mTH" );
	vars.push_back( "mTW" );
	vars.push_back( "mt2tag" );
	vars.push_back( "mTlmax" );
	vars.push_back( "mTlmin" );
	vars.push_back( "METrel" );
	vars.push_back( "mll" );
	vars.push_back( "distmZ" );
	vars.push_back( "HT" );
//	vars.push_back( "HT1" );
//	vars.push_back( "HT2" );
//	vars.push_back( "HT3" );
	vars.push_back( "LT" );
	vars.push_back( "mjj" );
	// DeltaAngll variables
	vars.push_back( "dRll" );
	vars.push_back( "dPhill" );
	vars.push_back( "dEtall" );
	// DeltaAngl1j variables
	vars.push_back( "dRl1j" );
	vars.push_back( "dPhil1j" );
	vars.push_back( "dEtal1j" );
	// DeltaAngl2j variables
	vars.push_back( "dRl2j" );
	vars.push_back( "dPhil2j" );
	vars.push_back( "dEtal2j" );
	// DeltaAngdlj variables
	vars.push_back( "dRdlj" );
	vars.push_back( "dPhidlj" );
	vars.push_back( "dEtadlj" );
	// DeltaAnglcj variables
	vars.push_back( "dRlcj" );
	vars.push_back( "dPhilcj" );
	vars.push_back( "dEtalcj" );
	// DeltaAngl1MET variables
	vars.push_back( "dRl1MET" );
	vars.push_back( "dPhil1MET" );
	vars.push_back( "dEtal1MET" );
	// DeltaAngl2MET variables
	vars.push_back( "dRl2MET" );
	vars.push_back( "dPhil2MET" );
	vars.push_back( "dEtal2MET" );
	// DeltaAngdlMET variables
	vars.push_back( "dRdlMET" );
	vars.push_back( "dPhidlMET" );
	vars.push_back( "dEtadlMET" );
	// DeltaAngjMET variables
	vars.push_back( "dRjMET" );
	vars.push_back( "dPhijMET" );
	vars.push_back( "dEtajMET" );

	// creating a new file and tree with the old and new branches
	dfPP.Snapshot( treeN.c_str(), (outpath+fileN).c_str(), vars );
	
	cout << "outfile created" << endl;
}



int main( int argc, char *argv[] ) {

	std::string treeN, inpath, fileN, sample, outpath, varnames, camp;
	varnames = "variables.conf";
	sample = "B";
	outpath = ".";

	if ( argc<2 ) 
 		exit(-1);
	else {
		int skip = -1;
		for ( int i=1 ; i<argc ; i++ ) {
			if ( i == skip ) { continue; skip = -1; }
    	std::string argi = std::string(argv[i]);
    	if      ( argi == "-b"	||	argi == "--background" )	sample = "B";
    	else if	( argi == "-s"	||	argi == "--signal" )	sample = "S";
    	else if	( argi == "-d"	||	argi == "--data" )	sample = "D";
    	else if ( argi == "-t"	||	argi == "--treename" )  { treeN = std::string( argv[i+1] ); skip = i+1; }
    	else if ( argi == "-n"	||	argi == "--filename" )  { fileN = std::string( argv[i+1] ); skip = i+1; }
    	else if ( argi == "-o"	||	argi == "--outpath" )		{ outpath = std::string( argv[i+1] ); skip = i+1; }
    	else if ( argi == "-i"	||	argi == "--inpath" )		{ inpath = std::string( argv[i+1] ); skip = i+1; }
    	else if ( argi == "-v"	||	argi == "--varnames" )	{ varnames = std::string( argv[i+1] ); skip = i+1; }
    	else if ( argi == "-c"	||	argi == "--mccampaign" ){ camp = std::string( argv[i+1] ); skip = i+1; }
    	//else    files.push_back( argi );
    }
	}

	Execute( treeN, inpath, fileN, sample, outpath, varnames, camp );

	return 0;
} 
