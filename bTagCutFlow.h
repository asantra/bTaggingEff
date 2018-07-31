//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 31 10:47:18 2018 by ROOT version 5.34/36
// from TTree AnaNtup/AnaNtup
// found on file: user.oducu.14520042._000001.output.root
//////////////////////////////////////////////////////////

#ifndef bTagCutFlow_h
#define bTagCutFlow_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class bTagCutFlow {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t          HLT_e24_lhmedium_nod0_ivarloose;
   Bool_t          HLT_e24_lhtight_nod0_ivarloose;
   Bool_t          HLT_e24_lhmedium_nod0_L1EM20VH;
   Bool_t          HLT_e26_lhtight_iloose;
   Bool_t          HLT_e26_lhtight_ivarloose;
   Bool_t          HLT_e26_lhtight_nod0_iloose;
   Bool_t          HLT_e26_lhtight_nod0_ivarloose;
   Bool_t          HLT_e60_lhmedium;
   Bool_t          HLT_e60_lhmedium_nod0;
   Bool_t          HLT_e120_lhloose_nod0;
   Bool_t          HLT_e140_lhloose_nod0;
   Bool_t          HLT_2e17_lhvloose_nod0;
   Bool_t          HLT_2e17_lhvloose_nod0_L12EM15VHI;
   Bool_t          HLT_2e15_lhvloose_nod0_L12EM13VH;
   Bool_t          HLT_2e24_lhvloose_nod0;
   Bool_t          HLT_e24_lhmedium_e9_lhmedium;
   Bool_t          HLT_e24_lhmedium_L1EM20VH;
   Bool_t          HLT_e24_lhmedium_iloose_L1EM20VH;
   Bool_t          HLT_e12_lhvloose_L12EM10VH;
   Bool_t          HLT_e17_lhloose_2e9_lhloose;
   Bool_t          HLT_mu24_ivarmedium;
   Bool_t          HLT_mu24_imedium;
   Bool_t          HLT_mu24_ivarloose;
   Bool_t          HLT_mu24_iloose;
   Bool_t          HLT_mu26_ivarmedium;
   Bool_t          HLT_mu20_ivarmedium_L1MU15;
   Bool_t          HLT_mu20_imedium_L1MU15;
   Bool_t          HLT_mu20_ivarloose_L1MU15;
   Bool_t          HLT_mu20_iloose_L1MU15;
   Bool_t          HLT_mu20_L1MU15;
   Bool_t          HLT_mu20_mu8noL1;
   Bool_t          HLT_mu22_mu8noL1;
   Bool_t          HLT_mu20_2mu4noL1;
   Bool_t          HLT_mu22_2mu4noL1;
   Bool_t          HLT_mu40;
   Bool_t          HLT_mu50;
   Bool_t          HLT_2mu10;
   Bool_t          HLT_2mu10_nomucomb;
   Bool_t          HLT_2mu14;
   Bool_t          HLT_2mu14_nomucomb;
   Bool_t          HLT_3mu6;
   Bool_t          HLT_3mu4;
   Bool_t          HLT_3mu6_msonly;
   Bool_t          HLT_xe100_L1XE50;
   Bool_t          HLT_xe80_mht_L1XE50;
   Bool_t          HLT_xe90_mht_L1XE50;
   Bool_t          HLT_xe90_pufit_L1XE50;
   Bool_t          HLT_xe100_mht_L1XE50;
   Bool_t          HLT_xe100_tc_em_L1XE50;
   Bool_t          HLT_xe100_pufit_L1XE55;
   Bool_t          HLT_xe100_pufit_L1XE50;
   Bool_t          HLT_xe110_mht_L1XE50;
   Bool_t          HLT_xe110_pufit_L1XE55;
   Bool_t          HLT_xe110_pufit_L1XE50;
   Bool_t          HLT_xe80_tc_lcw_L1XE50;
   Bool_t          HLT_xe90_tc_lcw_L1XE50;
   Bool_t          HLT_xe80_tc_lcw_wEFMu_L1XE50;
   Bool_t          HLT_e7_lhmedium_mu24;
   Bool_t          HLT_e7_lhmedium_nod0_mu24;
   Bool_t          HLT_e17_lhloose_nod0_mu14;
   Bool_t          HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;
   Bool_t          HLT_e26_lhmedium_nod0_mu8noL1;
   Bool_t          HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1;
   Bool_t          HLT_2e12_lhloose_L12EM10VH;
   Bool_t          HLT_e17_lhloose_mu14;
   Bool_t          HLT_mu18_mu8noL1;
   Bool_t          HLT_xe70;
   ULong64_t       EventNumber;
   Int_t           ChannelNumber;
   Double_t        AvgMu;
   Double_t        EventWeight;
   Double_t        PRWWeight;
   Double_t        TriggerDileptonSF;
   Int_t           bcid;
   Int_t           LB;
   Int_t           passGRL;
   Int_t           RunNb;
   UInt_t          PRWrandomRunNumber;
   Int_t           DetError;
   Int_t           Nvtx;
   Float_t         PV_z;
   Int_t           NMu;
   vector<double>  *Mu_eta;
   vector<double>  *Mu_phi;
   vector<double>  *Mu_pT;
   vector<int>     *Mu_charge;
   vector<double>  *Mu_sigd0;
   vector<double>  *Mu_d0pvtx;
   vector<double>  *Mu_d0pvtxerr;
   vector<double>  *Mu_z0sinTheta;
   vector<bool>    *Mu_isBad;
   vector<bool>    *Mu_isCosmic;
   vector<bool>    *Mu_passOR;
   vector<bool>    *Mu_isTight;
   vector<bool>    *Mu_isSig2017;
   vector<float>   *Mu_PromptLepTaggerBDTweight;
   vector<int>     *Mu_type;
   vector<int>     *Mu_origin;
   vector<double>  *Mu_ptcone20;
   vector<double>  *Mu_ptcone20_TightTTVA_pt500;
   vector<double>  *Mu_ptcone20_TightTTVA_pt1000;
   vector<double>  *Mu_ptcone30;
   vector<double>  *Mu_ptcone40;
   vector<double>  *Mu_ptvarcone20;
   vector<double>  *Mu_ptvarcone30;
   vector<double>  *Mu_ptvarcone30_TightTTVA_pt500;
   vector<double>  *Mu_ptvarcone30_TightTTVA_pt1000;
   vector<double>  *Mu_ptvarcone40;
   vector<double>  *Mu_topoetcone20;
   vector<double>  *Mu_topoetcone30;
   vector<double>  *Mu_topoetcone40;
   vector<double>  *Mu_neflowisol20;
   vector<bool>    *Mu_passIsoLooseTO;
   vector<bool>    *Mu_passIsoLoose;
   vector<bool>    *Mu_passIsoTight;
   vector<bool>    *Mu_passIsoGrad;
   vector<bool>    *Mu_passIsoGradCustomTight;
   vector<bool>    *Mu_passIsoGradCustom;
   vector<bool>    *Mu_passIsoGradLoose;
   vector<double>  *Mu_SFw;
   vector<float>   *Mu_IsoSFw;
   vector<float>   *Mu_StatUncReco;
   vector<float>   *Mu_SystUncReco;
   vector<float>   *Mu_StatUncReco_LOWPT;
   vector<float>   *Mu_SystUncReco_LOWPT;
   vector<float>   *Mu_StatUncISO;
   vector<float>   *Mu_SystUncISO;
   vector<bool>    *Mu_trigMatch_mu26_ivarmedium;
   vector<bool>    *Mu_trigMatch_mu20_iloose_L1MU15;
   vector<bool>    *Mu_trigMatch_mu40;
   vector<bool>    *Mu_trigMatch_mu50;
   vector<bool>    *Mu_trigMatch_mu8noL1;
   vector<bool>    *Mu_trigMatch_mu14;
   vector<bool>    *Mu_trigMatch_mu18;
   vector<bool>    *Mu_trigMatch_mu18_mu8noL1;
   vector<bool>    *Mu_trigMatch_e17_lhloose_mu14;
   vector<bool>    *Mu_trigMatch_e17_lhloose_nod0_mu14;
   vector<bool>    *Mu_trigMatch_mu20_mu8noL1;
   vector<bool>    *Mu_trigMatch_mu22_mu8noL1;
   vector<bool>    *Mu_trigMatch_mu24_iloose;
   vector<bool>    *Mu_trigMatch_mu24_ivarloose;
   vector<bool>    *Mu_trigMatch_mu24_iloose_L1MU15;
   vector<bool>    *Mu_trigMatch_mu24_ivarloose_L1MU15;
   vector<vector<bool> > *Mu_trigMatchPair_mu18_mu8noL1;
   vector<vector<bool> > *Mu_trigMatchPair_mu20_mu8noL1;
   vector<vector<bool> > *Mu_trigMatchPair_mu22_mu8noL1;
   Int_t           NEl;
   vector<double>  *El_eta;
   vector<double>  *El_etaclus;
   vector<double>  *El_phi;
   vector<double>  *El_pT;
   vector<double>  *El_E;
   vector<int>     *El_charge;
   vector<double>  *El_sigd0;
   vector<double>  *El_d0pvtx;
   vector<double>  *El_d0pvtxerr;
   vector<double>  *El_z0sinTheta;
   vector<bool>    *El_isLooseAndBLayerLH_baseline;
   vector<bool>    *El_isLooseAndBLayerLH_fromTool;
   vector<bool>    *El_isMediumLH;
   vector<bool>    *El_isTightLH;
   vector<bool>    *El_isSigNoCFT2017;
   vector<int>     *El_nBLayerHits;
   vector<int>     *El_expectBLayerHit;
   vector<bool>    *El_passOR;
   vector<float>   *El_passChargeFlipTaggerBDTmedium;
   vector<float>   *El_passChargeFlipTaggerBDTloose;
   vector<float>   *El_PromptLepTaggerBDTweight;
   vector<int>     *El_truthType;
   vector<int>     *El_truthOrigin;
   vector<int>     *El_truthPdgId;
   vector<int>     *El_bkgTruthType;
   vector<int>     *El_bkgTruthOrigin;
   vector<int>     *El_bkgMotherPdgId;
   vector<int>     *El_firstEgMotherTruthType;
   vector<int>     *El_firstEgMotherTruthOrigin;
   vector<int>     *El_firstEgMotherPdgId;
   vector<int>     *El_lastEgMotherTruthType;
   vector<int>     *El_lastEgMotherTruthOrigin;
   vector<int>     *El_lastEgMotherPdgId;
   vector<int>     *El_chFlip;
   vector<double>  *El_ptcone20;
   vector<double>  *El_ptcone20_TightTTVA_pt500;
   vector<double>  *El_ptvarcone20_TightTTVA_pt1000;
   vector<double>  *El_ptcone30;
   vector<double>  *El_ptcone40;
   vector<double>  *El_ptvarcone20;
   //vector<double>  *El_ptvarcone20_TightTTVA_pt1000;
   vector<double>  *El_ptvarcone30;
   vector<double>  *El_ptvarcone30_TightTTVA_pt1000;
   vector<double>  *El_ptvarcone40;
   vector<double>  *El_topoetcone20;
   vector<double>  *El_topoetcone30;
   vector<double>  *El_topoetcone40;
   vector<double>  *El_neflowisol20;
   vector<bool>    *El_passIsoLooseTO;
   vector<bool>    *El_passIsoLoose;
   vector<bool>    *El_passIsoTight;
   vector<bool>    *El_passIsoGrad;
   vector<bool>    *El_passIsoGradCustomTight;
   vector<bool>    *El_passIsoGradCustom;
   vector<bool>    *El_passIsoGradLoose;
   vector<double>  *El_SFwUncReco;
   vector<double>  *El_SFwTightLH;
   vector<double>  *El_SFwMediumLH;
   vector<double>  *El_SFwUncMediumLH;
   vector<double>  *El_SFwLooseAndBLayerLH;
   vector<double>  *El_SFweightCFT;
   vector<double>  *El_SFUncweightCFT;
   vector<double>  *El_SFweightCFID;
   vector<double>  *El_SFStatweightCFID;
   vector<double>  *El_SFSystweightCFID;
   vector<float>   *El_IsoSFwMediumLH;
   vector<float>   *El_IsoSFwUncMediumLH;
   vector<bool>    *El_trigMatch_e12_lhloose_L1EM10VH;
   vector<bool>    *El_trigMatch_e17_lhloose;
   vector<bool>    *El_trigMatch_e24_lhmedium_L1EM20VH;
   vector<bool>    *El_trigMatch_e24_lhmedium_iloose_L1EM20VH;
   vector<bool>    *El_trigMatch_e24_lhmedium_nod0_ivarloose;
   vector<bool>    *El_trigMatch_e24_lhtight_nod0_ivarloose;
   vector<bool>    *El_trigMatch_e26_lhtight_nod0_ivarloose;
   vector<bool>    *El_trigMatch_e60_lhmedium;
   vector<bool>    *El_trigMatch_e60_lhmedium_nod0;
   vector<bool>    *El_trigMatch_2e12_lhloose_L12EM10VH;
   vector<bool>    *El_trigMatch_2e15_lhloose_L12EM10VH;
   vector<bool>    *El_trigMatch_2e15_lhvloose_L12EM13VH;
   vector<bool>    *El_trigMatch_2e15_lhvloose_nod0_L12EM13VH;
   vector<bool>    *El_trigMatch_2e17_lhvloose_nod0;
   vector<bool>    *El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI;
   vector<bool>    *El_TrigMatch_2e24_lhvloose_nod0;
   vector<bool>    *El_trigMatch_e17_lhloose_mu14;
   vector<bool>    *El_trigMatch_e17_lhloose_nod0_mu14;
   Int_t           NJet;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_pT;
   vector<double>  *Jet_E;
   vector<int>     *Jet_nTrk;
   vector<double>  *Jet_quality;
   vector<float>   *Jet_EMFrac;
   vector<float>   *Jet_HECFrac;
   vector<double>  *Jet_JVT;
   vector<double>  *Jet_JVTsf;
   Float_t         totalJVTsf;
   vector<double>  *Jet_MV2c10;
   vector<double>  *Jet_SFw;
   vector<bool>    *Jet_passOR;
   vector<int>     *Jet_ConeTruthLabel;
   vector<int>     *Jet_PartonTruthLabel;
   vector<int>     *Jet_HadronConeExclTruthLabel;
   vector<double>  *Jet_deltaR;
   Float_t         Etmiss_TST_Etx;
   Float_t         Etmiss_TST_Ety;
   Float_t         Etmiss_TST_Et;
   Float_t         Etmiss_Truth_Etx;
   Float_t         Etmiss_Truth_Ety;
   Float_t         Etmiss_Truth_Et;
   Int_t           NTruthAntiktJet;
   vector<double>  *TruthAntiktJet_eta;
   vector<double>  *TruthAntiktJet_phi;
   vector<double>  *TruthAntiktJet_pT;
   vector<double>  *TruthAntiktJet_E;
   vector<int>     *TruthAntiktJet_ConeTruthLabel;
   vector<int>     *TruthAntiktJet_HadronConeExclTruthLabelID;
   vector<int>     *TruthAntiktJet_PartonTruthLabel;
   vector<int>     *TruthAntiktJetJet_ClassHF;
   Int_t           extraB;
   Int_t           extraC;
   Int_t           NTruthJet;
   vector<double>  *TruthJet_eta;
   vector<double>  *TruthJet_phi;
   vector<double>  *TruthJet_pT;
   vector<double>  *TruthJet_E;
   vector<double>  *TruthJet_id;
   vector<double>  *TruthJet_origin;
   Int_t           NTruthRealL;
   vector<double>  *TruthRealL_eta;
   vector<double>  *TruthRealL_phi;
   vector<double>  *TruthRealL_pT;
   vector<int>     *TruthRealL_id;
   vector<int>     *TruthRealL_origin;
   Int_t           SUSY_Spart_pdgId1;
   Int_t           SUSY_Spart_pdgId2;
   Int_t           SUSY_Gluino_decay1;
   Int_t           SUSY_Gluino_decay2;
   Float_t         GenFiltHT;
   Float_t         GenFiltMET;
   Float_t         TruthX1;
   Float_t         TruthX2;
   Float_t         TruthQ;
   Int_t           TruthPDGID1;
   Int_t           TruthPDGID2;
   Float_t         SherpaNjetWeight;

   // List of branches
   TBranch        *b_HLT_e24_lhmedium_nod0_ivarloose;   //!
   TBranch        *b_HLT_e24_lhtight_nod0_ivarloose;   //!
   TBranch        *b_HLT_e24_lhmedium_nod0_L1EM20VH;   //!
   TBranch        *b_HLT_e26_lhtight_iloose;   //!
   TBranch        *b_HLT_e26_lhtight_ivarloose;   //!
   TBranch        *b_HLT_e26_lhtight_nod0_iloose;   //!
   TBranch        *b_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_HLT_e60_lhmedium;   //!
   TBranch        *b_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_HLT_e120_lhloose_nod0;   //!
   TBranch        *b_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_HLT_2e17_lhvloose_nod0;   //!
   TBranch        *b_HLT_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_HLT_2e15_lhvloose_nod0_L12EM13VH;   //!
   TBranch        *b_HLT_2e24_lhvloose_nod0;   //!
   TBranch        *b_HLT_e24_lhmedium_e9_lhmedium;   //!
   TBranch        *b_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_HLT_e24_lhmedium_iloose_L1EM20VH;   //!
   TBranch        *b_HLT_e12_lhvloose_L12EM10VH;   //!
   TBranch        *b_HLT_e17_lhloose_2e9_lhloose;   //!
   TBranch        *b_HLT_mu24_ivarmedium;   //!
   TBranch        *b_HLT_mu24_imedium;   //!
   TBranch        *b_HLT_mu24_ivarloose;   //!
   TBranch        *b_HLT_mu24_iloose;   //!
   TBranch        *b_HLT_mu26_ivarmedium;   //!
   TBranch        *b_HLT_mu20_ivarmedium_L1MU15;   //!
   TBranch        *b_HLT_mu20_imedium_L1MU15;   //!
   TBranch        *b_HLT_mu20_ivarloose_L1MU15;   //!
   TBranch        *b_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_HLT_mu20_L1MU15;   //!
   TBranch        *b_HLT_mu20_mu8noL1;   //!
   TBranch        *b_HLT_mu22_mu8noL1;   //!
   TBranch        *b_HLT_mu20_2mu4noL1;   //!
   TBranch        *b_HLT_mu22_2mu4noL1;   //!
   TBranch        *b_HLT_mu40;   //!
   TBranch        *b_HLT_mu50;   //!
   TBranch        *b_HLT_2mu10;   //!
   TBranch        *b_HLT_2mu10_nomucomb;   //!
   TBranch        *b_HLT_2mu14;   //!
   TBranch        *b_HLT_2mu14_nomucomb;   //!
   TBranch        *b_HLT_3mu6;   //!
   TBranch        *b_HLT_3mu4;   //!
   TBranch        *b_HLT_3mu6_msonly;   //!
   TBranch        *b_HLT_xe100_L1XE50;   //!
   TBranch        *b_HLT_xe80_mht_L1XE50;   //!
   TBranch        *b_HLT_xe90_mht_L1XE50;   //!
   TBranch        *b_HLT_xe90_pufit_L1XE50;   //!
   TBranch        *b_HLT_xe100_mht_L1XE50;   //!
   TBranch        *b_HLT_xe100_tc_em_L1XE50;   //!
   TBranch        *b_HLT_xe100_pufit_L1XE55;   //!
   TBranch        *b_HLT_xe100_pufit_L1XE50;   //!
   TBranch        *b_HLT_xe110_mht_L1XE50;   //!
   TBranch        *b_HLT_xe110_pufit_L1XE55;   //!
   TBranch        *b_HLT_xe110_pufit_L1XE50;   //!
   TBranch        *b_HLT_xe80_tc_lcw_L1XE50;   //!
   TBranch        *b_HLT_xe90_tc_lcw_L1XE50;   //!
   TBranch        *b_HLT_xe80_tc_lcw_wEFMu_L1XE50;   //!
   TBranch        *b_HLT_e7_lhmedium_mu24;   //!
   TBranch        *b_HLT_e7_lhmedium_nod0_mu24;   //!
   TBranch        *b_HLT_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1;   //!
   TBranch        *b_HLT_e26_lhmedium_nod0_mu8noL1;   //!
   TBranch        *b_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1;   //!
   TBranch        *b_HLT_2e12_lhloose_L12EM10VH;   //!
   TBranch        *b_HLT_e17_lhloose_mu14;   //!
   TBranch        *b_HLT_mu18_mu8noL1;   //!
   TBranch        *b_HLT_xe70;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_ChannelNumber;   //!
   TBranch        *b_AvgMu;   //!
   TBranch        *b_EventWeight;   //!
   TBranch        *b_PRWWeight;   //!
   TBranch        *b_TriggerDileptonSF;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_LB;   //!
   TBranch        *b_passGRL;   //!
   TBranch        *b_RunNb;   //!
   TBranch        *b_PRWrandomRunNumber;   //!
   TBranch        *b_DetError;   //!
   TBranch        *b_Nvtx;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_NMu;   //!
   TBranch        *b_Mu_eta;   //!
   TBranch        *b_Mu_phi;   //!
   TBranch        *b_Mu_pT;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_Mu_sigd0;   //!
   TBranch        *b_Mu_d0pvtx;   //!
   TBranch        *b_Mu_d0pvtxerr;   //!
   TBranch        *b_Mu_z0sinTheta;   //!
   TBranch        *b_Mu_isBad;   //!
   TBranch        *b_Mu_isCosmic;   //!
   TBranch        *b_Mu_passOR;   //!
   TBranch        *b_Mu_isTight;   //!
   TBranch        *b_Mu_isSig2017;   //!
   TBranch        *b_Mu_PromptLepTaggerBDTweight;   //!
   TBranch        *b_Mu_type;   //!
   TBranch        *b_Mu_origin;   //!
   TBranch        *b_Mu_ptcone20;   //!
   TBranch        *b_Mu_ptcone20_TightTTVA_pt500;   //!
   TBranch        *b_Mu_ptcone20_TightTTVA_pt1000;   //!
   TBranch        *b_Mu_ptcone30;   //!
   TBranch        *b_Mu_ptcone40;   //!
   TBranch        *b_Mu_ptvarcone20;   //!
   TBranch        *b_Mu_ptvarcone30;   //!
   TBranch        *b_Mu_ptvarcone30_TightTTVA_pt500;   //!
   TBranch        *b_Mu_ptvarcone30_TightTTVA_pt1000;   //!
   TBranch        *b_Mu_ptvarcone40;   //!
   TBranch        *b_Mu_topoetcone20;   //!
   TBranch        *b_Mu_topoetcone30;   //!
   TBranch        *b_Mu_topoetcone40;   //!
   TBranch        *b_Mu_neflowisol20;   //!
   TBranch        *b_Mu_passIsoLooseTO;   //!
   TBranch        *b_Mu_passIsoLoose;   //!
   TBranch        *b_Mu_passIsoTight;   //!
   TBranch        *b_Mu_passIsoGrad;   //!
   TBranch        *b_Mu_passIsoGradCustomTight;   //!
   TBranch        *b_Mu_passIsoGradCustom;   //!
   TBranch        *b_Mu_passIsoGradLoose;   //!
   TBranch        *b_Mu_SFw;   //!
   TBranch        *b_Mu_IsoSFw;   //!
   TBranch        *b_Mu_StatUncReco;   //!
   TBranch        *b_Mu_SystUncReco;   //!
   TBranch        *b_Mu_StatUncReco_LOWPT;   //!
   TBranch        *b_Mu_SystUncReco_LOWPT;   //!
   TBranch        *b_Mu_StatUncISO;   //!
   TBranch        *b_Mu_SystUncISO;   //!
   TBranch        *b_Mu_trigMatch_mu26_ivarmedium;   //!
   TBranch        *b_Mu_trigMatch_mu20_iloose_L1MU15;   //!
   TBranch        *b_Mu_trigMatch_mu40;   //!
   TBranch        *b_Mu_trigMatch_mu50;   //!
   TBranch        *b_Mu_trigMatch_mu8noL1;   //!
   TBranch        *b_Mu_trigMatch_mu14;   //!
   TBranch        *b_Mu_trigMatch_mu18;   //!
   TBranch        *b_Mu_trigMatch_mu18_mu8noL1;   //!
   TBranch        *b_Mu_trigMatch_e17_lhloose_mu14;   //!
   TBranch        *b_Mu_trigMatch_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_Mu_trigMatch_mu20_mu8noL1;   //!
   TBranch        *b_Mu_trigMatch_mu22_mu8noL1;   //!
   TBranch        *b_Mu_trigMatch_mu24_iloose;   //!
   TBranch        *b_Mu_trigMatch_mu24_ivarloose;   //!
   TBranch        *b_Mu_trigMatch_mu24_iloose_L1MU15;   //!
   TBranch        *b_Mu_trigMatch_mu24_ivarloose_L1MU15;   //!
   TBranch        *b_Mu_trigMatchPair_mu18_mu8noL1;   //!
   TBranch        *b_Mu_trigMatchPair_mu20_mu8noL1;   //!
   TBranch        *b_Mu_trigMatchPair_mu22_mu8noL1;   //!
   TBranch        *b_NEl;   //!
   TBranch        *b_El_eta;   //!
   TBranch        *b_El_etaclus;   //!
   TBranch        *b_El_phi;   //!
   TBranch        *b_El_pT;   //!
   TBranch        *b_El_E;   //!
   TBranch        *b_El_charge;   //!
   TBranch        *b_El_sigd0;   //!
   TBranch        *b_El_d0pvtx;   //!
   TBranch        *b_El_d0pvtxerr;   //!
   TBranch        *b_El_z0sinTheta;   //!
   TBranch        *b_El_isLooseAndBLayerLH_baseline;   //!
   TBranch        *b_El_isLooseAndBLayerLH_fromTool;   //!
   TBranch        *b_El_isMediumLH;   //!
   TBranch        *b_El_isTightLH;   //!
   TBranch        *b_El_isSigNoCFT2017;   //!
   TBranch        *b_El_nBLayerHits;   //!
   TBranch        *b_El_expectBLayerHit;   //!
   TBranch        *b_El_passOR;   //!
   TBranch        *b_El_passChargeFlipTaggerBDTmedium;   //!
   TBranch        *b_El_passChargeFlipTaggerBDTloose;   //!
   TBranch        *b_El_PromptLepTaggerBDTweight;   //!
   TBranch        *b_El_truthType;   //!
   TBranch        *b_El_truthOrigin;   //!
   TBranch        *b_El_truthPdgId;   //!
   TBranch        *b_El_bkgTruthType;   //!
   TBranch        *b_El_bkgTruthOrigin;   //!
   TBranch        *b_El_bkgMotherPdgId;   //!
   TBranch        *b_El_firstEgMotherTruthType;   //!
   TBranch        *b_El_firstEgMotherTruthOrigin;   //!
   TBranch        *b_El_firstEgMotherPdgId;   //!
   TBranch        *b_El_lastEgMotherTruthType;   //!
   TBranch        *b_El_lastEgMotherTruthOrigin;   //!
   TBranch        *b_El_lastEgMotherPdgId;   //!
   TBranch        *b_El_chFlip;   //!
   TBranch        *b_El_ptcone20;   //!
   TBranch        *b_El_ptcone20_TightTTVA_pt500;   //!
   TBranch        *b_El_ptvarcone20_TightTTVA_pt1000;   //!
   TBranch        *b_El_ptcone30;   //!
   TBranch        *b_El_ptcone40;   //!
   TBranch        *b_El_ptvarcone20;   //!
   //TBranch        *b_El_ptvarcone20_TightTTVA_pt1000;   //!
   TBranch        *b_El_ptvarcone30;   //!
   TBranch        *b_El_ptvarcone30_TightTTVA_pt1000;   //!
   TBranch        *b_El_ptvarcone40;   //!
   TBranch        *b_El_topoetcone20;   //!
   TBranch        *b_El_topoetcone30;   //!
   TBranch        *b_El_topoetcone40;   //!
   TBranch        *b_El_neflowisol20;   //!
   TBranch        *b_El_passIsoLooseTO;   //!
   TBranch        *b_El_passIsoLoose;   //!
   TBranch        *b_El_passIsoTight;   //!
   TBranch        *b_El_passIsoGrad;   //!
   TBranch        *b_El_passIsoGradCustomTight;   //!
   TBranch        *b_El_passIsoGradCustom;   //!
   TBranch        *b_El_passIsoGradLoose;   //!
   TBranch        *b_El_SFwUncReco;   //!
   TBranch        *b_El_SFwTightLH;   //!
   TBranch        *b_El_SFwMediumLH;   //!
   TBranch        *b_El_SFwUncMediumLH;   //!
   TBranch        *b_El_SFwLooseAndBLayerLH;   //!
   TBranch        *b_El_SFweightCFT;   //!
   TBranch        *b_El_SFUncweightCFT;   //!
   TBranch        *b_El_SFweightCFID;   //!
   TBranch        *b_El_SFStatweightCFID;   //!
   TBranch        *b_El_SFSystweightCFID;   //!
   TBranch        *b_El_IsoSFwMediumLH;   //!
   TBranch        *b_El_IsoSFwUncMediumLH;   //!
   TBranch        *b_El_trigMatch_e12_lhloose_L1EM10VH;   //!
   TBranch        *b_El_trigMatch_e17_lhloose;   //!
   TBranch        *b_El_trigMatch_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_El_trigMatch_e24_lhmedium_iloose_L1EM20VH;   //!
   TBranch        *b_El_trigMatch_e24_lhmedium_nod0_ivarloose;   //!
   TBranch        *b_El_trigMatch_e24_lhtight_nod0_ivarloose;   //!
   TBranch        *b_El_trigMatch_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_El_trigMatch_e60_lhmedium;   //!
   TBranch        *b_El_trigMatch_e60_lhmedium_nod0;   //!
   TBranch        *b_El_trigMatch_2e12_lhloose_L12EM10VH;   //!
   TBranch        *b_El_trigMatch_2e15_lhloose_L12EM10VH;   //!
   TBranch        *b_El_trigMatch_2e15_lhvloose_L12EM13VH;   //!
   TBranch        *b_El_trigMatch_2e15_lhvloose_nod0_L12EM13VH;   //!
   TBranch        *b_El_trigMatch_2e17_lhvloose_nod0;   //!
   TBranch        *b_El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_El_TrigMatch_2e24_lhvloose_nod0;   //!
   TBranch        *b_El_trigMatch_e17_lhloose_mu14;   //!
   TBranch        *b_El_trigMatch_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_NJet;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pT;   //!
   TBranch        *b_Jet_E;   //!
   TBranch        *b_Jet_nTrk;   //!
   TBranch        *b_Jet_quality;   //!
   TBranch        *b_Jet_EMFrac;   //!
   TBranch        *b_Jet_HECFrac;   //!
   TBranch        *b_Jet_JVT;   //!
   TBranch        *b_Jet_JVTsf;   //!
   TBranch        *b_totalJVTsf;   //!
   TBranch        *b_Jet_MV2c10;   //!
   TBranch        *b_Jet_SFw;   //!
   TBranch        *b_Jet_passOR;   //!
   TBranch        *b_Jet_ConeTruthLabel;   //!
   TBranch        *b_Jet_PartonTruthLabel;   //!
   TBranch        *b_Jet_HadronConeExclTruthLabel;   //!
   TBranch        *b_Jet_deltaR;   //!
   TBranch        *b_Etmiss_TST_Etx;   //!
   TBranch        *b_Etmiss_TST_Ety;   //!
   TBranch        *b_Etmiss_TST_Et;   //!
   TBranch        *b_Etmiss_Truth_Etx;   //!
   TBranch        *b_Etmiss_Truth_Ety;   //!
   TBranch        *b_Etmiss_Truth_Et;   //!
   TBranch        *b_NTruthAntiktJet;   //!
   TBranch        *b_TruthAntiktJet_eta;   //!
   TBranch        *b_TruthAntiktJet_phi;   //!
   TBranch        *b_TruthAntiktJet_pT;   //!
   TBranch        *b_TruthAntiktJet_E;   //!
   TBranch        *b_TruthAntiktJet_ConeTruthLabel;   //!
   TBranch        *b_TruthAntiktJet_HadronConeExclTruthLabelID;   //!
   TBranch        *b_TruthAntiktJet_PartonTruthLabel;   //!
   TBranch        *b_TruthAntiktJetJet_ClassHF;   //!
   TBranch        *b_extraB;   //!
   TBranch        *b_extraC;   //!
   TBranch        *b_NTruthJet;   //!
   TBranch        *b_TruthJet_eta;   //!
   TBranch        *b_TruthJet_phi;   //!
   TBranch        *b_TruthJet_pT;   //!
   TBranch        *b_TruthJet_E;   //!
   TBranch        *b_TruthJet_id;   //!
   TBranch        *b_TruthJet_origin;   //!
   TBranch        *b_NTruthRealL;   //!
   TBranch        *b_TruthRealL_eta;   //!
   TBranch        *b_TruthRealL_phi;   //!
   TBranch        *b_TruthRealL_pT;   //!
   TBranch        *b_TruthRealL_id;   //!
   TBranch        *b_TruthRealL_origin;   //!
   TBranch        *b_SUSY_Spart_pdgId1;   //!
   TBranch        *b_SUSY_Spart_pdgId2;   //!
   TBranch        *b_SUSY_Gluino_decay1;   //!
   TBranch        *b_SUSY_Gluino_decay2;   //!
   TBranch        *b_GenFiltHT;   //!
   TBranch        *b_GenFiltMET;   //!
   TBranch        *b_TruthX1;   //!
   TBranch        *b_TruthX2;   //!
   TBranch        *b_TruthQ;   //!
   TBranch        *b_TruthPDGID1;   //!
   TBranch        *b_TruthPDGID2;   //!
   TBranch        *b_SherpaNjetWeight;   //!

   bTagCutFlow(TTree *tree=0);
   virtual ~bTagCutFlow();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bTagCutFlow_cxx
bTagCutFlow::bTagCutFlow(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("user.oducu.14520042._000001.output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("user.oducu.14520042._000001.output.root");
      }
      f->GetObject("AnaNtup",tree);

   }
   Init(tree);
}

bTagCutFlow::~bTagCutFlow()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bTagCutFlow::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bTagCutFlow::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bTagCutFlow::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Mu_eta = 0;
   Mu_phi = 0;
   Mu_pT = 0;
   Mu_charge = 0;
   Mu_sigd0 = 0;
   Mu_d0pvtx = 0;
   Mu_d0pvtxerr = 0;
   Mu_z0sinTheta = 0;
   Mu_isBad = 0;
   Mu_isCosmic = 0;
   Mu_passOR = 0;
   Mu_isTight = 0;
   Mu_isSig2017 = 0;
   Mu_PromptLepTaggerBDTweight = 0;
   Mu_type = 0;
   Mu_origin = 0;
   Mu_ptcone20 = 0;
   Mu_ptcone20_TightTTVA_pt500 = 0;
   Mu_ptcone20_TightTTVA_pt1000 = 0;
   Mu_ptcone30 = 0;
   Mu_ptcone40 = 0;
   Mu_ptvarcone20 = 0;
   Mu_ptvarcone30 = 0;
   Mu_ptvarcone30_TightTTVA_pt500 = 0;
   Mu_ptvarcone30_TightTTVA_pt1000 = 0;
   Mu_ptvarcone40 = 0;
   Mu_topoetcone20 = 0;
   Mu_topoetcone30 = 0;
   Mu_topoetcone40 = 0;
   Mu_neflowisol20 = 0;
   Mu_passIsoLooseTO = 0;
   Mu_passIsoLoose = 0;
   Mu_passIsoTight = 0;
   Mu_passIsoGrad = 0;
   Mu_passIsoGradCustomTight = 0;
   Mu_passIsoGradCustom = 0;
   Mu_passIsoGradLoose = 0;
   Mu_SFw = 0;
   Mu_IsoSFw = 0;
   Mu_StatUncReco = 0;
   Mu_SystUncReco = 0;
   Mu_StatUncReco_LOWPT = 0;
   Mu_SystUncReco_LOWPT = 0;
   Mu_StatUncISO = 0;
   Mu_SystUncISO = 0;
   Mu_trigMatch_mu26_ivarmedium = 0;
   Mu_trigMatch_mu20_iloose_L1MU15 = 0;
   Mu_trigMatch_mu40 = 0;
   Mu_trigMatch_mu50 = 0;
   Mu_trigMatch_mu8noL1 = 0;
   Mu_trigMatch_mu14 = 0;
   Mu_trigMatch_mu18 = 0;
   Mu_trigMatch_mu18_mu8noL1 = 0;
   Mu_trigMatch_e17_lhloose_mu14 = 0;
   Mu_trigMatch_e17_lhloose_nod0_mu14 = 0;
   Mu_trigMatch_mu20_mu8noL1 = 0;
   Mu_trigMatch_mu22_mu8noL1 = 0;
   Mu_trigMatch_mu24_iloose = 0;
   Mu_trigMatch_mu24_ivarloose = 0;
   Mu_trigMatch_mu24_iloose_L1MU15 = 0;
   Mu_trigMatch_mu24_ivarloose_L1MU15 = 0;
   Mu_trigMatchPair_mu18_mu8noL1 = 0;
   Mu_trigMatchPair_mu20_mu8noL1 = 0;
   Mu_trigMatchPair_mu22_mu8noL1 = 0;
   El_eta = 0;
   El_etaclus = 0;
   El_phi = 0;
   El_pT = 0;
   El_E = 0;
   El_charge = 0;
   El_sigd0 = 0;
   El_d0pvtx = 0;
   El_d0pvtxerr = 0;
   El_z0sinTheta = 0;
   El_isLooseAndBLayerLH_baseline = 0;
   El_isLooseAndBLayerLH_fromTool = 0;
   El_isMediumLH = 0;
   El_isTightLH = 0;
   El_isSigNoCFT2017 = 0;
   El_nBLayerHits = 0;
   El_expectBLayerHit = 0;
   El_passOR = 0;
   El_passChargeFlipTaggerBDTmedium = 0;
   El_passChargeFlipTaggerBDTloose = 0;
   El_PromptLepTaggerBDTweight = 0;
   El_truthType = 0;
   El_truthOrigin = 0;
   El_truthPdgId = 0;
   El_bkgTruthType = 0;
   El_bkgTruthOrigin = 0;
   El_bkgMotherPdgId = 0;
   El_firstEgMotherTruthType = 0;
   El_firstEgMotherTruthOrigin = 0;
   El_firstEgMotherPdgId = 0;
   El_lastEgMotherTruthType = 0;
   El_lastEgMotherTruthOrigin = 0;
   El_lastEgMotherPdgId = 0;
   El_chFlip = 0;
   El_ptcone20 = 0;
   El_ptcone20_TightTTVA_pt500 = 0;
   El_ptvarcone20_TightTTVA_pt1000 = 0;
   El_ptcone30 = 0;
   El_ptcone40 = 0;
   El_ptvarcone20 = 0;
   El_ptvarcone20_TightTTVA_pt1000 = 0;
   El_ptvarcone30 = 0;
   El_ptvarcone30_TightTTVA_pt1000 = 0;
   El_ptvarcone40 = 0;
   El_topoetcone20 = 0;
   El_topoetcone30 = 0;
   El_topoetcone40 = 0;
   El_neflowisol20 = 0;
   El_passIsoLooseTO = 0;
   El_passIsoLoose = 0;
   El_passIsoTight = 0;
   El_passIsoGrad = 0;
   El_passIsoGradCustomTight = 0;
   El_passIsoGradCustom = 0;
   El_passIsoGradLoose = 0;
   El_SFwUncReco = 0;
   El_SFwTightLH = 0;
   El_SFwMediumLH = 0;
   El_SFwUncMediumLH = 0;
   El_SFwLooseAndBLayerLH = 0;
   El_SFweightCFT = 0;
   El_SFUncweightCFT = 0;
   El_SFweightCFID = 0;
   El_SFStatweightCFID = 0;
   El_SFSystweightCFID = 0;
   El_IsoSFwMediumLH = 0;
   El_IsoSFwUncMediumLH = 0;
   El_trigMatch_e12_lhloose_L1EM10VH = 0;
   El_trigMatch_e17_lhloose = 0;
   El_trigMatch_e24_lhmedium_L1EM20VH = 0;
   El_trigMatch_e24_lhmedium_iloose_L1EM20VH = 0;
   El_trigMatch_e24_lhmedium_nod0_ivarloose = 0;
   El_trigMatch_e24_lhtight_nod0_ivarloose = 0;
   El_trigMatch_e26_lhtight_nod0_ivarloose = 0;
   El_trigMatch_e60_lhmedium = 0;
   El_trigMatch_e60_lhmedium_nod0 = 0;
   El_trigMatch_2e12_lhloose_L12EM10VH = 0;
   El_trigMatch_2e15_lhloose_L12EM10VH = 0;
   El_trigMatch_2e15_lhvloose_L12EM13VH = 0;
   El_trigMatch_2e15_lhvloose_nod0_L12EM13VH = 0;
   El_trigMatch_2e17_lhvloose_nod0 = 0;
   El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI = 0;
   El_TrigMatch_2e24_lhvloose_nod0 = 0;
   El_trigMatch_e17_lhloose_mu14 = 0;
   El_trigMatch_e17_lhloose_nod0_mu14 = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_pT = 0;
   Jet_E = 0;
   Jet_nTrk = 0;
   Jet_quality = 0;
   Jet_EMFrac = 0;
   Jet_HECFrac = 0;
   Jet_JVT = 0;
   Jet_JVTsf = 0;
   Jet_MV2c10 = 0;
   Jet_SFw = 0;
   Jet_passOR = 0;
   Jet_ConeTruthLabel = 0;
   Jet_PartonTruthLabel = 0;
   Jet_HadronConeExclTruthLabel = 0;
   Jet_deltaR = 0;
   TruthAntiktJet_eta = 0;
   TruthAntiktJet_phi = 0;
   TruthAntiktJet_pT = 0;
   TruthAntiktJet_E = 0;
   TruthAntiktJet_ConeTruthLabel = 0;
   TruthAntiktJet_HadronConeExclTruthLabelID = 0;
   TruthAntiktJet_PartonTruthLabel = 0;
   TruthAntiktJetJet_ClassHF = 0;
   TruthJet_eta = 0;
   TruthJet_phi = 0;
   TruthJet_pT = 0;
   TruthJet_E = 0;
   TruthJet_id = 0;
   TruthJet_origin = 0;
   TruthRealL_eta = 0;
   TruthRealL_phi = 0;
   TruthRealL_pT = 0;
   TruthRealL_id = 0;
   TruthRealL_origin = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_ivarloose", &HLT_e24_lhmedium_nod0_ivarloose, &b_HLT_e24_lhmedium_nod0_ivarloose);
   fChain->SetBranchAddress("HLT_e24_lhtight_nod0_ivarloose", &HLT_e24_lhtight_nod0_ivarloose, &b_HLT_e24_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_L1EM20VH", &HLT_e24_lhmedium_nod0_L1EM20VH, &b_HLT_e24_lhmedium_nod0_L1EM20VH);
   fChain->SetBranchAddress("HLT_e26_lhtight_iloose", &HLT_e26_lhtight_iloose, &b_HLT_e26_lhtight_iloose);
   fChain->SetBranchAddress("HLT_e26_lhtight_ivarloose", &HLT_e26_lhtight_ivarloose, &b_HLT_e26_lhtight_ivarloose);
   fChain->SetBranchAddress("HLT_e26_lhtight_nod0_iloose", &HLT_e26_lhtight_nod0_iloose, &b_HLT_e26_lhtight_nod0_iloose);
   fChain->SetBranchAddress("HLT_e26_lhtight_nod0_ivarloose", &HLT_e26_lhtight_nod0_ivarloose, &b_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("HLT_e60_lhmedium", &HLT_e60_lhmedium, &b_HLT_e60_lhmedium);
   fChain->SetBranchAddress("HLT_e60_lhmedium_nod0", &HLT_e60_lhmedium_nod0, &b_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("HLT_e120_lhloose_nod0", &HLT_e120_lhloose_nod0, &b_HLT_e120_lhloose_nod0);
   fChain->SetBranchAddress("HLT_e140_lhloose_nod0", &HLT_e140_lhloose_nod0, &b_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("HLT_2e17_lhvloose_nod0", &HLT_2e17_lhvloose_nod0, &b_HLT_2e17_lhvloose_nod0);
   fChain->SetBranchAddress("HLT_2e17_lhvloose_nod0_L12EM15VHI", &HLT_2e17_lhvloose_nod0_L12EM15VHI, &b_HLT_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("HLT_2e15_lhvloose_nod0_L12EM13VH", &HLT_2e15_lhvloose_nod0_L12EM13VH, &b_HLT_2e15_lhvloose_nod0_L12EM13VH);
   fChain->SetBranchAddress("HLT_2e24_lhvloose_nod0", &HLT_2e24_lhvloose_nod0, &b_HLT_2e24_lhvloose_nod0);
   fChain->SetBranchAddress("HLT_e24_lhmedium_e9_lhmedium", &HLT_e24_lhmedium_e9_lhmedium, &b_HLT_e24_lhmedium_e9_lhmedium);
   fChain->SetBranchAddress("HLT_e24_lhmedium_L1EM20VH", &HLT_e24_lhmedium_L1EM20VH, &b_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("HLT_e24_lhmedium_iloose_L1EM20VH", &HLT_e24_lhmedium_iloose_L1EM20VH, &b_HLT_e24_lhmedium_iloose_L1EM20VH);
   fChain->SetBranchAddress("HLT_e12_lhvloose_L12EM10VH", &HLT_e12_lhvloose_L12EM10VH, &b_HLT_e12_lhvloose_L12EM10VH);
   fChain->SetBranchAddress("HLT_e17_lhloose_2e9_lhloose", &HLT_e17_lhloose_2e9_lhloose, &b_HLT_e17_lhloose_2e9_lhloose);
   fChain->SetBranchAddress("HLT_mu24_ivarmedium", &HLT_mu24_ivarmedium, &b_HLT_mu24_ivarmedium);
   fChain->SetBranchAddress("HLT_mu24_imedium", &HLT_mu24_imedium, &b_HLT_mu24_imedium);
   fChain->SetBranchAddress("HLT_mu24_ivarloose", &HLT_mu24_ivarloose, &b_HLT_mu24_ivarloose);
   fChain->SetBranchAddress("HLT_mu24_iloose", &HLT_mu24_iloose, &b_HLT_mu24_iloose);
   fChain->SetBranchAddress("HLT_mu26_ivarmedium", &HLT_mu26_ivarmedium, &b_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("HLT_mu20_ivarmedium_L1MU15", &HLT_mu20_ivarmedium_L1MU15, &b_HLT_mu20_ivarmedium_L1MU15);
   fChain->SetBranchAddress("HLT_mu20_imedium_L1MU15", &HLT_mu20_imedium_L1MU15, &b_HLT_mu20_imedium_L1MU15);
   fChain->SetBranchAddress("HLT_mu20_ivarloose_L1MU15", &HLT_mu20_ivarloose_L1MU15, &b_HLT_mu20_ivarloose_L1MU15);
   fChain->SetBranchAddress("HLT_mu20_iloose_L1MU15", &HLT_mu20_iloose_L1MU15, &b_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("HLT_mu20_L1MU15", &HLT_mu20_L1MU15, &b_HLT_mu20_L1MU15);
   fChain->SetBranchAddress("HLT_mu20_mu8noL1", &HLT_mu20_mu8noL1, &b_HLT_mu20_mu8noL1);
   fChain->SetBranchAddress("HLT_mu22_mu8noL1", &HLT_mu22_mu8noL1, &b_HLT_mu22_mu8noL1);
   fChain->SetBranchAddress("HLT_mu20_2mu4noL1", &HLT_mu20_2mu4noL1, &b_HLT_mu20_2mu4noL1);
   fChain->SetBranchAddress("HLT_mu22_2mu4noL1", &HLT_mu22_2mu4noL1, &b_HLT_mu22_2mu4noL1);
   fChain->SetBranchAddress("HLT_mu40", &HLT_mu40, &b_HLT_mu40);
   fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
   fChain->SetBranchAddress("HLT_2mu10", &HLT_2mu10, &b_HLT_2mu10);
   fChain->SetBranchAddress("HLT_2mu10_nomucomb", &HLT_2mu10_nomucomb, &b_HLT_2mu10_nomucomb);
   fChain->SetBranchAddress("HLT_2mu14", &HLT_2mu14, &b_HLT_2mu14);
   fChain->SetBranchAddress("HLT_2mu14_nomucomb", &HLT_2mu14_nomucomb, &b_HLT_2mu14_nomucomb);
   fChain->SetBranchAddress("HLT_3mu6", &HLT_3mu6, &b_HLT_3mu6);
   fChain->SetBranchAddress("HLT_3mu4", &HLT_3mu4, &b_HLT_3mu4);
   fChain->SetBranchAddress("HLT_3mu6_msonly", &HLT_3mu6_msonly, &b_HLT_3mu6_msonly);
   fChain->SetBranchAddress("HLT_xe100_L1XE50", &HLT_xe100_L1XE50, &b_HLT_xe100_L1XE50);
   fChain->SetBranchAddress("HLT_xe80_mht_L1XE50", &HLT_xe80_mht_L1XE50, &b_HLT_xe80_mht_L1XE50);
   fChain->SetBranchAddress("HLT_xe90_mht_L1XE50", &HLT_xe90_mht_L1XE50, &b_HLT_xe90_mht_L1XE50);
   fChain->SetBranchAddress("HLT_xe90_pufit_L1XE50", &HLT_xe90_pufit_L1XE50, &b_HLT_xe90_pufit_L1XE50);
   fChain->SetBranchAddress("HLT_xe100_mht_L1XE50", &HLT_xe100_mht_L1XE50, &b_HLT_xe100_mht_L1XE50);
   fChain->SetBranchAddress("HLT_xe100_tc_em_L1XE50", &HLT_xe100_tc_em_L1XE50, &b_HLT_xe100_tc_em_L1XE50);
   fChain->SetBranchAddress("HLT_xe100_pufit_L1XE55", &HLT_xe100_pufit_L1XE55, &b_HLT_xe100_pufit_L1XE55);
   fChain->SetBranchAddress("HLT_xe100_pufit_L1XE50", &HLT_xe100_pufit_L1XE50, &b_HLT_xe100_pufit_L1XE50);
   fChain->SetBranchAddress("HLT_xe110_mht_L1XE50", &HLT_xe110_mht_L1XE50, &b_HLT_xe110_mht_L1XE50);
   fChain->SetBranchAddress("HLT_xe110_pufit_L1XE55", &HLT_xe110_pufit_L1XE55, &b_HLT_xe110_pufit_L1XE55);
   fChain->SetBranchAddress("HLT_xe110_pufit_L1XE50", &HLT_xe110_pufit_L1XE50, &b_HLT_xe110_pufit_L1XE50);
   fChain->SetBranchAddress("HLT_xe80_tc_lcw_L1XE50", &HLT_xe80_tc_lcw_L1XE50, &b_HLT_xe80_tc_lcw_L1XE50);
   fChain->SetBranchAddress("HLT_xe90_tc_lcw_L1XE50", &HLT_xe90_tc_lcw_L1XE50, &b_HLT_xe90_tc_lcw_L1XE50);
   fChain->SetBranchAddress("HLT_xe80_tc_lcw_wEFMu_L1XE50", &HLT_xe80_tc_lcw_wEFMu_L1XE50, &b_HLT_xe80_tc_lcw_wEFMu_L1XE50);
   fChain->SetBranchAddress("HLT_e7_lhmedium_mu24", &HLT_e7_lhmedium_mu24, &b_HLT_e7_lhmedium_mu24);
   fChain->SetBranchAddress("HLT_e7_lhmedium_nod0_mu24", &HLT_e7_lhmedium_nod0_mu24, &b_HLT_e7_lhmedium_nod0_mu24);
   fChain->SetBranchAddress("HLT_e17_lhloose_nod0_mu14", &HLT_e17_lhloose_nod0_mu14, &b_HLT_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1", &HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1, &b_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1);
   fChain->SetBranchAddress("HLT_e26_lhmedium_nod0_mu8noL1", &HLT_e26_lhmedium_nod0_mu8noL1, &b_HLT_e26_lhmedium_nod0_mu8noL1);
   fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", &HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1, &b_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1);
   fChain->SetBranchAddress("HLT_2e12_lhloose_L12EM10VH", &HLT_2e12_lhloose_L12EM10VH, &b_HLT_2e12_lhloose_L12EM10VH);
   fChain->SetBranchAddress("HLT_e17_lhloose_mu14", &HLT_e17_lhloose_mu14, &b_HLT_e17_lhloose_mu14);
   fChain->SetBranchAddress("HLT_mu18_mu8noL1", &HLT_mu18_mu8noL1, &b_HLT_mu18_mu8noL1);
   fChain->SetBranchAddress("HLT_xe70", &HLT_xe70, &b_HLT_xe70);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("ChannelNumber", &ChannelNumber, &b_ChannelNumber);
   fChain->SetBranchAddress("AvgMu", &AvgMu, &b_AvgMu);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
   fChain->SetBranchAddress("PRWWeight", &PRWWeight, &b_PRWWeight);
   fChain->SetBranchAddress("TriggerDileptonSF", &TriggerDileptonSF, &b_TriggerDileptonSF);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("LB", &LB, &b_LB);
   fChain->SetBranchAddress("passGRL", &passGRL, &b_passGRL);
   fChain->SetBranchAddress("RunNb", &RunNb, &b_RunNb);
   fChain->SetBranchAddress("PRWrandomRunNumber", &PRWrandomRunNumber, &b_PRWrandomRunNumber);
   fChain->SetBranchAddress("DetError", &DetError, &b_DetError);
   fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("NMu", &NMu, &b_NMu);
   fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
   fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
   fChain->SetBranchAddress("Mu_pT", &Mu_pT, &b_Mu_pT);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("Mu_sigd0", &Mu_sigd0, &b_Mu_sigd0);
   fChain->SetBranchAddress("Mu_d0pvtx", &Mu_d0pvtx, &b_Mu_d0pvtx);
   fChain->SetBranchAddress("Mu_d0pvtxerr", &Mu_d0pvtxerr, &b_Mu_d0pvtxerr);
   fChain->SetBranchAddress("Mu_z0sinTheta", &Mu_z0sinTheta, &b_Mu_z0sinTheta);
   fChain->SetBranchAddress("Mu_isBad", &Mu_isBad, &b_Mu_isBad);
   fChain->SetBranchAddress("Mu_isCosmic", &Mu_isCosmic, &b_Mu_isCosmic);
   fChain->SetBranchAddress("Mu_passOR", &Mu_passOR, &b_Mu_passOR);
   fChain->SetBranchAddress("Mu_isTight", &Mu_isTight, &b_Mu_isTight);
   fChain->SetBranchAddress("Mu_isSig2017", &Mu_isSig2017, &b_Mu_isSig2017);
   fChain->SetBranchAddress("Mu_PromptLepTaggerBDTweight", &Mu_PromptLepTaggerBDTweight, &b_Mu_PromptLepTaggerBDTweight);
   fChain->SetBranchAddress("Mu_type", &Mu_type, &b_Mu_type);
   fChain->SetBranchAddress("Mu_origin", &Mu_origin, &b_Mu_origin);
   fChain->SetBranchAddress("Mu_ptcone20", &Mu_ptcone20, &b_Mu_ptcone20);
   fChain->SetBranchAddress("Mu_ptcone20_TightTTVA_pt500", &Mu_ptcone20_TightTTVA_pt500, &b_Mu_ptcone20_TightTTVA_pt500);
   fChain->SetBranchAddress("Mu_ptcone20_TightTTVA_pt1000", &Mu_ptcone20_TightTTVA_pt1000, &b_Mu_ptcone20_TightTTVA_pt1000);
   fChain->SetBranchAddress("Mu_ptcone30", &Mu_ptcone30, &b_Mu_ptcone30);
   fChain->SetBranchAddress("Mu_ptcone40", &Mu_ptcone40, &b_Mu_ptcone40);
   fChain->SetBranchAddress("Mu_ptvarcone20", &Mu_ptvarcone20, &b_Mu_ptvarcone20);
   fChain->SetBranchAddress("Mu_ptvarcone30", &Mu_ptvarcone30, &b_Mu_ptvarcone30);
   fChain->SetBranchAddress("Mu_ptvarcone30_TightTTVA_pt500", &Mu_ptvarcone30_TightTTVA_pt500, &b_Mu_ptvarcone30_TightTTVA_pt500);
   fChain->SetBranchAddress("Mu_ptvarcone30_TightTTVA_pt1000", &Mu_ptvarcone30_TightTTVA_pt1000, &b_Mu_ptvarcone30_TightTTVA_pt1000);
   fChain->SetBranchAddress("Mu_ptvarcone40", &Mu_ptvarcone40, &b_Mu_ptvarcone40);
   fChain->SetBranchAddress("Mu_topoetcone20", &Mu_topoetcone20, &b_Mu_topoetcone20);
   fChain->SetBranchAddress("Mu_topoetcone30", &Mu_topoetcone30, &b_Mu_topoetcone30);
   fChain->SetBranchAddress("Mu_topoetcone40", &Mu_topoetcone40, &b_Mu_topoetcone40);
   fChain->SetBranchAddress("Mu_neflowisol20", &Mu_neflowisol20, &b_Mu_neflowisol20);
   fChain->SetBranchAddress("Mu_passIsoLooseTO", &Mu_passIsoLooseTO, &b_Mu_passIsoLooseTO);
   fChain->SetBranchAddress("Mu_passIsoLoose", &Mu_passIsoLoose, &b_Mu_passIsoLoose);
   fChain->SetBranchAddress("Mu_passIsoTight", &Mu_passIsoTight, &b_Mu_passIsoTight);
   fChain->SetBranchAddress("Mu_passIsoGrad", &Mu_passIsoGrad, &b_Mu_passIsoGrad);
   fChain->SetBranchAddress("Mu_passIsoGradCustomTight", &Mu_passIsoGradCustomTight, &b_Mu_passIsoGradCustomTight);
   fChain->SetBranchAddress("Mu_passIsoGradCustom", &Mu_passIsoGradCustom, &b_Mu_passIsoGradCustom);
   fChain->SetBranchAddress("Mu_passIsoGradLoose", &Mu_passIsoGradLoose, &b_Mu_passIsoGradLoose);
   fChain->SetBranchAddress("Mu_SFw", &Mu_SFw, &b_Mu_SFw);
   fChain->SetBranchAddress("Mu_IsoSFw", &Mu_IsoSFw, &b_Mu_IsoSFw);
   fChain->SetBranchAddress("Mu_StatUncReco", &Mu_StatUncReco, &b_Mu_StatUncReco);
   fChain->SetBranchAddress("Mu_SystUncReco", &Mu_SystUncReco, &b_Mu_SystUncReco);
   fChain->SetBranchAddress("Mu_StatUncReco_LOWPT", &Mu_StatUncReco_LOWPT, &b_Mu_StatUncReco_LOWPT);
   fChain->SetBranchAddress("Mu_SystUncReco_LOWPT", &Mu_SystUncReco_LOWPT, &b_Mu_SystUncReco_LOWPT);
   fChain->SetBranchAddress("Mu_StatUncISO", &Mu_StatUncISO, &b_Mu_StatUncISO);
   fChain->SetBranchAddress("Mu_SystUncISO", &Mu_SystUncISO, &b_Mu_SystUncISO);
   fChain->SetBranchAddress("Mu_trigMatch_mu26_ivarmedium", &Mu_trigMatch_mu26_ivarmedium, &b_Mu_trigMatch_mu26_ivarmedium);
   fChain->SetBranchAddress("Mu_trigMatch_mu20_iloose_L1MU15", &Mu_trigMatch_mu20_iloose_L1MU15, &b_Mu_trigMatch_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("Mu_trigMatch_mu40", &Mu_trigMatch_mu40, &b_Mu_trigMatch_mu40);
   fChain->SetBranchAddress("Mu_trigMatch_mu50", &Mu_trigMatch_mu50, &b_Mu_trigMatch_mu50);
   fChain->SetBranchAddress("Mu_trigMatch_mu8noL1", &Mu_trigMatch_mu8noL1, &b_Mu_trigMatch_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatch_mu14", &Mu_trigMatch_mu14, &b_Mu_trigMatch_mu14);
   fChain->SetBranchAddress("Mu_trigMatch_mu18", &Mu_trigMatch_mu18, &b_Mu_trigMatch_mu18);
   fChain->SetBranchAddress("Mu_trigMatch_mu18_mu8noL1", &Mu_trigMatch_mu18_mu8noL1, &b_Mu_trigMatch_mu18_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatch_e17_lhloose_mu14", &Mu_trigMatch_e17_lhloose_mu14, &b_Mu_trigMatch_e17_lhloose_mu14);
   fChain->SetBranchAddress("Mu_trigMatch_e17_lhloose_nod0_mu14", &Mu_trigMatch_e17_lhloose_nod0_mu14, &b_Mu_trigMatch_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("Mu_trigMatch_mu20_mu8noL1", &Mu_trigMatch_mu20_mu8noL1, &b_Mu_trigMatch_mu20_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatch_mu22_mu8noL1", &Mu_trigMatch_mu22_mu8noL1, &b_Mu_trigMatch_mu22_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatch_mu24_iloose", &Mu_trigMatch_mu24_iloose, &b_Mu_trigMatch_mu24_iloose);
   fChain->SetBranchAddress("Mu_trigMatch_mu24_ivarloose", &Mu_trigMatch_mu24_ivarloose, &b_Mu_trigMatch_mu24_ivarloose);
   fChain->SetBranchAddress("Mu_trigMatch_mu24_iloose_L1MU15", &Mu_trigMatch_mu24_iloose_L1MU15, &b_Mu_trigMatch_mu24_iloose_L1MU15);
   fChain->SetBranchAddress("Mu_trigMatch_mu24_ivarloose_L1MU15", &Mu_trigMatch_mu24_ivarloose_L1MU15, &b_Mu_trigMatch_mu24_ivarloose_L1MU15);
   fChain->SetBranchAddress("Mu_trigMatchPair_mu18_mu8noL1", &Mu_trigMatchPair_mu18_mu8noL1, &b_Mu_trigMatchPair_mu18_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatchPair_mu20_mu8noL1", &Mu_trigMatchPair_mu20_mu8noL1, &b_Mu_trigMatchPair_mu20_mu8noL1);
   fChain->SetBranchAddress("Mu_trigMatchPair_mu22_mu8noL1", &Mu_trigMatchPair_mu22_mu8noL1, &b_Mu_trigMatchPair_mu22_mu8noL1);
   fChain->SetBranchAddress("NEl", &NEl, &b_NEl);
   fChain->SetBranchAddress("El_eta", &El_eta, &b_El_eta);
   fChain->SetBranchAddress("El_etaclus", &El_etaclus, &b_El_etaclus);
   fChain->SetBranchAddress("El_phi", &El_phi, &b_El_phi);
   fChain->SetBranchAddress("El_pT", &El_pT, &b_El_pT);
   fChain->SetBranchAddress("El_E", &El_E, &b_El_E);
   fChain->SetBranchAddress("El_charge", &El_charge, &b_El_charge);
   fChain->SetBranchAddress("El_sigd0", &El_sigd0, &b_El_sigd0);
   fChain->SetBranchAddress("El_d0pvtx", &El_d0pvtx, &b_El_d0pvtx);
   fChain->SetBranchAddress("El_d0pvtxerr", &El_d0pvtxerr, &b_El_d0pvtxerr);
   fChain->SetBranchAddress("El_z0sinTheta", &El_z0sinTheta, &b_El_z0sinTheta);
   fChain->SetBranchAddress("El_isLooseAndBLayerLH_baseline", &El_isLooseAndBLayerLH_baseline, &b_El_isLooseAndBLayerLH_baseline);
   fChain->SetBranchAddress("El_isLooseAndBLayerLH_fromTool", &El_isLooseAndBLayerLH_fromTool, &b_El_isLooseAndBLayerLH_fromTool);
   fChain->SetBranchAddress("El_isMediumLH", &El_isMediumLH, &b_El_isMediumLH);
   fChain->SetBranchAddress("El_isTightLH", &El_isTightLH, &b_El_isTightLH);
   fChain->SetBranchAddress("El_isSigNoCFT2017", &El_isSigNoCFT2017, &b_El_isSigNoCFT2017);
   fChain->SetBranchAddress("El_nBLayerHits", &El_nBLayerHits, &b_El_nBLayerHits);
   fChain->SetBranchAddress("El_expectBLayerHit", &El_expectBLayerHit, &b_El_expectBLayerHit);
   fChain->SetBranchAddress("El_passOR", &El_passOR, &b_El_passOR);
   fChain->SetBranchAddress("El_passChargeFlipTaggerBDTmedium", &El_passChargeFlipTaggerBDTmedium, &b_El_passChargeFlipTaggerBDTmedium);
   fChain->SetBranchAddress("El_passChargeFlipTaggerBDTloose", &El_passChargeFlipTaggerBDTloose, &b_El_passChargeFlipTaggerBDTloose);
   fChain->SetBranchAddress("El_PromptLepTaggerBDTweight", &El_PromptLepTaggerBDTweight, &b_El_PromptLepTaggerBDTweight);
   fChain->SetBranchAddress("El_truthType", &El_truthType, &b_El_truthType);
   fChain->SetBranchAddress("El_truthOrigin", &El_truthOrigin, &b_El_truthOrigin);
   fChain->SetBranchAddress("El_truthPdgId", &El_truthPdgId, &b_El_truthPdgId);
   fChain->SetBranchAddress("El_bkgTruthType", &El_bkgTruthType, &b_El_bkgTruthType);
   fChain->SetBranchAddress("El_bkgTruthOrigin", &El_bkgTruthOrigin, &b_El_bkgTruthOrigin);
   fChain->SetBranchAddress("El_bkgMotherPdgId", &El_bkgMotherPdgId, &b_El_bkgMotherPdgId);
   fChain->SetBranchAddress("El_firstEgMotherTruthType", &El_firstEgMotherTruthType, &b_El_firstEgMotherTruthType);
   fChain->SetBranchAddress("El_firstEgMotherTruthOrigin", &El_firstEgMotherTruthOrigin, &b_El_firstEgMotherTruthOrigin);
   fChain->SetBranchAddress("El_firstEgMotherPdgId", &El_firstEgMotherPdgId, &b_El_firstEgMotherPdgId);
   fChain->SetBranchAddress("El_lastEgMotherTruthType", &El_lastEgMotherTruthType, &b_El_lastEgMotherTruthType);
   fChain->SetBranchAddress("El_lastEgMotherTruthOrigin", &El_lastEgMotherTruthOrigin, &b_El_lastEgMotherTruthOrigin);
   fChain->SetBranchAddress("El_lastEgMotherPdgId", &El_lastEgMotherPdgId, &b_El_lastEgMotherPdgId);
   fChain->SetBranchAddress("El_chFlip", &El_chFlip, &b_El_chFlip);
   fChain->SetBranchAddress("El_ptcone20", &El_ptcone20, &b_El_ptcone20);
   fChain->SetBranchAddress("El_ptcone20_TightTTVA_pt500", &El_ptcone20_TightTTVA_pt500, &b_El_ptcone20_TightTTVA_pt500);
   fChain->SetBranchAddress("El_ptvarcone20_TightTTVA_pt1000", &El_ptvarcone20_TightTTVA_pt1000, &b_El_ptvarcone20_TightTTVA_pt1000);
   fChain->SetBranchAddress("El_ptcone30", &El_ptcone30, &b_El_ptcone30);
   fChain->SetBranchAddress("El_ptcone40", &El_ptcone40, &b_El_ptcone40);
   fChain->SetBranchAddress("El_ptvarcone20", &El_ptvarcone20, &b_El_ptvarcone20);
//    fChain->SetBranchAddress("El_ptvarcone20_TightTTVA_pt1000", &El_ptvarcone20_TightTTVA_pt1000, &b_El_ptvarcone20_TightTTVA_pt1000);
   fChain->SetBranchAddress("El_ptvarcone30", &El_ptvarcone30, &b_El_ptvarcone30);
   fChain->SetBranchAddress("El_ptvarcone30_TightTTVA_pt1000", &El_ptvarcone30_TightTTVA_pt1000, &b_El_ptvarcone30_TightTTVA_pt1000);
   fChain->SetBranchAddress("El_ptvarcone40", &El_ptvarcone40, &b_El_ptvarcone40);
   fChain->SetBranchAddress("El_topoetcone20", &El_topoetcone20, &b_El_topoetcone20);
   fChain->SetBranchAddress("El_topoetcone30", &El_topoetcone30, &b_El_topoetcone30);
   fChain->SetBranchAddress("El_topoetcone40", &El_topoetcone40, &b_El_topoetcone40);
   fChain->SetBranchAddress("El_neflowisol20", &El_neflowisol20, &b_El_neflowisol20);
   fChain->SetBranchAddress("El_passIsoLooseTO", &El_passIsoLooseTO, &b_El_passIsoLooseTO);
   fChain->SetBranchAddress("El_passIsoLoose", &El_passIsoLoose, &b_El_passIsoLoose);
   fChain->SetBranchAddress("El_passIsoTight", &El_passIsoTight, &b_El_passIsoTight);
   fChain->SetBranchAddress("El_passIsoGrad", &El_passIsoGrad, &b_El_passIsoGrad);
   fChain->SetBranchAddress("El_passIsoGradCustomTight", &El_passIsoGradCustomTight, &b_El_passIsoGradCustomTight);
   fChain->SetBranchAddress("El_passIsoGradCustom", &El_passIsoGradCustom, &b_El_passIsoGradCustom);
   fChain->SetBranchAddress("El_passIsoGradLoose", &El_passIsoGradLoose, &b_El_passIsoGradLoose);
   fChain->SetBranchAddress("El_SFwUncReco", &El_SFwUncReco, &b_El_SFwUncReco);
   fChain->SetBranchAddress("El_SFwTightLH", &El_SFwTightLH, &b_El_SFwTightLH);
   fChain->SetBranchAddress("El_SFwMediumLH", &El_SFwMediumLH, &b_El_SFwMediumLH);
   fChain->SetBranchAddress("El_SFwUncMediumLH", &El_SFwUncMediumLH, &b_El_SFwUncMediumLH);
   fChain->SetBranchAddress("El_SFwLooseAndBLayerLH", &El_SFwLooseAndBLayerLH, &b_El_SFwLooseAndBLayerLH);
   fChain->SetBranchAddress("El_SFweightCFT", &El_SFweightCFT, &b_El_SFweightCFT);
   fChain->SetBranchAddress("El_SFUncweightCFT", &El_SFUncweightCFT, &b_El_SFUncweightCFT);
   fChain->SetBranchAddress("El_SFweightCFID", &El_SFweightCFID, &b_El_SFweightCFID);
   fChain->SetBranchAddress("El_SFStatweightCFID", &El_SFStatweightCFID, &b_El_SFStatweightCFID);
   fChain->SetBranchAddress("El_SFSystweightCFID", &El_SFSystweightCFID, &b_El_SFSystweightCFID);
   fChain->SetBranchAddress("El_IsoSFwMediumLH", &El_IsoSFwMediumLH, &b_El_IsoSFwMediumLH);
   fChain->SetBranchAddress("El_IsoSFwUncMediumLH", &El_IsoSFwUncMediumLH, &b_El_IsoSFwUncMediumLH);
   fChain->SetBranchAddress("El_trigMatch_e12_lhloose_L1EM10VH", &El_trigMatch_e12_lhloose_L1EM10VH, &b_El_trigMatch_e12_lhloose_L1EM10VH);
   fChain->SetBranchAddress("El_trigMatch_e17_lhloose", &El_trigMatch_e17_lhloose, &b_El_trigMatch_e17_lhloose);
   fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_L1EM20VH", &El_trigMatch_e24_lhmedium_L1EM20VH, &b_El_trigMatch_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_iloose_L1EM20VH", &El_trigMatch_e24_lhmedium_iloose_L1EM20VH, &b_El_trigMatch_e24_lhmedium_iloose_L1EM20VH);
   fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_nod0_ivarloose", &El_trigMatch_e24_lhmedium_nod0_ivarloose, &b_El_trigMatch_e24_lhmedium_nod0_ivarloose);
   fChain->SetBranchAddress("El_trigMatch_e24_lhtight_nod0_ivarloose", &El_trigMatch_e24_lhtight_nod0_ivarloose, &b_El_trigMatch_e24_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("El_trigMatch_e26_lhtight_nod0_ivarloose", &El_trigMatch_e26_lhtight_nod0_ivarloose, &b_El_trigMatch_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("El_trigMatch_e60_lhmedium", &El_trigMatch_e60_lhmedium, &b_El_trigMatch_e60_lhmedium);
   fChain->SetBranchAddress("El_trigMatch_e60_lhmedium_nod0", &El_trigMatch_e60_lhmedium_nod0, &b_El_trigMatch_e60_lhmedium_nod0);
   fChain->SetBranchAddress("El_trigMatch_2e12_lhloose_L12EM10VH", &El_trigMatch_2e12_lhloose_L12EM10VH, &b_El_trigMatch_2e12_lhloose_L12EM10VH);
   fChain->SetBranchAddress("El_trigMatch_2e15_lhloose_L12EM10VH", &El_trigMatch_2e15_lhloose_L12EM10VH, &b_El_trigMatch_2e15_lhloose_L12EM10VH);
   fChain->SetBranchAddress("El_trigMatch_2e15_lhvloose_L12EM13VH", &El_trigMatch_2e15_lhvloose_L12EM13VH, &b_El_trigMatch_2e15_lhvloose_L12EM13VH);
   fChain->SetBranchAddress("El_trigMatch_2e15_lhvloose_nod0_L12EM13VH", &El_trigMatch_2e15_lhvloose_nod0_L12EM13VH, &b_El_trigMatch_2e15_lhvloose_nod0_L12EM13VH);
   fChain->SetBranchAddress("El_trigMatch_2e17_lhvloose_nod0", &El_trigMatch_2e17_lhvloose_nod0, &b_El_trigMatch_2e17_lhvloose_nod0);
   fChain->SetBranchAddress("El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI", &El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI, &b_El_TrigMatch_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("El_TrigMatch_2e24_lhvloose_nod0", &El_TrigMatch_2e24_lhvloose_nod0, &b_El_TrigMatch_2e24_lhvloose_nod0);
   fChain->SetBranchAddress("El_trigMatch_e17_lhloose_mu14", &El_trigMatch_e17_lhloose_mu14, &b_El_trigMatch_e17_lhloose_mu14);
   fChain->SetBranchAddress("El_trigMatch_e17_lhloose_nod0_mu14", &El_trigMatch_e17_lhloose_nod0_mu14, &b_El_trigMatch_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("NJet", &NJet, &b_NJet);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pT", &Jet_pT, &b_Jet_pT);
   fChain->SetBranchAddress("Jet_E", &Jet_E, &b_Jet_E);
   fChain->SetBranchAddress("Jet_nTrk", &Jet_nTrk, &b_Jet_nTrk);
   fChain->SetBranchAddress("Jet_quality", &Jet_quality, &b_Jet_quality);
   fChain->SetBranchAddress("Jet_EMFrac", &Jet_EMFrac, &b_Jet_EMFrac);
   fChain->SetBranchAddress("Jet_HECFrac", &Jet_HECFrac, &b_Jet_HECFrac);
   fChain->SetBranchAddress("Jet_JVT", &Jet_JVT, &b_Jet_JVT);
   fChain->SetBranchAddress("Jet_JVTsf", &Jet_JVTsf, &b_Jet_JVTsf);
   fChain->SetBranchAddress("totalJVTsf", &totalJVTsf, &b_totalJVTsf);
   fChain->SetBranchAddress("Jet_MV2c10", &Jet_MV2c10, &b_Jet_MV2c10);
   fChain->SetBranchAddress("Jet_SFw", &Jet_SFw, &b_Jet_SFw);
   fChain->SetBranchAddress("Jet_passOR", &Jet_passOR, &b_Jet_passOR);
   fChain->SetBranchAddress("Jet_ConeTruthLabel", &Jet_ConeTruthLabel, &b_Jet_ConeTruthLabel);
   fChain->SetBranchAddress("Jet_PartonTruthLabel", &Jet_PartonTruthLabel, &b_Jet_PartonTruthLabel);
   fChain->SetBranchAddress("Jet_HadronConeExclTruthLabel", &Jet_HadronConeExclTruthLabel, &b_Jet_HadronConeExclTruthLabel);
   fChain->SetBranchAddress("Jet_deltaR", &Jet_deltaR, &b_Jet_deltaR);
   fChain->SetBranchAddress("Etmiss_TST_Etx", &Etmiss_TST_Etx, &b_Etmiss_TST_Etx);
   fChain->SetBranchAddress("Etmiss_TST_Ety", &Etmiss_TST_Ety, &b_Etmiss_TST_Ety);
   fChain->SetBranchAddress("Etmiss_TST_Et", &Etmiss_TST_Et, &b_Etmiss_TST_Et);
   fChain->SetBranchAddress("Etmiss_Truth_Etx", &Etmiss_Truth_Etx, &b_Etmiss_Truth_Etx);
   fChain->SetBranchAddress("Etmiss_Truth_Ety", &Etmiss_Truth_Ety, &b_Etmiss_Truth_Ety);
   fChain->SetBranchAddress("Etmiss_Truth_Et", &Etmiss_Truth_Et, &b_Etmiss_Truth_Et);
   fChain->SetBranchAddress("NTruthAntiktJet", &NTruthAntiktJet, &b_NTruthAntiktJet);
   fChain->SetBranchAddress("TruthAntiktJet_eta", &TruthAntiktJet_eta, &b_TruthAntiktJet_eta);
   fChain->SetBranchAddress("TruthAntiktJet_phi", &TruthAntiktJet_phi, &b_TruthAntiktJet_phi);
   fChain->SetBranchAddress("TruthAntiktJet_pT", &TruthAntiktJet_pT, &b_TruthAntiktJet_pT);
   fChain->SetBranchAddress("TruthAntiktJet_E", &TruthAntiktJet_E, &b_TruthAntiktJet_E);
   fChain->SetBranchAddress("TruthAntiktJet_ConeTruthLabel", &TruthAntiktJet_ConeTruthLabel, &b_TruthAntiktJet_ConeTruthLabel);
   fChain->SetBranchAddress("TruthAntiktJet_HadronConeExclTruthLabelID", &TruthAntiktJet_HadronConeExclTruthLabelID, &b_TruthAntiktJet_HadronConeExclTruthLabelID);
   fChain->SetBranchAddress("TruthAntiktJet_PartonTruthLabel", &TruthAntiktJet_PartonTruthLabel, &b_TruthAntiktJet_PartonTruthLabel);
   fChain->SetBranchAddress("TruthAntiktJetJet_ClassHF", &TruthAntiktJetJet_ClassHF, &b_TruthAntiktJetJet_ClassHF);
   fChain->SetBranchAddress("extraB", &extraB, &b_extraB);
   fChain->SetBranchAddress("extraC", &extraC, &b_extraC);
   fChain->SetBranchAddress("NTruthJet", &NTruthJet, &b_NTruthJet);
   fChain->SetBranchAddress("TruthJet_eta", &TruthJet_eta, &b_TruthJet_eta);
   fChain->SetBranchAddress("TruthJet_phi", &TruthJet_phi, &b_TruthJet_phi);
   fChain->SetBranchAddress("TruthJet_pT", &TruthJet_pT, &b_TruthJet_pT);
   fChain->SetBranchAddress("TruthJet_E", &TruthJet_E, &b_TruthJet_E);
   fChain->SetBranchAddress("TruthJet_id", &TruthJet_id, &b_TruthJet_id);
   fChain->SetBranchAddress("TruthJet_origin", &TruthJet_origin, &b_TruthJet_origin);
   fChain->SetBranchAddress("NTruthRealL", &NTruthRealL, &b_NTruthRealL);
   fChain->SetBranchAddress("TruthRealL_eta", &TruthRealL_eta, &b_TruthRealL_eta);
   fChain->SetBranchAddress("TruthRealL_phi", &TruthRealL_phi, &b_TruthRealL_phi);
   fChain->SetBranchAddress("TruthRealL_pT", &TruthRealL_pT, &b_TruthRealL_pT);
   fChain->SetBranchAddress("TruthRealL_id", &TruthRealL_id, &b_TruthRealL_id);
   fChain->SetBranchAddress("TruthRealL_origin", &TruthRealL_origin, &b_TruthRealL_origin);
   fChain->SetBranchAddress("SUSY_Spart_pdgId1", &SUSY_Spart_pdgId1, &b_SUSY_Spart_pdgId1);
   fChain->SetBranchAddress("SUSY_Spart_pdgId2", &SUSY_Spart_pdgId2, &b_SUSY_Spart_pdgId2);
   fChain->SetBranchAddress("SUSY_Gluino_decay1", &SUSY_Gluino_decay1, &b_SUSY_Gluino_decay1);
   fChain->SetBranchAddress("SUSY_Gluino_decay2", &SUSY_Gluino_decay2, &b_SUSY_Gluino_decay2);
   fChain->SetBranchAddress("GenFiltHT", &GenFiltHT, &b_GenFiltHT);
   fChain->SetBranchAddress("GenFiltMET", &GenFiltMET, &b_GenFiltMET);
   fChain->SetBranchAddress("TruthX1", &TruthX1, &b_TruthX1);
   fChain->SetBranchAddress("TruthX2", &TruthX2, &b_TruthX2);
   fChain->SetBranchAddress("TruthQ", &TruthQ, &b_TruthQ);
   fChain->SetBranchAddress("TruthPDGID1", &TruthPDGID1, &b_TruthPDGID1);
   fChain->SetBranchAddress("TruthPDGID2", &TruthPDGID2, &b_TruthPDGID2);
   fChain->SetBranchAddress("SherpaNjetWeight", &SherpaNjetWeight, &b_SherpaNjetWeight);
   Notify();
}

Bool_t bTagCutFlow::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bTagCutFlow::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bTagCutFlow::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bTagCutFlow_cxx
