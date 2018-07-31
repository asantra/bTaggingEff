#define bTagCutFlow_cxx
#include "bTagCutFlow.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;
void bTagCutFlow::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L bTagCutFlow.C
//      Root > bTagCutFlow t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   
   Long64_t beforeCut(0), primaryVertex(0), grl(0), trigger(0);
   Long64_t badMuonNo(0), jetPass(0), badJetNo(0), signalJetNo(0);
   Long64_t cosmicMuonNo(0), baselineLepton(0), signalLepton(0);
   Long64_t channelEE(0), triggerMatchEE(0), bJetEE(0), fourJetEE(0);
   Long64_t metEE(0);
   Long64_t channelEM(0), triggerMatchEM(0), bJetEM(0), fourJetEM(0);
   Long64_t metEM(0);
   Long64_t channelMM(0), triggerMatchMM(0), bJetMM(0), fourJetMM(0);
   Long64_t metMM(0);
   
   
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      beforeCut++;
      if(jentry%10000==0)std::cout << "Events processed: " << jentry << std::endl;
      // pass GRL
      if(passGRL){
        grl++;
        // pass primary vertex
        if(PV_z > -999){
          primaryVertex++;
          // pass trigger
          bool triggerOutput = false;
          if(PRWrandomRunNumber < 290000){
            triggerOutput = false;
          }
          else if(PRWrandomRunNumber < 320000){
            if(Etmiss_TST_Et < 250){
              triggerOutput = (HLT_2e17_lhvloose_nod0 || HLT_e17_lhloose_nod0_mu14 || HLT_mu22_mu8noL1 || HLT_xe100_mht_L1XE50 || HLT_xe110_mht_L1XE50); 
            }
            else{
              triggerOutput = (HLT_2e17_lhvloose_nod0 || HLT_e17_lhloose_nod0_mu14 || HLT_mu22_mu8noL1 || HLT_xe100_mht_L1XE50 || HLT_xe110_mht_L1XE50 || HLT_xe90_mht_L1XE50);
            }
          }
          else{
            if(Etmiss_TST_Et < 250){
               triggerOutput =  (HLT_2e17_lhvloose_nod0_L12EM15VHI || HLT_2e24_lhvloose_nod0 || HLT_e17_lhloose_nod0_mu14 || HLT_xe110_pufit_L1XE55);
            }
            else{
                triggerOutput =  (HLT_2e17_lhvloose_nod0_L12EM15VHI || HLT_2e24_lhvloose_nod0 || HLT_e17_lhloose_nod0_mu14 || HLT_xe110_pufit_L1XE55 || HLT_mu22_mu8noL1);
            }
                //, HLT_2e24_lhvloose_nod0 (use 2e17_lhloose||2e24_lhvloose for RUNS when 2e17_lhvloose is unprescaled, use 2e24_lhvloose for runs when 2e17_lhvloose is prescaled), , , (only for MET>250 GeV) 
          }
          if(triggerOutput){
            trigger++;
            // pass bad muon
            bool badMuon(false);
            for (size_t i=0; i < Mu_isBad->size(); ++i){
                if(Mu_isBad->at(i)){
                    badMuon = true;
                    break;
                }
            }
            if(!badMuon){
                badMuonNo++;
                // pass jetPassOR
                bool jetPassOR = false;
                for (size_t i=0; i < Jet_passOR->size(); ++i){
                    if(Jet_passOR->at(i)){
                        jetPassOR=true;
                        break;
                    }
                }
                if(jetPassOR){
                    jetPass++;
                    /// jet pass quality
                    bool jetPassQuality = false;
                    for (size_t i=0; i < Jet_quality->size(); ++i){
                        if(Jet_quality->at(i)){
                            jetPassQuality=true;
                            break;
                        }
                    }
                    if(!jetPassQuality){
                        // bad jet
                        badJetNo++;
                        int goodJetNo(0);
                        for(size_t i=0; i < Jet_eta->size(); ++i){
                            if(Jet_pT->at(i) > 20 && fabs(Jet_eta->at(i)) < 2.8){
                                if(!(Jet_pT->at(i) < 60 && fabs(Jet_eta->at(i)) < 2.4 && Jet_JVT->at(i) < 0.59))goodJetNo++;
                            }
                        }
                        if(goodJetNo>0){
                            // signal jet
                            signalJetNo++;
                            // cosmic muon
                            bool cosmicMuon = false;
                            for(size_t i = 0; i < Mu_isCosmic->size(); ++i){
                                if(Mu_isCosmic->at(i)){
                                    cosmicMuon = true;
                                    break;
                                }
                            }
                            if(!cosmicMuon){
                                cosmicMuonNo++;
                                // baselinelepton
                                bool baselineLeptonBool(false);
                                if(NMu > 1 && Mu_pT->at(0) > 10 && Mu_pT->at(1) > 10 && fabs(Mu_eta->at(0)) < 2.5 && fabs(Mu_eta->at(1)) < 2.5 && Mu_z0sinTheta->at(0) < 0.5 && Mu_z0sinTheta->at(1) < 0.5){
                                    baselineLeptonBool = true;
                                }
                                else if(NEl > 1 && El_pT->at(0) > 10 && El_pT->at(1) > 10 && fabs(El_eta->at(0)) < 2.47 && fabs(El_eta->at(1)) < 2.47 && El_z0sinTheta->at(0) < 0.5 && El_z0sinTheta->at(1) < 0.5 && El_sigd0->at(0) < 5 && El_sigd0->at(1) < 5){
                                    if(!((fabs(El_eta->at(0)) > 1.37 && fabs(El_eta->at(0)) < 1.52) || (fabs(El_eta->at(1)) > 1.37 && fabs(El_eta->at(1)) < 1.52)))baselineLeptonBool = true;
                                }
                                else if(NEl > 0 && El_pT->at(0) > 10 && fabs(El_eta->at(0)) < 2.47 && El_z0sinTheta->at(0) < 0.5 && El_sigd0->at(0) < 5 && fabs(El_eta->at(0)) > 1.52 && fabs(El_eta->at(0)) < 1.37 
                                    && NMu > 0 && Mu_pT->at(0) > 10 && fabs(Mu_eta->at(0)) < 2.5 && Mu_z0sinTheta->at(0) < 0.5){
                                    baselineLeptonBool = true;
                                }
                                if(baselineLeptonBool){
                                    baselineLepton++;
                                    // signal lepton
                                    bool signalLeptonBool(false);
                                    if(NMu > 1 && Mu_ptvarcone30->at(0)/Mu_pT->at(0) < 0.06 && Mu_ptvarcone30->at(1)/Mu_pT->at(1) < 0.06 && Mu_sigd0->at(0) < 3 && Mu_sigd0->at(1) < 3 && Mu_pT->at(0)>20 && Mu_pT->at(1)>20){
                                      signalLeptonBool = true;
                                    }
                                    else if(NEl > 1 && El_pT->at(0) > 20 && El_pT->at(1) > 20 && fabs(El_eta->at(0)) < 2.0 && fabs(El_eta->at(1)) < 2.0 && El_ptvarcone20->at(0)/El_pT->at(0) < 0.06 && El_ptvarcone20->at(1)/El_pT->at(1) < 0.06 && El_topoetcone20->at(0)/El_pT->at(0) < 0.06 && El_topoetcone20->at(1)/El_pT->at(1) < 0.06 && El_chFlip->at(0)>-0.28087 && El_chFlip->at(1)>-0.28087){
                                      signalLeptonBool = true;
                                    }
                                    else if(NEl > 0 && El_pT->at(0) > 20 && fabs(El_eta->at(0)) < 2.0 && El_ptvarcone20->at(0)/El_pT->at(0) < 0.06 && El_topoetcone20->at(0)/El_pT->at(0) < 0.06 && El_chFlip->at(0)>-0.28087 
                                        && NMu > 0 && Mu_ptvarcone30->at(0)/Mu_pT->at(0) < 0.06 && Mu_sigd0->at(0) < 3 && Mu_pT->at(0)>20 ){
                                        signalLeptonBool = true;
                                    }
                                    if(signalLeptonBool){
                                        signalLepton++;
                                        
                                        // channel separation
                                        // MM
                                        if(ChannelNumber==0){
                                            if(Mu_charge->at(0) == Mu_charge->at(1)){
                                                channelMM++;
                                                int numberBJet(0);
                                                for(size_t i=0; i < Jet_pT->size(); ++i){
                                                    if(Jet_pT->at(i) > 20 && fabs(Jet_eta->at(i)) < 2.5 && !(Jet_pT->at(i) < 60 && fabs(Jet_eta->at(i)) < 2.4 && Jet_JVT->at(i) < 0.59) && Jet_MV2c10->at(i) > 0.83)numberBJet++;
                                                }
                                                if(numberBJet > 0){
                                                    bJetMM++;
                                                    if(NJet > 3 && Jet_pT->at(0) > 50 && Jet_pT->at(1) > 50 && Jet_pT->at(2) > 50 && Jet_pT->at(3) > 50){
                                                        fourJetMM++;
                                                        if(Etmiss_TST_Et > 125){
                                                            metMM++;
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                        // EE
                                        else if(ChannelNumber==1){
                                            if(El_charge->at(0) == El_charge->at(1)){
                                                channelEE++;
                                                int numberBJet(0);
                                                for(size_t i=0; i < Jet_pT->size(); ++i){
                                                    if(Jet_pT->at(i) > 20 && fabs(Jet_eta->at(i)) < 2.5 && !(Jet_pT->at(i) < 60 && fabs(Jet_eta->at(i)) < 2.4 && Jet_JVT->at(i) < 0.59) && Jet_MV2c10->at(i) > 0.83)numberBJet++;
                                                }
                                                if(numberBJet > 0){
                                                    bJetEE++;
                                                    if(NJet > 3 && Jet_pT->at(0) > 50 && Jet_pT->at(1) > 50 && Jet_pT->at(2) > 50 && Jet_pT->at(3) > 50){
                                                        fourJetEE++;
                                                        if(Etmiss_TST_Et > 125){
                                                            metEE++;
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                        // EM or ME
                                        else if(ChannelNumber==2 || ChannelNumber==3){
                                            if(Mu_charge->at(0) == El_charge->at(0)){
                                                channelEM++;
                                                int numberBJet(0);
                                                for(size_t i=0; i < Jet_pT->size(); ++i){
                                                    if(Jet_pT->at(i) > 20 && fabs(Jet_eta->at(i)) < 2.5 && !(Jet_pT->at(i) < 60 && fabs(Jet_eta->at(i)) < 2.4 && Jet_JVT->at(i) < 0.59) && Jet_MV2c10->at(i) > 0.83)numberBJet++;
                                                }
                                                if(numberBJet > 0){
                                                    bJetEM++;
                                                    if(NJet > 3 && Jet_pT->at(0) > 50 && Jet_pT->at(1) > 50 && Jet_pT->at(2) > 50 && Jet_pT->at(3) > 50){
                                                        fourJetEM++;
                                                        if(Etmiss_TST_Et > 125){
                                                            metEM++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
             }
          }
        }
      }
   }
   std::cout << "events before cut: " << beforeCut << std::endl;
   std::cout << "passGRL          : " << grl << std::endl;
   std::cout << "primaryVertex    : " << primaryVertex << std::endl;
   std::cout << "triggerMatch     : " << trigger << std::endl;
   std::cout << "Bad Muon         : " << badMuonNo << std::endl;
   std::cout << "Jet pass or      : " << jetPass << std::endl;
   std::cout << "bad jet          : " << badJetNo << std::endl;
   std::cout << "signal jet       : " << signalJetNo << std::endl;
   std::cout << "cosmic muon      : " << cosmicMuonNo << std::endl;
   std::cout << "baselineLepton   : " << baselineLepton << std::endl;
   std::cout << "signalLepton     : " << signalLepton << std::endl;
   std::cout << "channelEE        : " << channelEE << std::endl;
   std::cout << "bJetEE           : " << bJetEE << std::endl;
   std::cout << "fourJetEE        : " << fourJetEE << std::endl;
   std::cout << "metEE            : " << metEE << std::endl;
   std::cout << "channelMM        : " << channelMM << std::endl;
   std::cout << "bJetMM           : " << bJetMM << std::endl;
   std::cout << "fourJetMM        : " << fourJetMM << std::endl;
   std::cout << "metMM            : " << metMM << std::endl;
   std::cout << "channelEM        : " << channelEM << std::endl;
   std::cout << "bJetEM           : " << bJetEM << std::endl;
   std::cout << "fourJetEM        : " << fourJetEM << std::endl;
   std::cout << "metEM            : " << metEM << std::endl;
   
}
