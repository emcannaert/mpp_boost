


// system include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// new includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include  "DataFormats/PatCandidates/interface/PackedCandidate.h"

//#include "​DataFormats/​MuonReco/​interface/​MuonSelectors.h"

#include <algorithm>   

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <string>
using namespace reco;
typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZVector Vector;

class ttbarAnalyzer : public edm::EDAnalyzer 
{
public:
   explicit ttbarAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   double calc_mag(double px,double py, double pz);
   bool isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles);
   bool isgoodAK4(const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles);
   bool isgoodmuon(const float isMediumMuon, const float eta, const float pt);
   const reco::Candidate* parse_chain(const reco::Candidate* cand);

   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 
   edm::EDGetTokenT<std::vector<pat::Jet>> PUPPI_AK4Token_;
   edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

   TTree * tree;

   int eventnum = 0;
   int nfatjets = 0;
   int raw_nfatjets;
   int nlead_AK4 = 3;
   int nAK4_20_int = 0;
   int nfatjets_minus1;
   int nmuons = 0;
   int nevents = 0;
   int nPass_iso = 0;
   int nPass_hpt = 0;
   int nbtagged_AK4 = 0;
   int hpt_muons = 0;
   double gen_top_eta, gen_top_phi, gen_top_pt; 
   double gen_W_eta, gen_W_phi, gen_W_pt; 
   double had_gen_bW_dr;
   double top_had_AK8_dr[100],top_lept_AK8_dr[100], AK4_20_phi[100],AK4_20_eta[100],jet_PUPPI_mass[100];
   double AK4_20_px[100],AK4_20_py[100],AK4_20_pz[100],AK8_gen_cone_mass_lept[100];
   double AK4_20_px_mpp[100],AK4_20_py_mpp[100],AK4_20_pz_mpp[100],AK4_20_E[100];
   double top_lept_eta,top_lept_phi, top_lept_pt;
   double jet_pt[100], jet_eta[100], jet_mass[100], jet_dr[100], raw_jet_mass[100],raw_jet_pt[100],raw_jet_phi[100], muon_pt[100],bjet_pt[100];
   double jet_beta[100], beta_T[100], AK4_mass_20[100],AK4_mass_30[100],AK4_mass_50[100],AK4_mass_70[100],AK4_mass_100[100],AK4_mass_150[100];
   double bjet_eta[100],bjet_phi[100], AK8_b_deltaR[100][100], bjet_muon_dr[100],AK4_bjet_angle_20[100],AK4_Wq_angle_20[100];
   double AK4_jet_mass_20[100],AK4_m12_20[100],AK4_m13_20[100],AK4_m23_20[100],AK8_gen_cone_mass[100];
   int jet_ndaughters[100], jet_nAK4[100],jet_nAK4_20[100],jet_nAK4_30[100],jet_nAK4_50[100],jet_nAK4_70[100],jet_nAK4_100[100],jet_nAK4_150[100];
   double had_gen_bjet_pt, had_gen_bjet_energy, had_gen_bjet_eta, had_gen_bjet_phi;
   double lept_gen_bjet_pt, lept_gen_bjet_energy, lept_gen_bjet_eta, lept_gen_bjet_phi;
   double had_b_px_mpp,had_b_py_mpp,had_b_pz_mpp,q1_px_mpp,q1_py_mpp,q1_pz_mpp,q2_px_mpp, q2_py_mpp, q2_pz_mpp;
   double had_b_E_mpp, q1_E_mpp, q2_E_mpp;
   double muon_eta,muon_phi;
};


ttbarAnalyzer::ttbarAnalyzer(const edm::ParameterSet& iConfig)
{
   genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));

   muonToken_ =    consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
   PUPPI_AK4Token_ =  consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("PUPPI_AK4Collection"));
   fatJetToken_ =    consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetCollection"));
   //triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));

   edm::Service<TFileService> fs;      

   tree = fs->make<TTree>("tree", "tree");

   tree->Branch("nfatjets", &nfatjets, "nfatjets/I");

   tree->Branch("nAK4_20_int", &nAK4_20_int, "nAK4_20_int/I");

   tree->Branch("nevents", &nevents, "nevents/I");
   tree->Branch("nPass_iso", &nPass_iso, "nPass_iso/I");
   tree->Branch("nPass_hpt", &nPass_hpt, "nPass_hpt/I");

   tree->Branch("jet_ndaughters", jet_ndaughters, "jet_ndaughters[nfatjets]/I");
   //tree->Branch("nfatjets_minus1",&nfatjets_minus1, "nfatjets_minus1/I");
   tree->Branch("jet_nAK4", jet_nAK4, "jet_nAK4[nfatjets]/I");
   tree->Branch("jet_nAK4_20", jet_nAK4_20, "jet_nAK4_20[nfatjets]/I");
   tree->Branch("jet_nAK4_30", jet_nAK4_30, "jet_nAK4_30[nfatjets]/I");
   tree->Branch("jet_nAK4_50", jet_nAK4_50, "jet_nAK4_50[nfatjets]/I");
   tree->Branch("jet_nAK4_70", jet_nAK4_70, "jet_nAK4_70[nfatjets]/I");
   tree->Branch("jet_nAK4_100", jet_nAK4_100, "jet_nAK4_100[nfatjets]/I");
   tree->Branch("jet_nAK4_150", jet_nAK4_150, "jet_nAK4_150[nfatjets]/I");

   tree->Branch("beta_T", beta_T, "beta_T[nfatjets]/D");

   tree->Branch("AK4_mass_20", AK4_mass_20, "AK4_mass_20[nfatjets]/D");
   tree->Branch("AK4_mass_30", AK4_mass_30, "AK4_mass_30[nfatjets]/D");
   tree->Branch("AK4_mass_50", AK4_mass_50, "AK4_mass_50[nfatjets]/D");
   tree->Branch("AK4_mass_70", AK4_mass_70, "AK4_mass_70[nfatjets]/D");
   tree->Branch("AK4_mass_100", AK4_mass_100, "AK4_mass_100[nfatjets]/D");
   tree->Branch("AK4_mass_150", AK4_mass_150, "AK4_mass_150[nfatjets]/D");

   tree->Branch("nbtagged_AK4", &nbtagged_AK4, "nbtagged_AK4/I");
   tree->Branch("hpt_muons", &hpt_muons, "hpt_muons/I");

   tree->Branch("bjet_eta", bjet_eta, "bjet_eta[nbtagged_AK4]/D");
   tree->Branch("bjet_phi", bjet_phi, "bjet_phi[nbtagged_AK4]/D");

   tree->Branch("AK8_b_deltaR", AK8_b_deltaR, "AK8_b_deltaR[nfatjets][nbtagged_AK4]/D");
   tree->Branch("bjet_muon_dr", bjet_muon_dr, "bjet_muon_dr[nbtagged_AK4]/D");



   tree->Branch("muon_pt", muon_pt, "muon_pt[hpt_muons]/D");
   tree->Branch("bjet_pt", bjet_pt, "bjet_pt[nbtagged_AK4]/D");
   tree->Branch("jet_beta", jet_beta, "jet_beta[nfatjets]/D");

   tree->Branch("jet_pt", jet_pt, "jet_pt[nfatjets]/D");
   tree->Branch("jet_eta", jet_eta, "jet_eta[nfatjets]/D");
   tree->Branch("jet_mass", jet_mass, "jet_mass[nfatjets]/D");
   tree->Branch("raw_nfatjets",&raw_nfatjets,"raw_nfatjets/I");
   tree->Branch("raw_jet_mass",raw_jet_mass,"raw_jet_mass[raw_nfatjets]/D");
   tree->Branch("raw_jet_pt",raw_jet_pt,"raw_jet_pt[raw_nfatjets]/D");
   tree->Branch("raw_jet_phi",raw_jet_phi,"raw_jet_phi[raw_nfatjets]/D");

   tree->Branch("AK4_jet_mass_20",AK4_jet_mass_20,"AK4_jet_mass_20[nAK4_20_int]/D");


   tree->Branch("AK4_bjet_angle_20",AK4_bjet_angle_20,"AK4_bjet_angle_20[nAK4_20_int]/D");
   tree->Branch("AK4_Wq_angle_20",AK4_Wq_angle_20,"AK4_Wq_angle_20[nAK4_20_int]/D");

   tree->Branch("AK4_20_eta",AK4_20_eta,"AK4_20_eta[nAK4_20_int]/D");
   tree->Branch("AK4_20_phi",AK4_20_phi,"AK4_20_phi[nAK4_20_int]/D");

   tree->Branch("AK4_20_px",AK4_20_px,"AK4_20_px[nAK4_20_int]/D");
   tree->Branch("AK4_20_py",AK4_20_py,"AK4_20_py[nAK4_20_int]/D");
   tree->Branch("AK4_20_pz",AK4_20_pz,"AK4_20_pz[nAK4_20_int]/D");
   tree->Branch("AK4_20_E",AK4_20_E,"AK4_20_E[nAK4_20_int]/D");



   tree->Branch("AK4_m12_20",AK4_m12_20,"AK4_m12_20[nfatjets]/D");
   tree->Branch("AK4_m13_20",AK4_m13_20,"AK4_m13_20[nfatjets]/D");
   tree->Branch("AK4_m23_20",AK4_m23_20,"AK4_m23_20[nfatjets]/D");

   tree->Branch("gen_top_eta",&gen_top_eta,"gen_top_eta/D");
   tree->Branch("gen_top_phi",&gen_top_phi,"gen_top_phi/D");
   tree->Branch("gen_top_pt",&gen_top_pt,"gen_top_pt/D");


   tree->Branch("gen_W_eta",&gen_W_eta,"gen_W_eta/D");
   tree->Branch("gen_W_phi",&gen_W_phi,"gen_W_phi/D");
   tree->Branch("gen_W_pt",&gen_W_pt,"gen_W_pt/D");

   tree->Branch("top_had_AK8_dr",top_had_AK8_dr,"top_had_AK8_dr[nfatjets]/D");
   //tree->Branch("jet_dr", jet_dr, "jet_dr[nfatjets_minus1]/D");
   tree->Branch("top_lept_AK8_dr",top_lept_AK8_dr,"top_lept_AK8_dr[nfatjets]/D");
   tree->Branch("AK8_gen_cone_mass",AK8_gen_cone_mass,"AK8_gen_cone_mass[nfatjets]/D");


   tree->Branch("top_lept_eta",&top_lept_eta,"top_lept_eta/D");
   tree->Branch("top_lept_phi",&top_lept_phi,"top_lept_phi/D");
   tree->Branch("top_lept_pt",&top_lept_pt,"top_lept_pt/D");

   tree->Branch("had_gen_bjet_pt",&had_gen_bjet_pt,"had_gen_bjet_pt/D");
   tree->Branch("had_gen_bjet_energy",&had_gen_bjet_energy,"had_gen_bjet_energy/D");
   tree->Branch("had_gen_bjet_eta",&had_gen_bjet_eta,"had_gen_bjet_eta/D");
   tree->Branch("had_gen_bjet_phi",&had_gen_bjet_phi,"had_gen_bjet_phi/D");


   tree->Branch("lept_gen_bjet_pt",&lept_gen_bjet_pt,"lept_gen_bjet_pt/D");
   tree->Branch("lept_gen_bjet_energy",&lept_gen_bjet_energy,"lept_gen_bjet_energy/D");
   tree->Branch("lept_gen_bjet_eta",&lept_gen_bjet_eta,"lept_gen_bjet_eta/D");
   tree->Branch("lept_gen_bjet_phi",&lept_gen_bjet_phi,"lept_gen_bjet_phi/D");


   tree->Branch("had_gen_bW_dr",&had_gen_bW_dr,"had_gen_bW_dr/D");

   tree->Branch("jet_PUPPI_mass",jet_PUPPI_mass,"jet_PUPPI_mass[nfatjets]/D");

   //information of boosted gen particlres

   tree->Branch("had_b_px_mpp",&had_b_px_mpp,"had_b_px_mpp/D");
   tree->Branch("had_b_py_mpp",&had_b_py_mpp,"had_b_py_mpp/D");
   tree->Branch("had_b_pz_mpp",&had_b_pz_mpp,"had_b_pz_mpp/D");
   tree->Branch("had_b_E_mpp",&had_b_E_mpp,"had_b_E_mpp/D");

   tree->Branch("q1_px_mpp",&q1_px_mpp,"q1_px_mpp/D");
   tree->Branch("q1_py_mpp",&q1_py_mpp,"q1_py_mpp/D");
   tree->Branch("q1_pz_mpp",&q1_pz_mpp,"q1_pz_mpp/D");
   tree->Branch("q1_E_mpp",&q1_E_mpp,"q1_E_mpp/D");

   tree->Branch("q2_px_mpp",&q2_px_mpp,"q2_px_mpp/D");
   tree->Branch("q2_py_mpp",&q2_py_mpp,"q2_py_mpp/D");
   tree->Branch("q2_pz_mpp",&q2_pz_mpp,"q2_pz_mpp/D");
   tree->Branch("q2_E_mpp",&q2_E_mpp,"q2_E_mpp/D");

   tree->Branch("AK8_gen_cone_mass_lept",AK8_gen_cone_mass_lept,"AK8_gen_cone_mass_lept[nfatjets]/D");
   tree->Branch("AK8_muon_dr",AK8_muon_dr,"AK8_muon_dr[nfatjets]/D");

   tree->Branch("muon_eta",&muon_eta,"muon_eta/D");
   tree->Branch("muon_phi",&muon_phi,"muon_phi/D");

}


double ttbarAnalyzer::calc_mag(double px,double py, double pz)
{
   return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
}

bool ttbarAnalyzer::isgoodjet(const float pt, const float fatjet_sd_mass, const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles)
{
   if( (abs(eta) > 2.4) || (pt < 300.) ||(fatjet_sd_mass < 30.0)) return false;

   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0) || (MUF > 0.8) || (CEMF > 0.8)) 
      {
         return false;
      }
   else{ return true;}
}
bool ttbarAnalyzer::isgoodAK4(const float eta, const float NHF,const float NEMF, const size_t NumConst,const float CHF,const int CHM, const float MUF, const float CEMF, const int NumNeutralParticles)
{
   if( (eta > 2.4)) return false;

   if ((NHF>0.9) || (NEMF>0.9) || (NumConst<1) || (CHF<0.) || (CHM<0) || (MUF > 0.8) || (CEMF > 0.8)) 
      {
         return false;
      }
   else{ return true;}
}
bool ttbarAnalyzer::isgoodmuon(const float isMediumMuon, const float eta, const float pt)
{
  
  if((isMediumMuon) && (abs(eta) < 2.4)  && ( pt >45.) ) return true;

  return false;
}

const reco::Candidate* ttbarAnalyzer::parse_chain(const reco::Candidate* cand)
{  
   for (unsigned int iii=0; iii<cand->numberOfDaughters(); iii++)
   {
      if(cand->daughter(iii)->pdgId() == cand->pdgId()) return parse_chain(cand->daughter(iii));
   }
   return cand;
}


void ttbarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

////////////////Triggers//////////////////////////
/*
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    nevents++;

   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   bool pass = false;
   std::string trigname  = "HLT_Mu50";
   std::string trigname2 = "HLT_OldMu100_v";
   std::string trigname3 = "HLT_TkMu100_v";

   for (unsigned int i = 0; i < triggerBits->size(); ++i) 
   {
      const std::string name = names.triggerName(i);
      const bool accept = triggerBits->accept(i);
      if (  ((name.find(trigname) != std::string::npos) &&(accept)) )
        {
          pass =true;
          nPass_hpt++;
        }
      else if ( ((name.find(trigname2) != std::string::npos) &&(accept)) ) pass = true;
      else if ( ((name.find(trigname3) != std::string::npos) &&(accept)) ) pass = true;
   }
   
   if(!pass) return;
*/
   nevents++;

////////////////////Muons////////////////////////////////////////////

 edm::Handle<std::vector<pat::Muon>> muons;
 iEvent.getByToken(muonToken_, muons);
 hpt_muons = 0;
 for(auto iMuon = muons->begin(); iMuon != muons->end(); iMuon++) 
 {
    if(!(isgoodmuon(iMuon->isMediumMuon(),iMuon->eta(),iMuon->pt() )) || !(iMuon->passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoMedium|reco::Muon::PFIsoLoose)) )continue;
    muon_pt[hpt_muons] = iMuon->pt(); 
    muon_eta = iMuon->eta();
    muon_phi = iMuon->phi();
    hpt_muons++;
 }


 if(hpt_muons!=1)return;


////////////////AK4 Jets//////////////////////////////////////

   edm::Handle<std::vector<pat::Jet> > PUPPI_AK4Jets;
   iEvent.getByToken(PUPPI_AK4Token_, PUPPI_AK4Jets);
   nbtagged_AK4 = 0;
   double min_bjet_dr = 999.;
   for(auto iJet = PUPPI_AK4Jets->begin(); iJet != PUPPI_AK4Jets->end(); iJet++)
   {
      double bdisc = iJet->bDiscriminator("pfDeepCSVJetTags:probb") + iJet->bDiscriminator("pfDeepCSVJetTags:probbb");
      if( !(isgoodAK4(iJet->eta(), iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(), iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(), iJet->chargedEmEnergyFraction(), iJet->neutralMultiplicity()  )) || (bdisc<0.1522 ) || (iJet->pt() < 50.) ) continue;
      bjet_pt[nbtagged_AK4] = iJet->pt();
      bjet_phi[nbtagged_AK4] = iJet->phi();
      bjet_eta[nbtagged_AK4] = iJet->eta();

      bjet_muon_dr[nbtagged_AK4] = sqrt(pow(iJet->eta()-muon_eta,2) + pow(iJet->phi()-muon_phi,2));


      if (bjet_muon_dr[nbtagged_AK4] < min_bjet_dr) min_bjet_dr = bjet_muon_dr[nbtagged_AK4];
      nbtagged_AK4++;
   }         

   if ((nbtagged_AK4<1)) return; //|| (min_bjet_dr>0.5))return;


////////////////Gen Particles///////////////////////////////
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
   iEvent.getByToken(genPartToken_, genParticles);

   TLorentzVector q1;
   TLorentzVector q2;
   TLorentzVector gen_b_lept;
   TLorentzVector gen_b_had;

   TLorentzVector top_had;
   TLorentzVector top_lept;
   TLorentzVector gen_W_had;
  

   int nb = 0;
   int nq = 0;
   for (auto iM = genParticles->begin(); iM != genParticles->end(); ++iM)  //////////Over Gen Particles
   {  
      if((abs(iM->pdgId())  == 5) && (iM->isLastCopy())) nb++;
      if((abs(iM->pdgId()) < 6) && (abs(iM->mother()->pdgId()) == 24))
      {
        if      (nq == 0) q1.SetPxPyPzE(iM->px(),iM->py(),iM->pz(),iM->energy());
        else if (nq == 1) q2.SetPxPyPzE(iM->px(),iM->py(),iM->pz(),iM->energy());
        nq++;
      }
      else if ((abs(iM->pdgId()) == 6) && (iM->isLastCopy()))
      {
        int nq_    = 0;
        int nlept  = 0;
        for(unsigned int iii = 0; iii<iM->numberOfDaughters();iii++)
        {
          if (abs(iM->daughter(iii)->pdgId()) == 24)
          {
            const reco::Candidate* gen_W = parse_chain(iM->daughter(iii)->clone());
            for(unsigned int jjj = 0; jjj< gen_W->numberOfDaughters();jjj++)
            {
              if ( abs(gen_W->daughter(jjj)->pdgId()) < 6 ) nq_++;
              else if ( (abs(gen_W->daughter(jjj)->pdgId()) < 19) && (abs(gen_W->daughter(jjj)->pdgId()) > 10 ) ) nlept++;
              if (nq_ == 2) gen_W_had.SetPxPyPzE(gen_W->px(),gen_W->py(),gen_W->pz(),gen_W->energy());
            }
          }
        }

        if ( nq_ == 2)
        {
          for(unsigned int iii = 0; iii<iM->numberOfDaughters();iii++)
          {
            if (abs(iM->daughter(iii)->pdgId()) == 5)
            {
              had_gen_bjet_eta = iM->daughter(iii)->eta();
              had_gen_bjet_phi = iM->daughter(iii)->phi();
              had_gen_bjet_pt  = iM->daughter(iii)->pt();
              had_gen_bjet_energy = iM->daughter(iii)->energy();
              gen_b_had.SetPxPyPzE(iM->daughter(iii)->px(),iM->daughter(iii)->py(),iM->daughter(iii)->pz(),iM->daughter(iii)->energy());

            }
          }

///////////////////////////////////////////////////////////////////
          const reco::Candidate* hadronic_top = parse_chain(iM->clone());
          top_had.SetPxPyPzE(hadronic_top->px(),hadronic_top->py(),hadronic_top->pz(),hadronic_top->energy());
          gen_top_eta = top_had.Eta();
          gen_top_phi = top_had.Phi();
          gen_top_pt  = top_had.Pt(); 

          gen_W_eta = gen_W_had.Eta(); 
          gen_W_phi = gen_W_had.Phi();
          gen_W_pt  = gen_W_had.Pt();

        }
        else if (nlept == 2) 
          {
            for(unsigned int iii = 0; iii<iM->numberOfDaughters();iii++)
            {
              if (abs(iM->daughter(iii)->pdgId()) == 5)
              {
                lept_gen_bjet_eta = iM->daughter(iii)->eta();
                lept_gen_bjet_phi = iM->daughter(iii)->phi();
                lept_gen_bjet_pt  = iM->daughter(iii)->pt();
                lept_gen_bjet_energy = iM->daughter(iii)->energy();
                gen_b_lept.SetPxPyPzE(iM->daughter(iii)->px(),iM->daughter(iii)->py(),iM->daughter(iii)->pz(),iM->daughter(iii)->energy());
              }

            }
            top_lept.SetPxPyPzE(iM->px(),iM->py(),iM->pz(),iM->energy());
            top_lept_pt = top_lept.Pt();
            top_lept_phi = top_lept.Phi();
            top_lept_eta = top_lept.Eta();

          }
      }
   }

   had_gen_bW_dr = sqrt(pow(gen_b_had.Eta()-gen_W_had.Eta(),2)+pow(gen_b_had.Phi()-gen_W_had.Phi(),2));

///////////////////////////////////////////////////////////////
////////////////Jets//////////////////////////////////////
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);

   nfatjets = 0;
   raw_nfatjets = 0;

   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   { 
      if(nfatjets>5)continue;
      raw_jet_pt[raw_nfatjets]   = iJet->pt();
      raw_jet_phi[raw_nfatjets]  = iJet->phi();
      raw_jet_mass[raw_nfatjets] = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      raw_nfatjets++;

      if(!(iJet->isPFJet())) return;
      if(!isgoodjet(iJet->pt(),iJet->userFloat("ak8PFJetsPuppiSoftDropMass"), iJet->eta(),iJet->neutralHadronEnergyFraction(), iJet->neutralEmEnergyFraction(),iJet->numberOfDaughters(),iJet->chargedHadronEnergyFraction(),iJet->chargedMultiplicity(),iJet->muonEnergyFraction(),iJet->chargedEmEnergyFraction(),iJet->neutralMultiplicity())) continue;
      jet_pt[nfatjets]         = iJet->pt();
      jet_PUPPI_mass[nfatjets] = iJet->userFloat("ak8PFJetsPuppiSoftDropMass");
      jet_eta[nfatjets]        = iJet->eta();
      jet_mass[nfatjets]       = sqrt(pow(iJet->energy(),2) - pow(iJet->p(),2));
      jet_ndaughters[nfatjets] = iJet->numberOfDaughters();
      
      top_lept_AK8_dr[nfatjets] = sqrt(pow(top_lept.Eta()-iJet->eta(),2)+pow(top_lept.Phi()-iJet->phi(),2));
      top_had_AK8_dr[nfatjets]  = sqrt(pow(top_had.Eta()-iJet->eta(),2)+pow(top_had.Phi()-iJet->phi(),2));
      double gen_cone_px = 0;
      double gen_cone_py = 0;
      double gen_cone_pz = 0;
      double gen_cone_E  = 0;
    
      if((top_had_AK8_dr[nfatjets] < 0.5))
      {
          for (auto iM = genParticles->begin(); iM != genParticles->end(); ++iM) ////over gen particles
          {
            if( (sqrt(pow(iM->eta()-iJet->eta(),2)+pow(iM->phi()-iJet->phi(),2))<0.8 ) && (iM->status()==1 ) ) 
            {
              gen_cone_px+= iM->px();
              gen_cone_py+= iM->py();
              gen_cone_pz+= iM->pz();
              gen_cone_E += iM->energy();
            }
          }
      }      

      AK8_gen_cone_mass[nfatjets] = sqrt(pow(gen_cone_E,2)-pow(gen_cone_px,2)-pow(gen_cone_py,2)-pow(gen_cone_pz,2));

/*
      double gen_cone_lept_px = 0;
      double gen_cone_lept_py = 0;
      double gen_cone_lept_pz = 0;
      double gen_cone_lept_E  = 0;


      if((top_lept_AK8_dr[nfatjets] < 0.5))
      {
          for (auto iM = genParticles->begin(); iM != genParticles->end(); ++iM) ////over gen particles
          {
            if( (sqrt(pow(iM->eta()-iJet->eta(),2)+pow(iM->phi()-iJet->phi(),2))<0.8 ) && (iM->status()==1 ) ) 
            {
            gen_cone_lept_px+= iM->px();
            gen_cone_lept_py+= iM->py();
            gen_cone_lept_pz+= iM->pz();
            gen_cone_lept_E+= iM->energy();

            }
          }
      }      

      AK8_gen_cone_mass_lept[nfatjets] = sqrt(pow(gen_cone_lept_E,2)-pow(gen_cone_lept_px,2)-pow(gen_cone_lept_py,2)-pow(gen_cone_lept_pz,2));
*/


      //double max_dr = -999.;

      for(int iii = 0;iii<nbtagged_AK4;iii++)
      {
        AK8_b_deltaR[nfatjets][iii] = sqrt(pow(iJet->phi()-bjet_phi[iii],2)+pow(iJet->eta()-bjet_eta[iii],2));
          //if(sqrt(pow(iJet->phi()-bjet_phi[iii],2)+pow(iJet->eta()-bjet_eta[iii],2)) > max_dr) max_dr = sqrt(pow(iJet->phi()-bjet_phi[iii],2)+pow(iJet->eta()-bjet_eta[iii],2));
      }
      AK8_muon_dr[nfatjets] = sqrt(pow(iJet->eta()-muon_eta,2) + pow(iJet->phi()-muon_phi,2)) ; 
      //if(max_dr < 1.4) continue;
      ////////////////calc mpp beta
      Vector jet_axis(iJet->px()/iJet->p(),iJet->py()/iJet->p(),iJet->pz()/iJet->p());
      TLorentzVector jet_particles(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      double min_pp = 999999999.;
      double min_boost = 0.;

      for(int iii=0;iii<10000;iii++)
      {
         TLorentzVector jet_particles_ = jet_particles;
         double beta_cand = iii/10000.;
         jet_particles_.Boost(-beta_cand*jet_axis.X(),-beta_cand*jet_axis.Y(),-beta_cand*jet_axis.Z());
         if(abs( ( jet_particles_.Px()*iJet->px()+jet_particles_.Py()*iJet->py() +jet_particles_.Pz()*iJet->py() )/iJet->p() ) < min_pp) 
            {
               min_boost = beta_cand; 
               min_pp = abs( ( jet_particles_.Px()*iJet->px()+jet_particles_.Py()*iJet->py() +jet_particles_.Pz()*iJet->py() )/iJet->p() ) ;
            }
      }
      double beta_mag = min_boost;
      ///////////////
      jet_beta[nfatjets] = beta_mag;
      beta_T[nfatjets] = beta_mag*sin(iJet->theta());


      ///Recluster//
      Vector beta(beta_mag*iJet->px()/iJet->p(),beta_mag*iJet->py()/iJet->p(),beta_mag*iJet->pz()/iJet->p());

      std::vector<fastjet::PseudoJet> cands_;

      for (unsigned int i=0; i<iJet->numberOfDaughters();i++)
      {
         const reco::Candidate* iJ = iJet->daughter(i);
         const pat::PackedCandidate* cand_begin = (pat::PackedCandidate*) iJ;
         double puppiweight = cand_begin->puppiWeight();
         TLorentzVector w(puppiweight*iJ->px(),puppiweight*iJ->py(),puppiweight*iJ->pz(),puppiweight*iJ->energy());
         w.Boost(-beta.X(),-beta.Y(),-beta.Z());
         cands_.push_back(fastjet::PseudoJet(w.Px(),w.Py(),w.Pz(),w.E()));

      }

      double R = 0.4;
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
      fastjet::ClusterSequence cs_jet(cands_, jet_def); 
      std::vector<fastjet::PseudoJet> jetsFJ_jet = fastjet::sorted_by_E(cs_jet.inclusive_jets(0.0));
      
      q1.Boost(-beta.X(),-beta.Y(),-beta.Z());
      q2.Boost(-beta.X(),-beta.Y(),-beta.Z());
      gen_b_had.Boost(-beta.X(),-beta.Y(),-beta.Z());

      q1_px_mpp = q1.Px();
      q1_py_mpp = q1.Py();
      q1_pz_mpp = q1.Pz();
      q1_E_mpp  = q1.E();

      q2_px_mpp =  q2.Px();
      q2_py_mpp =  q2.Py();
      q2_pz_mpp =  q2.Pz();
      q2_E_mpp =   q2.E();

      had_b_px_mpp = gen_b_had.Px();
      had_b_py_mpp = gen_b_had.Py();
      had_b_pz_mpp = gen_b_had.Pz();
      had_b_E_mpp = gen_b_had.E();

      jet_nAK4[nfatjets]    = 0;
      jet_nAK4_20[nfatjets] = 0;
      jet_nAK4_30[nfatjets] = 0;
      jet_nAK4_50[nfatjets] = 0;
      jet_nAK4_70[nfatjets] = 0;
      jet_nAK4_100[nfatjets] = 0;
      jet_nAK4_150[nfatjets] = 0;
      nAK4_20_int = 0;

      double jet_px_20 =0, jet_py_20 =0,jet_pz_20 =0,jet_E_20 =0;
      double jet_px_30 =0, jet_py_30 =0,jet_pz_30 =0,jet_E_30 =0;
      double jet_px_50 =0, jet_py_50 =0,jet_pz_50 =0,jet_E_50 =0;
      double jet_px_70 =0, jet_py_70 =0,jet_pz_70 =0,jet_E_70 =0;
      double jet_px_100 =0, jet_py_100 =0,jet_pz_100 =0,jet_E_100 =0;
      double jet_px_150 =0, jet_py_150 =0,jet_pz_150 =0,jet_E_150 =0;

      double jet_px_12 = 0, jet_py_12=0,jet_pz_12 = 0,jet_E_12=0;
      double jet_px_13 = 0, jet_py_13=0,jet_pz_13 = 0,jet_E_13=0;
      double jet_px_23 = 0, jet_py_23=0,jet_pz_23 = 0,jet_E_23=0;

      for (auto i=jetsFJ_jet.begin(); i<jetsFJ_jet.end(); i++)                             
      {
         if (i->E() > 20.)
         {

            jet_nAK4_20[nfatjets]++;
            jet_px_20+= i->px(); jet_py_20+= i->py();jet_pz_20+= i->pz();jet_E_20+= i->E();

            if(nfatjets < 1)
            {
                TLorentzVector AK4_20(i->px(),i->py(),i->pz(),i->E());
                AK4_jet_mass_20[nAK4_20_int] = i->m(); //sqrt(pow(i->E(),2)-pow(i->px(),2) -pow(i->py(),2)-pow(i->pz(),2));
                AK4_bjet_angle_20[nAK4_20_int] = AK4_20.Angle(gen_b_had.Vect());//sqrt(pow(i->eta()-gen_bjet.Eta(),2)+ pow(i->phi()-gen_bjet.Phi(),2));
                AK4_Wq_angle_20[nAK4_20_int]   = std::min(AK4_20.Angle(q1.Vect()), AK4_20.Angle(q2.Vect()));
            
              if(nAK4_20_int == 0)
              {
                jet_px_12+= i->px();
                jet_py_12+= i->py();
                jet_pz_12+= i->pz();
                jet_E_12 += i->E();
                jet_px_13+= i->px();
                jet_py_13+= i->py();
                jet_pz_13+= i->pz();
                jet_E_13 += i->E();
              }
              else if(nAK4_20_int == 1)
              {
                jet_px_12+= i->px();
                jet_py_12+= i->py();
                jet_pz_12+= i->pz();
                jet_E_12 += i->E();
                jet_px_23+= i->px();
                jet_py_23+= i->py();
                jet_pz_23+= i->pz();
                jet_E_23 += i->E();              
              }
              else if(nAK4_20_int == 2)
              {
                jet_px_13+= i->px();
                jet_py_13+= i->py();
                jet_pz_13+= i->pz();
                jet_E_13 += i->E();
                jet_px_23+= i->px();
                jet_py_23+= i->py();
                jet_pz_23+= i->pz();
                jet_E_23 += i->E();                
              }
              AK4_20_phi[nAK4_20_int] = i->phi();
              AK4_20_eta[nAK4_20_int] = i->eta();
              AK4_20_px[nAK4_20_int]  = i->px();
              AK4_20_py[nAK4_20_int]  = i->py();
              AK4_20_pz[nAK4_20_int]  = i->pz();
              AK4_20_E[nAK4_20_int]   = i->E();
              nAK4_20_int++;
            }
         }
         if(i->E() > 30.) 
         {
            jet_nAK4_30[nfatjets]++;
            jet_px_30+= i->px();jet_py_30+= i->py();jet_pz_30+= i->pz();jet_E_30+= i->E();
         }
         if (i->E() > 50.)
         {
            jet_nAK4_50[nfatjets]++;
            jet_px_50+= i->px();jet_py_50+= i->py();jet_pz_50+= i->pz();jet_E_50+= i->E();
         }
         if (i->E() > 70.)
         {
            jet_nAK4_70[nfatjets]++;
            jet_px_70+= i->px();jet_py_70+= i->py();jet_pz_70+= i->pz();jet_E_70+= i->E();
         } 
         if (i->E() > 100.)
         {
            jet_nAK4_100[nfatjets]++;
            jet_px_100+= i->px();jet_py_100+= i->py();jet_pz_100+= i->pz();jet_E_100+= i->E();
         }    
         if (i->E() > 150.)
         {
            jet_nAK4_150[nfatjets]++;
            jet_px_150+= i->px();jet_py_150+= i->py();jet_pz_150+= i->pz();jet_E_150+= i->E();
         }       
         jet_nAK4[nfatjets]++;

         AK4_m12_20[nfatjets] = 0;
         AK4_m12_20[nfatjets] = 0;
         AK4_m12_20[nfatjets] = 0;
         if(nAK4_20_int>2)
         {
          AK4_m12_20[nfatjets]=sqrt(pow(jet_E_12,2)-pow(jet_px_12,2)-pow(jet_py_12,2)-pow(jet_pz_12,2));
          AK4_m13_20[nfatjets]=sqrt(pow(jet_E_13,2)-pow(jet_px_13,2)-pow(jet_py_13,2)-pow(jet_pz_13,2));
          AK4_m23_20[nfatjets]=sqrt(pow(jet_E_23,2)-pow(jet_px_23,2)-pow(jet_py_23,2)-pow(jet_pz_23,2));
         }
         AK4_mass_20[nfatjets] = sqrt(pow(jet_E_20,2)-pow(jet_px_20,2) - pow(jet_py_20,2)-pow(jet_pz_20,2));
         AK4_mass_30[nfatjets] = sqrt(pow(jet_E_30,2)-pow(jet_px_30,2) - pow(jet_py_30,2)-pow(jet_pz_30,2));
         AK4_mass_50[nfatjets] = sqrt(pow(jet_E_50,2)-pow(jet_px_50,2) - pow(jet_py_50,2)-pow(jet_pz_50,2));
         AK4_mass_70[nfatjets] = sqrt(pow(jet_E_70,2)-pow(jet_px_70,2) - pow(jet_py_70,2)-pow(jet_pz_70,2));
         AK4_mass_100[nfatjets] = sqrt(pow(jet_E_100,2)-pow(jet_px_100,2) - pow(jet_py_100,2)-pow(jet_pz_100,2));
         AK4_mass_150[nfatjets] = sqrt(pow(jet_E_150,2)-pow(jet_px_150,2) - pow(jet_py_150,2)-pow(jet_pz_150,2));
      }      
      nfatjets++;
   }

   if(nfatjets < 1) return;

   eventnum++;
   tree->Fill();
}   
DEFINE_FWK_MODULE(ttbarAnalyzer);



