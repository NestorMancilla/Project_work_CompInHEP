#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace Pythia8;

int main() {
  Pythia pythia;

  // ---------- beam & tune ----------
  pythia.readString("Beams:idA = 2212");             // proton
  pythia.readString("Beams:idB = 2212");             // proton
  pythia.readString("Beams:eCM = 13600.");           // 13.6 TeV

  // ---------- physics process ----------
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 13 -13");        // Z → μ+ μ–

  // ---------- initialise ----------
  pythia.init();

  // ---------- trigger definition ----------
  const double pTmin = 20.0;
  const double etaMax = 2.1;

  // ---------- ROOT TTree setup ----------
    TFile *f = new TFile("zmm_signal.root", "RECREATE");
    TTree *tree = new TTree("Events", "PYTHIA8 Z->mumu events");
  
  // ---------- Variables to store in TTree ----------

  // scalars
  UInt_t   run = 1;
  // Float_t evt_w = 1.0;                                  // event weight

  tree->Branch("run",&run,"run/i");

  // muon vectors
  std::vector<float> mu_pt, mu_eta, mu_phi, mu_E;
  // muon charge and PDG ID for OS pair selection and MC truth checks
  std::vector<int> mu_charge;
  std::vector<float> pi_pT, pi_eta, pi_phi;      // charged pions

  tree->Branch("mu_pt",&mu_pt);
  tree->Branch("mu_eta",&mu_eta);
  tree->Branch("mu_phi",&mu_phi);
  tree->Branch("mu_E",&mu_E);
  tree->Branch("mu_charge", &mu_charge);

  tree->Branch("pi_pT",  &pi_pT);
  tree->Branch("pi_eta", &pi_eta);
  tree->Branch("pi_phi", &pi_phi);

    // ------------------------------------------------------------------
  // 3.  Event loop
  // ------------------------------------------------------------------

    const int nTotal = 1000000; //10 000 000

    for (int iEvent = 0; iEvent < nTotal; ++iEvent) {
        if (!pythia.next()) continue;
    
        mu_pt.clear(); mu_eta.clear(); mu_phi.clear(); mu_E.clear();
        mu_charge.clear();
        pi_pT.clear(); pi_eta.clear(); pi_phi.clear();

        int passHLTmu = 0;
    
        for (int i = 0; i < pythia.event.size(); ++i) {
          const Particle& p = pythia.event[i];
          if (!p.isFinal()) continue;

          if (p.pT() > pTmin && std::abs(p.eta()) < etaMax){
            // ---- store muons ----
            if (std::abs(p.id()) == 13) {
                mu_pt .push_back(p.pT());
                mu_eta.push_back(p.eta());
                mu_phi.push_back(p.phi());
                mu_E  .push_back(p.e());
                mu_charge.push_back(p.charge());
                // count how many muons satisfy the HLT kinematics
                ++passHLTmu;
            }

            // ---- store charged pions ----
            if (std::abs(p.id()) == 211) {
                pi_pT .push_back(p.pT());
                pi_eta.push_back(p.eta());
                pi_phi.push_back(p.phi());
            }
          }
        }
        if (passHLTmu >= 2) {          // ≥2 muons pass HLT → keep event
            tree->Fill();
        }
    }

  // ------------------------------------------------------------------
  // 4.  Cross-section & per-event weight
  // ------------------------------------------------------------------
  pythia.stat();
  // sigmaGen() returns cross-section in mb, convert to fb
  Float_t xsec = pythia.info.sigmaGen() * 1e12;    // fb
  Float_t weight   = xsec / tree->GetEntries(); // fb
  
  // store metadata
  TParameter<float>("weight", weight).Write();
  TParameter<float>("xsec",    xsec).Write();
  
  std::cout << "signal : tried " << nTotal
            << ", kept " << tree->GetEntries()
            << ", xsec = " << xsec << " fb\n";
  
  tree->Write();
  f->Close();
  return 0;
}