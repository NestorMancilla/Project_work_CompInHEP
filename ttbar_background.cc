#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
using namespace Pythia8;

int main() {
  Pythia pythia;
  pythia.readString("Beams:eCM = 13600.");
  pythia.readString("Top:gg2ttbar = on");
  pythia.init();

  TFile *outFile = new TFile("ttbar_bkg.root", "RECREATE");
  TTree *tree = new TTree("Events", "ttbar events");

  Int_t event_id;
  std::vector<float> mu_pt, mu_eta, mu_phi, mu_e;
  std::vector<float> pi_pt, pi_eta, pi_phi, pi_e;

  tree->Branch("event_id", &event_id);
  tree->Branch("mu_pt", &mu_pt);
  tree->Branch("mu_eta", &mu_eta);
  tree->Branch("mu_phi", &mu_phi);
  tree->Branch("mu_e", &mu_e);
  tree->Branch("pi_pt", &pi_pt);
  tree->Branch("pi_eta", &pi_eta);
  tree->Branch("pi_phi", &pi_phi);
  tree->Branch("pi_e", &pi_e);

  int nEvents = 500000;
  for (int i = 0; i < nEvents; ++i) {
    if (!pythia.next()) continue;

    mu_pt.clear(); mu_eta.clear(); mu_phi.clear(); mu_e.clear();
    pi_pt.clear(); pi_eta.clear(); pi_phi.clear(); pi_e.clear();

  // id
  // 6 top
  // 13 mu-
  // 14 nu_mu-
  // 5 b

    for (int j = 0; j < pythia.event.size(); ++j) {
      const Particle& p = pythia.event[j];
      if (p.idAbs() == 13 && p.isFinal() && p.pT() > 20. && std::abs(p.eta()) < 2.1) {
        mu_pt.push_back(p.pT());
        mu_eta.push_back(p.eta());
        mu_phi.push_back(p.phi());
        mu_e.push_back(p.e());
      }
      if (p.idAbs() == 211 && p.isFinal() && p.charge() != 0) {
        pi_pt.push_back(p.pT());
        pi_eta.push_back(p.eta());
        pi_phi.push_back(p.phi());
        pi_e.push_back(p.e());
      }
    }

    if (mu_pt.size() < 2) continue; // trigger condition

    event_id = i;
    tree->Fill();
  }

  TParameter<double> xsec("xsec_fb", pythia.info.sigmaGen() * 1e3); // pb to fb
  xsec.Write();

  tree->Write();
  outFile->Close();
  pythia.stat();
  return 0;
}

