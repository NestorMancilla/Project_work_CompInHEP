#!/usr/bin/env python

import ROOT
import numpy as np
import datetime

ROOT.gROOT.SetBatch(True)
ROOT.gRandom.SetSeed(42)


# cross section. Signal and background have different xsec format.
def get_xsec(file_path, is_signal):
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        raise OSError(f"Could not open file {file_path}")

    xsec_param = f.Get("xsec")
    if not xsec_param:
        raise ValueError(f"File {file_path} missing 'xsec' parameter!")

    try:
        # First try to get the value as a TParameter (works for both signal and background)
        xsec = xsec_param.GetVal()
    except AttributeError:
        try:
            xsec = float(xsec_param)
        except (TypeError, ValueError):
            raise ValueError(f"File {file_path} has invalid 'xsec' type - neither TParameter nor float!")
    #if not is_signal: # the xsec from the background is not in fb
        #xsec = xsec * 1e9 
    if is_signal:
        xsec = 62 * 1e3 # fb, numbers by Sami
    else:
        xsec = 924 * 1e3 # fb, numbers by Sami
    f.Close()
    return xsec

# 2a smearing theta -> eta, phi and the pT.
def smear_muon(pt, eta, phi):
    smeared_pt = pt * (1.0 + ROOT.gRandom.Gaus(0, 0.01))
    theta = 2.0 * np.arctan(np.exp(-eta))
    theta_smeared = theta + ROOT.gRandom.Gaus(0, 0.002)
    smeared_eta = -np.log(np.tan(theta_smeared / 2.0))
    smeared_phi = phi + ROOT.gRandom.Gaus(0, 0.002)
    return smeared_pt, smeared_eta, smeared_phi

# 2b check for isolation.
def is_isolated(mu_eta, mu_phi, pi_pts, pi_etas, pi_phis):
    sum_pi_pt = 0.0
    for pi_pt, pi_eta, pi_phi in zip(pi_pts, pi_etas, pi_phis):
        deta = mu_eta - pi_eta
        dphi = ROOT.TVector2.Phi_mpi_pi(mu_phi - pi_phi)
        dr = (deta**2 + dphi**2)**0.5
        if dr < 0.3:
            sum_pi_pt += pi_pt
    return sum_pi_pt < 1.5

def process_file(file_path, is_signal, hist, label):
    xsec = get_xsec(file_path, is_signal)
    f = ROOT.TFile.Open(file_path)
    tree = f.Get("Events")
    n_passed = 0
    n_total = tree.GetEntries()

    for event in tree:
        mu_pts, mu_etas, mu_phis, mu_Es = [], [], [], []

        for i in range(len(event.mu_pt)):
            smeared_pt, smeared_eta, smeared_phi = smear_muon(event.mu_pt[i], event.mu_eta[i], event.mu_phi[i])
            if smeared_pt > 30 and is_isolated(smeared_eta, smeared_phi, event.pi_pt, event.pi_eta, event.pi_phi):
                mu_pts.append(smeared_pt)
                mu_etas.append(smeared_eta)
                mu_phis.append(smeared_phi)
                mu_Es.append(event.mu_E[i])

        if len(mu_pts) >= 2:
            n_passed += 1
            mu1 = ROOT.TLorentzVector()
            mu2 = ROOT.TLorentzVector()
            mu1.SetPtEtaPhiE(mu_pts[0], mu_etas[0], mu_phis[0], mu_Es[0])
            mu2.SetPtEtaPhiE(mu_pts[1], mu_etas[1], mu_phis[1], mu_Es[1])

            mumu_mass = (mu1 + mu2).M()
            weight = xsec / n_total
            hist.Fill(mumu_mass, weight)
    print(f"{label}: the cross section is {xsec}")
    print(f"{label}: {n_passed} events passed selection out of {n_total}")
    f.Close()

def main():
    hist_sig = ROOT.TH1F("h_sig", "Invariant Mass;M_{#mu#mu} [GeV];Events (fb)", 100, 60, 120)
    hist_bkg = ROOT.TH1F("h_bkg", "Invariant Mass;M_{#mu#mu} [GeV];Events (fb)", 100, 60, 120)

    process_file("zmm_signal_10mil.root", is_signal=True, hist=hist_sig, label="Signal")
    process_file("ttbar_bkg_10mil.root", is_signal=False, hist=hist_bkg, label="Background")

    hist_total = hist_sig.Clone("h_total")
    hist_total.Add(hist_bkg)

    canvas = ROOT.TCanvas("c", "Invariant Mass", 800, 600)
    hist_sig.SetLineColor(ROOT.kRed)
    hist_sig.SetLineStyle(10) # 10 dashed line to distinguish from total hist
    hist_sig.SetLineWidth(2)
    hist_bkg.SetLineColor(ROOT.kBlue)
    hist_bkg.SetLineWidth(2)
    hist_total.SetLineColor(ROOT.kBlack)
    hist_total.SetLineStyle(2)  # 2 = dashed line
    hist_total.SetLineWidth(2)
    
    hist_total.Draw("HIST")
    hist_sig.Draw("HIST SAME")
    hist_bkg.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.65, 0.75, 0.85, 0.85)
    legend.AddEntry(hist_sig, "Signal", "l")
    legend.AddEntry(hist_bkg, "Background", "l")
    legend.AddEntry(hist_total, "Total", "l")
    legend.SetBorderSize(0)
    legend.Draw()

    canvas.SaveAs("Invariant_mumu_mass.pdf")
    canvas.SaveAs("Invariant_mumu_mass.png")

    f_out = ROOT.TFile("output.root", "RECREATE")
    hist_sig.Write()
    hist_bkg.Write()
    hist_total.Write()

    # timestamp
    days = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
    now = datetime.datetime.now()
    m = f"Produced: {days[now.weekday()]} {now}"
    timestamp = ROOT.TNamed(m, "")
    timestamp.Write()

    # --- 2c) Gaussian + linear background fit and significance (with Voigtian) ---
    try:
        from array import array
        # Define fit range
        fitMin, fitMax = 60.0, 120.0
        # Sideband background fit (pol1) outside [80,100]
        fbkg = ROOT.TF1("fbkg", "pol1", fitMin, fitMax)
        hist_total.Fit(fbkg, "R", "", fitMin, 80.0)
        hist_total.Fit(fbkg, "R+", "", 100.0, fitMax)

        # Gaussian + linear background fit
        fGauss = ROOT.TF1("fGauss",
            "[0]*TMath::Gaus(x,[1],[2]) + [3] + [4]*x",
            fitMin, fitMax)
        # initial parameters: amplitude, mean, sigma, background const, slope
        fGauss.SetParameters(
            hist_total.GetMaximum(),  # amplitude
            91.15,                   # mean
            2.5,                     # sigma
            0.1*hist_total.GetMaximum(),
            0.0                     # slope
        )
        hist_total.Fit(fGauss, "R+")
        fGauss.SetLineColor(ROOT.kMagenta)
        fGauss.Draw("same")

        # Compute integrals in signal window [80,100]
        winMin, winMax = 80.0, 100.0

        # Draw and save combined plot
        c3 = ROOT.TCanvas("c3", "Normalized invariant mass", 800, 600)
        hist_total.SetLineColor(ROOT.kBlack)
        hist_sig.SetLineColor(ROOT.kRed)
        hist_bkg.SetLineColor(ROOT.kBlue)
        hist_total.Draw("hist")
        hist_sig.Draw("hist SAME")
        hist_bkg.Draw("hist SAME")
        fGauss.Draw("same")

        # Voigtian + linear background fit
        fVoigt = ROOT.TF1("fVoigt",
            "[0]*TMath::Voigt(x-[1],[2],[3]) + [4] + [5]*x",
            fitMin, fitMax)
        # initial parameters: norm, mean approx 91.15, sigma=2.5, gamma=2.5, background const and slope
        fVoigt.SetParameters(
            hist_total.GetMaximum(),  # normalization
            91.15,                    # mean
            2.5,                      # sigma (Gaussian width)
            2.5,                      # gamma (Lorentzian width)
            0.1 * hist_total.GetMaximum(),
            0.0                       # linear slope
        )
        hist_total.Fit(fVoigt, "R+")
        fVoigt.SetLineColor(ROOT.kGreen+2)
        fVoigt.Draw("same")
        # add Voigtian to legend
        legend = ROOT.TLegend(0.65, 0.75, 0.85, 0.85)
        legend.AddEntry(hist_sig, "Signal", "l")
        legend.AddEntry(hist_bkg, "Background", "l")
        legend.AddEntry(hist_total, "Total", "l")
        legend.AddEntry(fGauss, "Gaussian + linear", "l")
        legend.AddEntry(fVoigt, "Voigtian + linear", "l")
        legend.SetBorderSize(0)
        legend.Draw()

        # Compute integrals and significance for both fits
        NB = fbkg.Integral(winMin, winMax)
        fullIntGauss = fGauss.Integral(winMin, winMax)
        NS_gauss = fullIntGauss - NB
        signif_gauss = NS_gauss/np.sqrt(NB) if NB > 0 else 0.0
        fullIntVoigt = fVoigt.Integral(winMin, winMax)
        NS_voigt = fullIntVoigt - NB
        signif_voigt = NS_voigt/np.sqrt(NB) if NB > 0 else 0.0

        # Determine peak positions
        xPeakGauss = fGauss.GetMaximumX(winMin, winMax)
        xPeakVoigt = fVoigt.GetMaximumX(winMin, winMax)

        # Print results for both fits
        L_req = (5.0/signif_gauss)**2 if signif_gauss > 0 else 0.0
        days = L_req/50.0*365.0 if signif_gauss > 0 else 0.0
        print(f"2c) Invariant mass window [{winMin},{winMax}] GeV")
        print(f"  Gaussian: NS={NS_gauss:.3f} fb, NB={NB:.3f} fb, significance={signif_gauss:.3f} sigma, peak={xPeakGauss:.3f} GeV")
        print(f"  Voigtian: NS={NS_voigt:.3f} fb, NB={NB:.3f} fb, significance={signif_voigt:.3f} sigma, peak={xPeakVoigt:.3f} GeV")
        
        c3.SaveAs("sum_and_fit.pdf")
    except Exception as e:
        print("2c block failed:", e)

    f_out.Close()

if __name__ == "__main__":
    main()