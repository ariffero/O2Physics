// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "PWGUD/DataModel/UDTables.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TMath.h"
#include "TRandom3.h"

namespace dimu
{
// dimuon
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(E, energy, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PhiAv, phiAv, float);
DECLARE_SOA_COLUMN(PhiCh, phiCh, float);
// tracks positive (p) and negative (n)
DECLARE_SOA_COLUMN(Ep, energyp, float);
DECLARE_SOA_COLUMN(Pxp, pxp, float);
DECLARE_SOA_COLUMN(Pyp, pyp, float);
DECLARE_SOA_COLUMN(Pzp, pzp, float);
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(En, energyn, float);
DECLARE_SOA_COLUMN(Pxn, pxn, float);
DECLARE_SOA_COLUMN(Pyn, pyn, float);
DECLARE_SOA_COLUMN(Pzn, pzn, float);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
// zn
DECLARE_SOA_COLUMN(Tzna, tzna, float);
DECLARE_SOA_COLUMN(Ezna, ezna, float);
DECLARE_SOA_COLUMN(Tznc, tznc, float);
DECLARE_SOA_COLUMN(Eznc, eznc, float);
DECLARE_SOA_COLUMN(Nclass, nclass, int);
} // namespace dimu

namespace o2::aod
{
DECLARE_SOA_TABLE(DiMu, "AOD", "DIMU",
                  dimu::M, dimu::E, dimu::Px, dimu::Py, dimu::Pz, dimu::Pt, dimu::Rap, dimu::Phi,
                  dimu::PhiAv, dimu::PhiCh,
                  dimu::Ep, dimu::Pxp, dimu::Pyp, dimu::Pzp, dimu::Ptp, dimu::Etap, dimu::Phip,
                  dimu::En, dimu::Pxn, dimu::Pyn, dimu::Pzn, dimu::Ptn, dimu::Etan, dimu::Phin,
                  dimu::Tzna, dimu::Ezna, dimu::Tznc, dimu::Eznc, dimu::Nclass);
}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// defining constants
double mMu = 0.10566; // mass of muon


// constants used in the track selection
const float kRAbsMin = 17.6;
const float kRAbsMid = 26.5;
const float kRAbsMax = 89.5;
const float kPDca1 = 200.;
const float kPDca2 = 200.;
const float kEtaMin = -4.0;
const float kEtaMax = -2.5;
const float kPtMin = 0.;

struct fwdMuonsUPC {

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteCandidatesFwd = soa::Join<CandidatesFwd, o2::aod::McTrackLabels>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::McTrackLabels>;

  Produces<o2::aod::DiMu> dimuSel;

  // defining histograms using histogram registry: different histos for the different process functions
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry reg0n0n{"reg0n0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXn0n{"regXn0n", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry regXnXn{"regXnXn", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry McGenRegistry{"McGenRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry McCorrReg{"McCorrReg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  
  

  // CONFIGURABLES
  // pT of muon pairs
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0., "lower limit in pT histo"};
  Configurable<float> highPt{"highPt", 2, "upper limit in pT histo"};
  // mass of muon pairs
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in mass histo"};
  Configurable<float> lowMass{"lowMass", 0., "lower limit in mass histo"};
  Configurable<float> highMass{"highMass", 10., "upper limit in mass histo"};
  // eta of muon pairs
  Configurable<int> nBinsEta{"nBinsEta", 600, "N bins in eta histo"};
  Configurable<float> lowEta{"lowEta", -10., "lower limit in eta histo"};
  Configurable<float> highEta{"highEta", -2., "upper limit in eta histo"};
  // rapidity of muon pairs
  Configurable<int> nBinsRapidity{"nBinsRapidity", 250, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4.5, "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", -2., "upper limit in rapidity histo"};
  // phi of muon pairs
  Configurable<int> nBinsPhi{"nBinsPhi", 600, "N bins in phi histo"};
  Configurable<float> lowPhi{"lowPhi", -TMath::Pi(), "lower limit in phi histo"};
  Configurable<float> highPhi{"highPhi", TMath::Pi(), "upper limit in phi histo"};
  // pT of single muons
  Configurable<int> nBinsPtSingle{"nBinsPtSingle", 500, "N bins in pT histo single muon"};
  Configurable<float> lowPtSingle{"lowPtSingle", 0., "lower limit in pT histo single muon"};
  Configurable<float> highPtSingle{"highPtSingle", 2., "upper limit in pT histo single muon"};
  // eta of single muons
  Configurable<int> nBinsEtaSingle{"nBinsEtaSingle", 250, "N bins in eta histo single muon"};
  Configurable<float> lowEtaSingle{"lowEtaSingle", -4.5, "lower limit in eta histo single muon"};
  Configurable<float> highEtaSingle{"highEtaSingle", -2., "upper limit in eta histo single muon"};
  // phi of single muons
  Configurable<int> nBinsPhiSingle{"nBinsPhiSingle", 600, "N bins in phi histo single muon"};
  Configurable<float> lowPhiSingle{"lowPhiSingle", -TMath::Pi(), "lower limit in phi histo single muon"};
  Configurable<float> highPhiSingle{"highPhiSingle", TMath::Pi(), "upper limit in phi histo single muon"};
  // ZDC
  Configurable<int> nBinsZDCen{"nBinsZDCen", 200, "N bins in ZN energy"};
  Configurable<float> lowEnZN{"lowEnZN", -50., "lower limit in ZN energy histo"};
  Configurable<float> highEnZN{"highEnZN", 250., "upper limit in ZN energy histo"};

  // configuarble rapidity cuts
  //Configurable<float> yCutLow{"yCutLow", -4, "Lower cut in pair rapidity"};
  //Configurable<float> yCutUp{"yCutUp", -2.5, "Upper cut in pair rapidity"};

  void init(InitContext& context)
  {
    // binning of pT axis fr fit
    std::vector<double> ptFitBinning = {
      0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
      0.11, 0.12, 0.13, 0.14, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50,
      0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50,
      3.00, 3.50};

    // axis
    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisPtFit = {ptFitBinning, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisMass{nBinsMass, lowMass, highMass, "m_{#mu#mu} GeV/#it{c}^{2}"};
    const AxisSpec axisEta{nBinsEta, lowEta, highEta, "#eta"};
    const AxisSpec axisRapidity{nBinsRapidity, lowRapidity, highRapidity, "Rapidity"};
    const AxisSpec axisPhi{nBinsPhi, lowPhi, highPhi, "#varphi"};
    const AxisSpec axisPtSingle{nBinsPtSingle, lowPtSingle, highPtSingle, "#it{p}_{T}_{ trk} GeV/#it{c}"};
    const AxisSpec axisTimeZN{200, -10, 10, "ZDC time (ns)"};
    const AxisSpec axisEnergyZNA{nBinsZDCen, lowEnZN, highEnZN, "ZNA energy (TeV)"};
    const AxisSpec axisEnergyZNC{nBinsZDCen, lowEnZN, highEnZN, "ZNC energy (TeV)"};
    const AxisSpec axisEtaSingle{nBinsEtaSingle, lowEtaSingle, highEtaSingle, "#eta_{trk}"};
    const AxisSpec axisPhiSingle{nBinsPhiSingle, lowPhiSingle, highPhiSingle, "#varphi_{trk}"};

    // histos
    // data and reco MC
    registry.add("hMass", "Ivariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPt", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPt});
    registry.add("hPtFit", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPtFit});
    registry.add("hEta", "Pseudorapidty of muon pairs;;#counts", kTH1D, {axisEta});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});
    registry.add("hPhi", "#varphi of muon pairs;;#counts", kTH1D, {axisPhi});
    registry.add("hCharge", "Charge;#it{charge};;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hContrib", "hContrib;;#counts", kTH1D, {{6, -0.5, 5.5}});
    registry.add("hEvSign", "Sum of the charges of all the tracks in each event;;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hPtTrkPos", "Pt of positive muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hPtTrkNeg", "Pt of negative muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hEtaTrkPos", "#eta of positive muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hEtaTrkNeg", "#eta of negative muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hPhiTrkPos", "#varphi of positive muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hPhiTrkNeg", "#varphi of negative muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hSameSign", "hSameSign;;#counts", kTH1D, {{6, -0.5, 5.5}});
    registry.add("hPhiCharge", "#phi #it{charge}", kTH1D,{axisPhi});
    registry.add("hPhiAverage", "#phi #it{average}", kTH1D,{axisPhi});

    // data
    registry.add("hTimeZNA", "ZNA Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hTimeZNC", "ZNC Times;;#counts", kTH1D, {axisTimeZN});
    registry.add("hEnergyZN", "ZNA vs ZNC energy", kTH2D, {axisEnergyZNA, axisEnergyZNC});

    reg0n0n.add("hMass", "Ivariant mass of muon pairs - 0n0n;;#counts", kTH1D, {axisMass});
    reg0n0n.add("hPt", "Transverse momentum of muon pairs - 0n0n;;#counts", kTH1D, {axisPt});
    reg0n0n.add("hEta", "Pseudorapidty of muon pairs - 0n0n;;#counts", kTH1D, {axisEta});
    reg0n0n.add("hRapidity", "Rapidty of muon pairs - 0n0n;;#counts", kTH1D, {axisRapidity});
    reg0n0n.add("hPtFit", "Transverse momentum of muon pairs - 0n0n;;#counts", kTH1D, {axisPtFit});

    regXn0n.add("hMass", "Ivariant mass of muon pairs - Xn0n;;#counts", kTH1D, {axisMass});
    regXn0n.add("hPt", "Transverse momentum of muon pairs - Xn0n;;#counts", kTH1D, {axisPt});
    regXn0n.add("hEta", "Pseudorapidty of muon pairs - Xn0n;;#counts", kTH1D, {axisEta});
    regXn0n.add("hRapidity", "Rapidty of muon pairs - Xn0n;;#counts", kTH1D, {axisRapidity});
    regXn0n.add("hPtFit", "Transverse momentum of muon pairs - Xn0n;;#counts", kTH1D, {axisPtFit});

    regXnXn.add("hMass", "Ivariant mass of muon pairs - XnXn;;#counts", kTH1D, {axisMass});
    regXnXn.add("hPt", "Transverse momentum of muon pairs - XnXn;;#counts", kTH1D, {axisPt});
    regXnXn.add("hEta", "Pseudorapidty of muon pairs - XnXn;;#counts", kTH1D, {axisEta});
    regXnXn.add("hRapidity", "Rapidty of muon pairs - XnXn;;#counts", kTH1D, {axisRapidity});
    regXnXn.add("hPtFit", "Transverse momentum of muon pairs - XnXn;;#counts", kTH1D, {axisPtFit});

    // gen MC
    McGenRegistry.add("hMass", "Ivariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    McGenRegistry.add("hPt", "Transverse momentum of muon pairs;;#counts", kTH1D, {axisPt});
    McGenRegistry.add("hEta", "Pseudorapidty of muon pairs;;#counts", kTH1D, {axisEta});
    McGenRegistry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});
    McGenRegistry.add("hPhi", "#varphi of muon pairs;;#counts", kTH1D, {axisPhi});
    McGenRegistry.add("hPtTrkPos", "Pt of positive muons;;#counts", kTH1D, {axisPtSingle});
    McGenRegistry.add("hPtTrkNeg", "Pt of negative muons;;#counts", kTH1D, {axisPtSingle});
    McGenRegistry.add("hEtaTrkPos", "#eta of positive muons;;#counts", kTH1D, {axisEtaSingle});
    McGenRegistry.add("hEtaTrkNeg", "#eta of negative muons;;#counts", kTH1D, {axisEtaSingle});
    McGenRegistry.add("hPhiTrkPos", "#varphi of positive muons;;#counts", kTH1D, {axisPhiSingle});
    McGenRegistry.add("hPhiTrkNeg", "#varphi of negative muons;;#counts", kTH1D, {axisPhiSingle});
    McGenRegistry.add("hPhiCharge", "#phi #it{charge}", kTH1D,{axisPhi});
    McGenRegistry.add("hPhiAverage", "#phi #it{average}", kTH1D,{axisPhi});

    //corr gen-reco
    McCorrReg.add("hPtcorr", "gen pT vs reco pT", kTH2D, {axisPt, axisPt});
  }

  // FUNCTIONS

  // template function that fills a map with the collision id of each udcollision as key
  // and a vector with the tracks
  // map == (key, element) == (udCollisionId, vector of trks)
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // template function that fills a map with the collision id of each udmccollision as key
  // and a vector with the tracks
  // map == (key, element) == (udMcCollisionId, vector of mc particles)
  template <typename TTracks>
  void collectMcCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udMcCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  template <typename TTracks>
  void collectMcGenRecoCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      if(tr.has_mcParticle()){
        int32_t candId = tr.udCollisionId();
        if (candId < 0) {
          continue;
        }
        tracksPerCand[candId].push_back(tr.globalIndex());
        
        auto mcPart = tr.mcParticle();
        tracksPerCand[candId].push_back(mcPart.globalIndex());
      }
    }
  }

  // struct used to store the ZDC info in a map
  struct ZDCinfo {
    float timeA;
    float timeC;
    float enA;
    float enC;
    int32_t id;
  };

  // function that fills a map with the collision id of each udcollision as key
  // and a ZDCinfo struct with the ZDC information
  void collectCandZDCInfo(std::unordered_map<int32_t, ZDCinfo>& zdcPerCand, o2::aod::UDZdcsReduced& ZDCs)
  {

    for (auto& zdc : ZDCs) {
      int32_t candId = zdc.udCollisionId();
      if (candId < 0) {
        continue;
      }

      zdcPerCand[candId].timeA = zdc.timeZNA();
      zdcPerCand[candId].timeC = zdc.timeZNC();
      zdcPerCand[candId].enA = zdc.energyCommonZNA();
      zdcPerCand[candId].enC = zdc.energyCommonZNC();

      // take care of the infinity
      if (std::isinf(zdcPerCand[candId].timeA))
        zdcPerCand[candId].timeA = -999;
      if (std::isinf(zdcPerCand[candId].timeC))
        zdcPerCand[candId].timeC = -999;
      if (std::isinf(zdcPerCand[candId].enA))
        zdcPerCand[candId].enA = -999;
      if (std::isinf(zdcPerCand[candId].enC))
        zdcPerCand[candId].enC = -999;
    }
  }

  // function to select muon tracks
  bool isMuonSelected(const ForwardTracks::iterator& fwdTrack)
  {
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();
    TLorentzVector p;
    p.SetXYZM(fwdTrack.px(), fwdTrack.py(), fwdTrack.pz(), mMu);
    float eta = p.Eta();
    float pt = p.Pt();
    float pDcaMax = rAbs < kRAbsMid ? kPDca1 : kPDca2;

    if (eta < kEtaMin || eta > kEtaMax)
      return false;
    if (pt < kPtMin)
      return false;
    if (rAbs < kRAbsMin || rAbs > kRAbsMax)
      return false;
    if (pDca > pDcaMax)
      return false;
    return true;
  }

  //compute phi for azimuth anisotropy
  void computePhiAnis(TLorentzVector p1, TLorentzVector p2, int sign1, float &phiAverage, float &phiCharge){
    
    TLorentzVector tSum, tDiffAv, tDiffCh;
    tSum = p1 + p2;
    if(sign1>0){
      tDiffCh = p1 - p2;
      if(gRandom->Rndm()>0.5)
        tDiffAv = p1 - p2;
      else tDiffAv = p2 - p1;
    }else{
      tDiffCh = p2 - p1;
      if(gRandom->Rndm()>0.5)
        tDiffAv = p2 - p1;
      else tDiffAv = p1 - p2;
    }

    //average
    phiAverage = tSum.DeltaPhi(tDiffAv);
    //charge
    phiCharge = tSum.DeltaPhi(tDiffCh);
  }
    

  // function that processes the candidates:
  // it applies V0 selection, trk selection, and fills the histograms
  void processCand(CandidatesFwd::iterator const& cand,
                   const ForwardTracks::iterator& tr1, const ForwardTracks::iterator& tr2,
                   ZDCinfo& zdc)
  {
    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > 100.)
          return;
      }
    }

    // track selection
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;
    if (!isMuonSelected(tr1))
      return;
    if (!isMuonSelected(tr2))
      return;

    //select MCH-MID tracks only
    //if(tr1.trackType() != 3) 
    //  return;
    //if(tr2.trackType() != 3) 
    //  return;

    // MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != 2)
      return;

    //cut on pair kinematics
    //select mass
    if(p.M() < lowMass)
      return;
    if(p.M() > highMass)
      return;
    //select pt
    if(p.Pt() < lowPt)
      return;
    if(p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    // select opposite charge events only
    if (cand.netCharge() != 0) {
      registry.fill(HIST("hSameSign"), cand.numContrib());
      return;
    }

    //compute phi for azimuth anisotropy
    TLorentzVector tSum, tDiffAv, tDiffCh;
    tSum = p1 + p2;
    if(tr1.sign()>0){
      tDiffCh = p1 - p2;
      if(gRandom->Rndm()>0.5)
        tDiffAv = p1 - p2;
      else tDiffAv = p2 - p1;
    }else{
      tDiffCh = p2 - p1;
      if(gRandom->Rndm()>0.5)
        tDiffAv = p2 - p1;
      else tDiffAv = p1 - p2;
    }
    //average
    float phiAverage = tSum.DeltaPhi(tDiffAv);
    float phiCharge = tSum.DeltaPhi(tDiffCh);

    // zdc info
    if (TMath::Abs(zdc.timeA) < 10)
      registry.fill(HIST("hTimeZNA"), zdc.timeA);
    if (TMath::Abs(zdc.timeC) < 10)
      registry.fill(HIST("hTimeZNC"), zdc.timeC);
    registry.fill(HIST("hEnergyZN"), zdc.enA, zdc.enC);

    // divide the events in neutron classes
    bool neutron_A = false;
    bool neutron_C = false;
    int znClass = -1;

    if (TMath::Abs(zdc.timeA) < 2)
      neutron_A = true;
    if (TMath::Abs(zdc.timeC) < 2)
      neutron_C = true;

    if (std::isinf(zdc.timeC))
      neutron_C = false;
    if (std::isinf(zdc.timeA))
      neutron_A = false;

    // fill the histos in neutron classes and assign neutron class label
    // 0n0n
    if (neutron_C == false && neutron_A == false) {
      znClass = 1;
      reg0n0n.fill(HIST("hMass"), p.M());
      reg0n0n.fill(HIST("hPt"), p.Pt());
      reg0n0n.fill(HIST("hPtFit"), p.Pt());
      reg0n0n.fill(HIST("hEta"), p.Eta());
      reg0n0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutron_A ^ neutron_C) { // Xn0n + 0nXn
      if (neutron_A)// select opposite charge events only
    if (cand.netCharge() != 0) {
      registry.fill(HIST("hSameSign"), cand.numContrib());
      return;
    }
        znClass = 3;
      regXn0n.fill(HIST("hMass"), p.M());
      regXn0n.fill(HIST("hPt"), p.Pt());
      regXn0n.fill(HIST("hPtFit"), p.Pt());
      regXn0n.fill(HIST("hEta"), p.Eta());
      regXn0n.fill(HIST("hRapidity"), p.Rapidity());
    } else if (neutron_A && neutron_C) { // XnXn
      znClass = 4;
      regXnXn.fill(HIST("hMass"), p.M());
      regXnXn.fill(HIST("hPt"), p.Pt());
      regXnXn.fill(HIST("hPtFit"), p.Pt());
      regXnXn.fill(HIST("hEta"), p.Eta());
      regXnXn.fill(HIST("hRapidity"), p.Rapidity());
    }

    // fill the histos without looking at neutron emission
    registry.fill(HIST("hContrib"), cand.numContrib());
    registry.fill(HIST("hPtTrkPos"), p1.Pt());
    registry.fill(HIST("hPtTrkNeg"), p2.Pt());
    registry.fill(HIST("hEtaTrkPos"), p1.Eta());
    registry.fill(HIST("hEtaTrkNeg"), p2.Eta());
    registry.fill(HIST("hPhiTrkPos"), p1.Phi());
    registry.fill(HIST("hPhiTrkNeg"), p2.Phi());
    registry.fill(HIST("hEvSign"), cand.netCharge());
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"), p.Pt());
    registry.fill(HIST("hPtFit"), p.Pt());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRapidity"), p.Rapidity());
    registry.fill(HIST("hPhi"), p.Phi());
    registry.fill(HIST("hCharge"), tr1.sign());
    registry.fill(HIST("hCharge"), tr2.sign());
    registry.fill(HIST("hPhiAverage"), phiAverage);
    registry.fill(HIST("hPhiCharge"), phiCharge);

    // store the event to save it into a tree
    if (tr1.sign() > 0) {
      dimuSel(p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.PseudoRapidity(), p1.Phi(),
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.PseudoRapidity(), p2.Phi(),
              zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, znClass);TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;
    } else {
      dimuSel(p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              p2.E(), p2.Px(), p2.Py(), p2.Pz(), p2.Pt(), p2.PseudoRapidity(), p2.Phi(),
              p1.E(), p1.Px(), p1.Py(), p1.Pz(), p1.Pt(), p1.PseudoRapidity(), p1.Phi(),
              zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, znClass);
    }
    
  }

  // function that processes the MC gen candidates:
  // it applies some kinematics cut and fills the histograms
  void processMcCand(aod::UDMcCollisions::iterator const& mcCand,
                     aod::UDMcParticles::iterator const& McPart1, aod::UDMcParticles::iterator const& McPart2){

    //check that all pairs are mu+mu-
    if(McPart1.pdgCode() + McPart2.pdgCode() != 0) LOGF(info,"PDG codes: %d | %d" ,McPart1.pdgCode(),McPart2.pdgCode());
    
    // create Lorentz vectors
    TLorentzVector p1, p2;
    p1.SetXYZM(McPart1.px(), McPart1.py(), McPart1.pz(), mMu);
    p2.SetXYZM(McPart2.px(), McPart2.py(), McPart2.pz(), mMu);
    TLorentzVector p = p1 + p2;

    //cut on pair kinematics
    //select mass
    if(p.M() < lowMass)
      return;
    if(p.M() > highMass)
      return;
    //select pt
    if(p.Pt() < lowPt)
      return;
    if(p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    //compute phi for azimuth anisotropy
    float phiAverage = 0;
    float phiCharge = 0;
    computePhiAnis(p1,p2,McPart1.pdgCode(),phiAverage, phiCharge);    

    // fill the histos
    McGenRegistry.fill(HIST("hPtTrkPos"), p1.Pt());
    McGenRegistry.fill(HIST("hPtTrkNeg"), p2.Pt());
    McGenRegistry.fill(HIST("hEtaTrkPos"), p1.Eta());
    McGenRegistry.fill(HIST("hEtaTrkNeg"), p2.Eta());
    McGenRegistry.fill(HIST("hPhiTrkPos"), p1.Phi());
    McGenRegistry.fill(HIST("hPhiTrkNeg"), p2.Phi());
    McGenRegistry.fill(HIST("hMass"), p.M());
    McGenRegistry.fill(HIST("hPt"), p.Pt());
    McGenRegistry.fill(HIST("hEta"), p.Eta());
    McGenRegistry.fill(HIST("hRapidity"), p.Rapidity());
    McGenRegistry.fill(HIST("hPhi"), p.Phi());
    McGenRegistry.fill(HIST("hPhiAverage"), phiAverage);
    McGenRegistry.fill(HIST("hPhiCharge"), phiCharge);

  }

  void processMcCandGenReco(CandidatesFwd::iterator const& cand,
                            const ForwardTracks::iterator& tr1, const ForwardTracks::iterator& tr2,
                            aod::UDMcParticles::iterator const& McPart1, aod::UDMcParticles::iterator const& McPart2)
  {
    // V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > 100.)
          return;
      }
    }

    // track selection
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;
    if (!isMuonSelected(tr1))
      return;
    if (!isMuonSelected(tr2))
      return;

    // MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != 2)
      return;

    // select opposite charge events only
    if (cand.netCharge() != 0) {
      registry.fill(HIST("hSameSign"), cand.numContrib());
      return;
    }

    // mc particle
    TLorentzVector p1Mc, p2Mc;
    p1Mc.SetXYZM(McPart1.px(), McPart1.py(), McPart1.pz(), mMu);
    p2Mc.SetXYZM(McPart2.px(), McPart2.py(), McPart2.pz(), mMu);
    TLorentzVector pMc = p1Mc + p2Mc;

    //cut on pair kinematics
    //select mass
    if(p.M() < lowMass)
      return;
    if(p.M() > highMass)
      return;
    //select pt
    if(p.Pt() < lowPt)
      return;
    if(p.Pt() > highPt)
      return;
    // select rapidity
    if (p.Rapidity() < lowRapidity)
      return;
    if (p.Rapidity() > highRapidity)
      return;

    McCorrReg.fill(HIST("hPtcorr"), p.Pt(),pMc.Pt());
  }

  // PROCESS FUNCTION
  void processDataAndReco(CandidatesFwd const& eventCandidates,
               o2::aod::UDZdcsReduced& ZDCs,
               ForwardTracks const& fwdTracks)
  {

    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // map with the ZDC info
    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    collectCandZDCInfo(zdcPerCand, ZDCs);

    // loop over the candidates
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      auto tr1 = fwdTracks.iteratorAt(trId1);
      auto tr2 = fwdTracks.iteratorAt(trId2);

      ZDCinfo zdc;

      if (zdcPerCand.count(candID) != 0) {
        zdc = zdcPerCand.at(candID);
      } else {
        zdc.timeA = -999;
        zdc.timeC = -999;
        zdc.enA = -999;
        zdc.enC = -999;
      }

      processCand(cand, tr1, tr2, zdc);
    }
  }

  PROCESS_SWITCH(fwdMuonsUPC, processDataAndReco, "", true);


  // process MC Truth
  void processGen(aod::UDMcCollisions const& mccollisions, aod::UDMcParticles const& McParts){
    
    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectMcCandIDs(tracksPerCand, McParts);

    // loop over the candidates
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = mccollisions.iteratorAt(candID);
      auto tr1 = McParts.iteratorAt(trId1);
      auto tr2 = McParts.iteratorAt(trId2);

      processMcCand(cand, tr1, tr2);
    }
  }
  PROCESS_SWITCH(fwdMuonsUPC, processGen, "", false);

/*
  void processGenAndReco(CandidatesFwd const& eventCandidates,
                       ForwardTracks const& fwdTracks,
                       aod::UDMcParticles const&)
  {
    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectMcGenRecoCandIDs(tracksPerCand, fwdTracks);

    // loop over the candidates
      for (const auto& item : tracksPerCand) {
        int32_t trId1 = item.second[0];
        int32_t trId2 = item.second[1];

        int32_t candID = item.first;
        auto cand = eventCandidates.iteratorAt(candID);
        auto tr1 = fwdTracks.iteratorAt(trId1);
        auto tr2 = fwdTracks.iteratorAt(trId2);

        int32_t trMcId1 = item.second[2];
        int32_t trMcId2 = item.second[3];
        auto trMc1 = fwdTracks.iteratorAt(trMcId1);
        auto trMc2 = fwdTracks.iteratorAt(trMcId2);

        processMcCandGenReco(cand, tr1, tr2, trMc1, trMc2);
      }
  }
  PROCESS_SWITCH(fwdMuonsUPC, processGenAndReco, "", false);
*/
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<fwdMuonsUPC>(cfgc),
  };
}
