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

#include "TLorentzVector.h"
#include "TSystem.h"


namespace mpt {
  DECLARE_SOA_COLUMN(M, m, float);
  DECLARE_SOA_COLUMN(Pt, pt, float);
}

namespace o2::aod {
  DECLARE_SOA_TABLE(MassPt, "AOD", "MASSPT",
		    mpt::M, mpt::Pt);
}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//defining constants
double mMu = 0.10566; //mass of muon

//constants used in the track selection
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
  using udZDCreduced = o2::aod::UDZdcsReduced::iterator;

  Produces<o2::aod::MassPt> mptSel;
  
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //CONFIGURABLES
  //pT of muon pairs
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0., "lower limit in pT histo"};
  Configurable<float> highPt{"highPt", 0.5, "upper limit in pT histo"};
  //mass of muon pairs
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in mass histo"};
  Configurable<float> lowMass{"lowMass", 0., "lower limit in mass histo"};
  Configurable<float> highMass{"highMass", 10., "upper limit in mass histo"};
  //eta of muon pairs
  Configurable<int> nBinsEta{"nBinsEta", 600, "N bins in eta histo"};
  Configurable<float> lowEta{"lowEta", -10., "lower limit in eta histo"};
  Configurable<float> highEta{"highEta", -2., "upper limit in eta histo"};
  //rapidity of muon pairs
  Configurable<int> nBinsRapidity{"nBinsRapidity", 250, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4., "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", -2.5, "upper limit in rapidity histo"};
  //phi of muon pairs
  Configurable<int> nBinsPhi{"nBinsPhi", 600, "N bins in phi histo"};
  Configurable<float> lowPhi{"lowPhi", -TMath::Pi(), "lower limit in phi histo"};
  Configurable<float> highPhi{"highPhi", TMath::Pi(), "upper limit in phi histo"};
  //pT of single muons
  Configurable<int> nBinsPtSingle{"nBinsPtSingle", 500, "N bins in pT histo single muon"};
  Configurable<float> lowPtSingle{"lowPtSingle", 0., "lower limit in pT histo single muon"};
  Configurable<float> highPtSingle{"highPtSingle", 2., "upper limit in pT histo single muon"};
  //eta of single muons
  Configurable<int> nBinsEtaSingle{"nBinsEtaSingle", 250, "N bins in eta histo single muon"};
  Configurable<float> lowEtaSingle{"lowEtaSingle", -4.5, "lower limit in eta histo single muon"};
  Configurable<float> highEtaSingle{"highEtaSingle", -2., "upper limit in eta histo single muon"};
  //phi of single muons
  Configurable<int> nBinsPhiSingle{"nBinsPhiSingle", 600, "N bins in phi histo single muon"};
  Configurable<float> lowPhiSingle{"lowPhiSingle", -TMath::Pi(), "lower limit in phi histo single muon"};
  Configurable<float> highPhiSingle{"highPhiSingle", TMath::Pi(), "upper limit in phi histo single muon"};
  //ZDC
  Configurable<int> nBinsZDCen{"nBinsZDCen",200,"N bins in ZN energy"};
  Configurable<float> lowEnZN{"lowEnZN",-50.,"lower limit in ZN energy histo"};
  Configurable<float> highEnZN{"highEnZN",250.,"upper limit in ZN energy histo"};


  void init(InitContext&)
  {
    //axis
    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisMass{nBinsMass, lowMass, highMass, "m_{#mu#mu} GeV/#it{c}^{2}"};
    const AxisSpec axisEta{nBinsEta, lowEta, highEta, "#eta"};
    const AxisSpec axisRapidity{nBinsRapidity, lowRapidity, highRapidity, "Rapidity"};
    const AxisSpec axisPhi{nBinsPhi, lowPhi, highPhi, "#varphi"};
    const AxisSpec axisPtSingle{nBinsPtSingle, lowPtSingle, highPtSingle, "#it{p}_{T}_{ trk} GeV/#it{c}"};
    const AxisSpec axisTimeZN{500,-5,5,"ZDC time (ns)"};
    const AxisSpec axisEnergyZNA{nBinsZDCen,lowEnZN,highEnZN,"ZNA energy (TeV)"};
    const AxisSpec axisEnergyZNC{nBinsZDCen,lowEnZN,highEnZN,"ZNC energy (TeV)"};
    const AxisSpec axisEtaSingle{nBinsEtaSingle, lowEtaSingle, highEtaSingle, "#eta_{trk}"};
    const AxisSpec axisPhiSingle{nBinsPhiSingle, lowPhiSingle, highPhiSingle, "#varphi_{trk}"};
    
    //histos
    registry.add("hMass", "Ivariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPt", "Transverse momentum mass of muon pairs;;#counts", kTH1D, {axisPt});
    registry.add("hEta", "Pseudorapidty of muon pairs;;#counts", kTH1D, {axisEta});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});
    registry.add("hPhi", "#varphi of muon pairs;;#counts", kTH1D, {axisPhi});
    registry.add("hCharge", "Charge;#it{charge};;#counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("hContrib", "hContrib;;#counts",kTH1D,{{6,-0.5,5.5}});
    registry.add("hEvSign","Sum of the charges of all the tracks in each event;;#counts",kTH1D,{{5,-2.5,2.5}});
    registry.add("hPtTrkPos", "Pt of positive muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hPtTrkNeg", "Pt of negative muons;;#counts", kTH1D, {axisPtSingle});
    registry.add("hEtaTrkPos", "#eta of positive muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hEtaTrkNeg", "#eta of negative muons;;#counts", kTH1D, {axisEtaSingle});
    registry.add("hPhiTrkPos", "#varphi of positive muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hPhiTrkNeg", "#varphi of negative muons;;#counts", kTH1D, {axisPhiSingle});
    registry.add("hTimeZNA","ZNA Times;;#counts",kTH1D,{axisTimeZN});
    registry.add("hTimeZNC","ZNC Times;;#counts",kTH1D,{axisTimeZN});
    registry.add("hEnergyZN","ZNA vs ZNC energy",kTH2D,{axisEnergyZNA,axisEnergyZNC});
  }


  //FUNCTIONS

  //template function that fills a map with the collision id of each udcollision as key
  //and a vector with the tracks
  //map == (key, element) == (udCollisionId, vector of trks)
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

  // struct used to store the ZDC info in a map
  struct ZDCinfo{
    float timeA;
    float timeC;
    float enA;
    float enC;
    int32_t id;
  };

  //function that fills a map with the collision id of each udcollision as key
  //and a ZDCinfo struct with the ZDC information 
  void collectCandZDCInfo(std::unordered_map<int32_t, ZDCinfo>& zdcPerCand, o2::aod::UDZdcsReduced& ZDCs){
    
    for(auto &zdc: ZDCs){
      int32_t candId = zdc.udCollisionId();
      if (candId < 0) {
        continue;
      }
      zdcPerCand[candId].timeA = zdc.timeZNA();
      zdcPerCand[candId].timeC = zdc.timeZNC();
      zdcPerCand[candId].enA = zdc.energyCommonZNA();
      zdcPerCand[candId].enC = zdc.energyCommonZNC();

      //LOGP(info, "zdc = {}", zdc.timeZNA());
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


  //function that processes the candidates:
  //it applies V0 selection, trk selection, and fills the histograms
  void processCand(CandidatesFwd::iterator const& cand,
                   const ForwardTracks::iterator& tr1, const ForwardTracks::iterator& tr2,
                   ZDCinfo& zdc){
    //V0 selection
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1){
        if(ampsV0A[i] > 100.) return;
      }
    }

    //select opposite charge events only
    if(cand.netCharge()!=0) return;

    //track selection
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;
    if (!isMuonSelected(tr1))
      return;
    if (!isMuonSelected(tr2))
      return;

    //MCH-MID match selection
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != 2)
      return;
    
    //fill the histos
    registry.fill(HIST("hContrib"), cand.numContrib());
    registry.fill(HIST("hPtTrkPos"), p1.Pt());
    registry.fill(HIST("hPtTrkNeg"), p2.Pt());
    registry.fill(HIST("hEtaTrkPos"), p1.Eta());
    registry.fill(HIST("hEtaTrkNeg"), p2.Eta());
    registry.fill(HIST("hPhiTrkPos"), p1.Phi());
    registry.fill(HIST("hPhiTrkNeg"), p2.Phi());
    registry.fill(HIST("hEvSign"),cand.netCharge());
    registry.fill(HIST("hMass"), p.M());
    registry.fill(HIST("hPt"),  p.Pt());
    registry.fill(HIST("hEta"), p.Eta());
    registry.fill(HIST("hRapidity"), p.Rapidity());
    registry.fill(HIST("hPhi"), p.Phi());
    registry.fill(HIST("hCharge"),tr1.sign());
    registry.fill(HIST("hCharge"),tr2.sign());

    //zdc info
    registry.fill(HIST("hTimeZNA"), zdc.timeA);
    registry.fill(HIST("hTimeZNC"), zdc.timeC);

    registry.fill(HIST("hEnergyZN"),zdc.enA,zdc.enC);

    // restric the kine
    if (p.M() > 2 && p.M() < 6 && p.Pt() < 5) mptSel(p.M(),p.Pt());
  }


  //PROCESS FUNCTION
  void process(CandidatesFwd const& eventCandidates,
               o2::aod::UDZdcsReduced& ZDCs,
               ForwardTracks const& fwdTracks)
  {

    // map with the tracks
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    //map with the ZDC info
    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    collectCandZDCInfo(zdcPerCand, ZDCs);

    int match = 0;
    int tot = 0;
    for (const auto& item : tracksPerCand) {
      
      int32_t candID = item.first;
      match = match + zdcPerCand.count(candID);
      tot++;
      
    }
    //LOGP(info,"matches = {}",match);
    //LOGP(info,"tot = {}",tot);

    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      auto tr1 = fwdTracks.iteratorAt(trId1);
      auto tr2 = fwdTracks.iteratorAt(trId2);
      if(zdcPerCand.count(candID)==0) continue;
      auto zdc = zdcPerCand.at(candID);
      
      processCand(cand, tr1, tr2, zdc);
      
    }

    

  }

  PROCESS_SWITCH(fwdMuonsUPC, process, "", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<fwdMuonsUPC>(cfgc)};
}
