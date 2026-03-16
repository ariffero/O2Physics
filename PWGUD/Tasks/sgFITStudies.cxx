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
//
// \Single Gap Event Studies
// \author Andrea Giovanni Riffero andrea.giovanni.riffero@cern.ch
// \since  March 2026

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/PIDResponseTOF.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


struct SGFITStudies {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  // ConfigurableAxis BCAxis{"BCAxis", {1000000000000, 0.5, 1000000000000.5}, "BC axis"};
  ConfigurableAxis BCAxis{"BCAxis", {100000000000, 500000000000.5, 600000000000.5}, "BC axis"};
  ConfigurableAxis etaAxis{"etaAxis", {300, -1.5, 1.5}, ""};
  ConfigurableAxis sigTPCAxis{"sigTPCAxis", {100, -100.0, 100.0}, ""};
  ConfigurableAxis sigTOFAxis{"sigTOFAxis", {100, -100.0, 100.0}, ""};
  ConfigurableAxis multAxis{"multAxis", {51, -.5, 50.5}, ""};
  ConfigurableAxis FitAxis{"FitAxis", {2000, -0.5, 3999.5}, ""};
  ConfigurableAxis ZDCAxis{"ZDCAxis", {1000, -2.5, 199.5}, ""};
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<std::string> outputFileName{"outputFileName", "AnalysisResults.root", "Output file name"};
  
  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}
  }; 

  void init(InitContext&){
    const AxisSpec axispt{ptAxis, "p_{T}"};
    const AxisSpec axismeanpt{ptAxis, "<p_{T}>"};
    const AxisSpec axisBC{BCAxis, "BC"};
    const AxisSpec axiseta{etaAxis, "#eta"};
    const AxisSpec axismult{multAxis, "N_{tracks}"};
    const AxisSpec axisfit{FitAxis, "FIT Amplitude"};
    const AxisSpec axiszdc{ZDCAxis, "ZDC Amplitude"};
    
    // collision histograms
    registry.add("collisions/GapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("collisions/TrueGapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{4, -1.5, 2.5}}});

    // track histograms
    registry.add("tracks/QCAll", "Track QC of all tracks; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/QCPVC", "Track QC of PV contributors; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/etaA", "track eta of PV contributors - A gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaC", "track eta of PV contributors - C gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaAC", "track eta of PV contributors - AC gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaApv", "track eta of PV contributors - A gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaCpv", "track eta of PV contributors - C gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaACpv", "track eta of PV contributors - AC gap; #eta; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etavsptA", "track eta versus pt of all tracks; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});
    registry.add("tracks/etavsptC", "track eta versus pt of PV contributors; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});
    registry.add("tracks/etavsptAC", "track eta versus pt of PV contributors; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});

    // FIT histograms
    // A gap
    registry.add("FIT/AFT0A", "Amplitude FT0A - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C", "Amplitude FT0C - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDA", "Amplitude FDDA - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC", "Amplitude FDDC - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A", "Amplitude FV0A - A gap", {HistType::kTH1F, {{axisfit}}});
    // C gap
    registry.add("FIT/CFT0A", "Amplitude FT0A - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0C", "Amplitude FT0C - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA", "Amplitude FDDA - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDC", "Amplitude FDDC - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A", "Amplitude FV0A - C gap", {HistType::kTH1F, {{axisfit}}});

    // FIT histos +
    // A gap, > 2 tracks
    registry.add("FIT/ApnFT0A", "Amplitude FT0A - > 2 tracks - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ApnFT0C", "Amplitude FT0C - > 2 tracks - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ApnFDDA", "Amplitude FDDA - > 2 tracks - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ApnFDDC", "Amplitude FDDC - > 2 tracks - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ApnFV0A", "Amplitude FV0A - > 2 tracks - A gap", {HistType::kTH1F, {{axisfit}}});
    // C gap, > 2 tracks
    registry.add("FIT/CpnFT0A", "Amplitude FT0A - > 2 tracks - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CpnFT0C", "Amplitude FT0C - > 2 tracks - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CpnFDDA", "Amplitude FDDA - > 2 tracks - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CpnFDDC", "Amplitude FDDC - > 2 tracks - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CpnFV0A", "Amplitude FV0A - > 2 tracks - C gap", {HistType::kTH1F, {{axisfit}}});
    // A gap, 2 tracks, coh
    registry.add("FIT/A2cFT0A", "Amplitude FT0A - 2 tracks coh - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2cFT0C", "Amplitude FT0C - 2 tracks coh - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2cFDDA", "Amplitude FDDA - 2 tracks coh - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2cFDDC", "Amplitude FDDC - 2 tracks coh - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2cFV0A", "Amplitude FV0A - 2 tracks coh - A gap", {HistType::kTH1F, {{axisfit}}});
    // C gap, 2 tracks, coh
    registry.add("FIT/C2cFT0A", "Amplitude FT0A - 2 tracks coh - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2cFT0C", "Amplitude FT0C - 2 tracks coh - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2cFDDA", "Amplitude FDDA - 2 tracks coh - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2cFDDC", "Amplitude FDDC - 2 tracks coh - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2cFV0A", "Amplitude FV0A - 2 tracks coh - C gap", {HistType::kTH1F, {{axisfit}}});
    // A gap, 2 tracks, incoh
    registry.add("FIT/A2iFT0A", "Amplitude FT0A - 2 tracks inc - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2iFT0C", "Amplitude FT0C - 2 tracks inc - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2iFDDA", "Amplitude FDDA - 2 tracks inc - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2iFDDC", "Amplitude FDDC - 2 tracks inc - A gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/A2iFV0A", "Amplitude FV0A - 2 tracks inc - A gap", {HistType::kTH1F, {{axisfit}}});
    // C gap, 2 tracks, incoh
    registry.add("FIT/C2iFT0A", "Amplitude FT0A - 2 tracks inc - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2iFT0C", "Amplitude FT0C - 2 tracks inc - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2iFDDA", "Amplitude FDDA - 2 tracks inc - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2iFDDC", "Amplitude FDDC - 2 tracks inc - C gap", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/C2iFV0A", "Amplitude FV0A - 2 tracks inc - C gap", {HistType::kTH1F, {{axisfit}}});

  // gap side histos cosmetics
  const char *gapSideLabels[3] = {"A", "C", "AC"};
  const char *trueGapSideLabels[4] = {"No Gap", "A", "C", "AC"};
  for (int i=0; i<3; i++) {
    registry.get<TH1>(HIST("collisions/GapSide"))->GetXaxis()->SetBinLabel(i+1, gapSideLabels[i]);
  }
  for (int i=0; i<4; i++) {
    registry.get<TH1>(HIST("collisions/TrueGapSide"))->GetXaxis()->SetBinLabel(i+1, trueGapSideLabels[i]);
  }

  }

  
  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    if (verbose) {
      LOGF(info, "DG candidate %d", dgcand.globalIndex());
    }

    // constants
    const float mpion = pdg->Mass(211);
    const float mmuon = pdg->Mass(13);
    
    // fill gap side histo
    registry.get<TH1>(HIST("collisions/GapSide"))->Fill(dgcand.gapSide(), 1.);

    //true gap side
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(dgcand, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    // fill true gap side histo
    registry.get<TH1>(HIST("collisions/TrueGapSide"))->Fill(truegapSide, 1.);

    // select PV contributors
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);

    // fill track histograms
    if (verbose) {
      LOGF(info, "Number of tracks %d", dgtracks.size());
      LOGF(info, "Number of PV contributors %d", PVContributors.size());
    }

    // parameters for track selection
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};

    std::vector<TLorentzVector> goodTracks;
    for(auto t : dgtracks) {
      if (trackselector(t, parameters)){
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), mpion);
        goodTracks.push_back(a);
      }
    } // end of loop on tracks

    // fill FIT histos
    if (truegapSide == o2::aod::sgselector::SingleGapA){
      registry.get<TH1>(HIST("FIT/AFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/AFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/AFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/AFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/AFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
    }
    else if(truegapSide == o2::aod::sgselector::SingleGapC){
      registry.get<TH1>(HIST("FIT/CFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/CFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/CFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/CFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/CFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
    }

    if(goodTracks.size() > 2){
      if (truegapSide == o2::aod::sgselector::SingleGapA){
        registry.get<TH1>(HIST("FIT/ApnFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ApnFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/ApnFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ApnFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/ApnFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
      }
      else if(truegapSide == o2::aod::sgselector::SingleGapC){
        registry.get<TH1>(HIST("FIT/CpnFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CpnFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/CpnFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CpnFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/CpnFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
      }
    }

    if(goodTracks.size() == 2){
      float pT = (goodTracks[0] + goodTracks[1]).Pt();
        if(pT < 0.2){
          if (truegapSide == o2::aod::sgselector::SingleGapA){
          registry.get<TH1>(HIST("FIT/A2cFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/A2cFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/A2cFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/A2cFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/A2cFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        }
        else if(truegapSide == o2::aod::sgselector::SingleGapC){
          registry.get<TH1>(HIST("FIT/C2cFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/C2cFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/C2cFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/C2cFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/C2cFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        }
      }
      else{
        if (truegapSide == o2::aod::sgselector::SingleGapA){
          registry.get<TH1>(HIST("FIT/A2iFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/A2iFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/A2iFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/A2iFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/A2iFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        }
        else if(truegapSide == o2::aod::sgselector::SingleGapC){
          registry.get<TH1>(HIST("FIT/C2iFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/C2iFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/C2iFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
          registry.get<TH1>(HIST("FIT/C2iFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
          registry.get<TH1>(HIST("FIT/C2iFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        }
      }
    }

  } // end of process function

};  //end of struct

  WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGFITStudies>(cfgc),
  };
}
