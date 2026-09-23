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

/// \file upcTrackVertexingQA.cxx
/// \brief task to study the performance of vertexing for low-multiplicity UPC collisions
/// \author Andrea Tavira Garcia a.tavira@cern.ch
/// \author Andrea Giovanni Riffero andrea.giovanni.riffero@cern.ch

///

#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>
#include <TH2.h>

#include <array>
#include <bit>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcTrackVertexingQA {
  // Configurables
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 0., "min. cand. pT (GeV/c)"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track pT (GeV/c)"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.9, "max. |eta| of tracks"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 3.f, "max. TPC N_sigma (pion)"};
  Configurable<int>   nMinTpcClusters{"nMinTpcClusters", 60, "min. number of TPC clusters"};
  Configurable<float> massMin{"massMin", 0.5, "min. inv. mass (GeV/c^2)"};
  Configurable<float> massMax{"massMax", 1.3, "max. inv. mass (GeV/c^2)"};

  // Name shortenings
  // passed* columns are in TrackSelectionExtension; isGlobalTrack* and trackCutFlag in TrackSelection
  using TracksExtraWPidPi = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                      aod::TrackSelection, aod::TrackSelectionExtension,
                                      aod::pidTPCFullPi>;

  static constexpr float MassPion = o2::constants::physics::MassPionCharged;

  // Ordered list of cuts. The cumulative cut flow follows this order,
  // so reorder CutLabels and getCutResults together to change what "tighter" means.
  static constexpr int NCuts = 8;
  static constexpr std::array<const char*, NCuts> CutLabels = {
    "all tracks", "hasTPC", "passedEtaRange",
    "passedTPCNCls", "passedTPCChi2NDF", "passedITSNCls",
    "passedITSChi2NDF", "hasITS"};

  // Axes
  ConfigurableAxis axisMass{"axisMass", {160, 0.4, 1.2}, "m_{#pi#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPt{"axisPt", {200, 0., 2.}, "#it{p}_{T} (GeV/#it{c})"};

  HistogramRegistry registry{
    "registry",
    {// candidate level
      {"Cand/hMass", ";m_{#pi#pi} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}}},
      {"Cand/hPt", ";#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}}},
      {"Cand/hRapidity", ";#it{y};entries", {HistType::kTH1F, {{100, -1., 1.}}}},

      // collision level
      {"Coll/hNContrib", ";N_{PV contributors};entries", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
      {"Coll/hBCid", ";BCid;entries", {HistType::kTH1F, {{1000, 0., 1000.}}}},
      {"Coll/hVtxZ", ";#it{z}_{vtx} (cm);entries", {HistType::kTH1F, {{200, -20., 20.}}}},
      {"Coll/hVtxX", ";#it{x}_{vtx} (cm);entries", {HistType::kTH1F, {{200, -0.05, 0.05}}}},
      {"Coll/hVtxY", ";#it{y}_{vtx} (cm);entries", {HistType::kTH1F, {{200, -0.05, 0.05}}}},
      {"Coll/hVtxChi2", ";#chi^{2} vtx;entries", {HistType::kTH1F, {{100, 0., 10.}}}},

      // track level (prongs of the candidate)
      {"Trk/hPt", ";#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}}},
      {"Trk/hEta", ";#eta;entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"Trk/hChi2NCl", ";#chi^{2}/N_{cls} TPC;entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"Trk/hTpcSignalVsP", ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x}", {HistType::kTH2F, {{200, 0., 2.}, {300, 0., 300.}}}},
      {"Trk/hNSigmaPiVsP", ";#it{p} (GeV/#it{c});n#sigma^{TPC}_{#pi}", {HistType::kTH2F, {{200, 0., 2.}, {100, -10., 10.}}}},
      {"Trk/hHasIts", ";has ITS;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"Trk/hIsPvContrib", ";is PV contributor;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"Trk/hTpcNClsFound", ";N_{cls} TPC;entries", {HistType::kTH1F, {{160, 0., 160.}}}},
      {"Trk/hItsChi2NCl", ";#chi^{2}/N_{cls} ITS;entries", {HistType::kTH1F, {{100, 0., 40.}}}},
      {"Trk/hItsNCls", ";N_{cls} ITS;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"Trk/hItsNClsInnerBarrel", ";N_{cls} ITS Inner Barrel;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"Trk/hDcaXY", ";DCA_{xy} (cm);entries", {HistType::kTH1F, {{200, -0.1, 0.1}}}},
      {"Trk/hDcaZ", ";DCA_{z} (cm);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}}}
    }
  };
     

  void init(InitContext&)
  {
    // Cut flows: cumulative, single cut, and tracks per collision vs. cut step
    registry.add("Cut/hCutFlowCumulative", "tracks surviving cuts applied in sequence;;entries",
                 HistType::kTH1F, {{NCuts, -0.5, NCuts - 0.5}});
    registry.add("Cut/hCutFlowSingle", "tracks passing each cut individually;;entries",
                 HistType::kTH1F, {{NCuts, -0.5, NCuts - 0.5}});
    registry.add("Cut/h2NTracksPerCollVsCut", "tracks per collision after each cumulative cut;;N_{tracks} / collision",
                 HistType::kTH2F, {{NCuts, -0.5, NCuts - 0.5}, {10, -0.5, 9.5}});

    auto hCum = registry.get<TH1>(HIST("Cut/hCutFlowCumulative"));
    auto hSingle = registry.get<TH1>(HIST("Cut/hCutFlowSingle"));
    auto hPerColl = registry.get<TH2>(HIST("Cut/h2NTracksPerCollVsCut"));
    for (int i = 0; i < NCuts; ++i) {
      hCum->GetXaxis()->SetBinLabel(i + 1, CutLabels[i]);
      hSingle->GetXaxis()->SetBinLabel(i + 1, CutLabels[i]);
      hPerColl->GetXaxis()->SetBinLabel(i + 1, CutLabels[i]);
    }
  }

  // Result of every cut for one track (same order as CutLabels)
  template <typename TTrack>
  std::array<bool, NCuts> getCutResults(TTrack const& track)
  {
    return {true,
            track.hasTPC(),
            track.passedTPCNCls(),
            track.passedEtaRange(),
            track.passedTPCChi2NDF(),
            track.passedITSNCls(),
            track.passedITSChi2NDF(),
            track.hasITS()};
  }

  // Fills the cut flows for all tracks of one collision
  template <typename TTracks>
  void fillCutFlow(TTracks const& tracks)
  {
    std::array<int, NCuts> nSurvivors{}; // tracks of this collision surviving up to step i

    for (auto const& track : tracks) {
      auto passed = getCutResults(track);
      bool passedSoFar = true;
      for (int i = 0; i < NCuts; ++i) {
        if (passed[i]) {
          registry.fill(HIST("Cut/hCutFlowSingle"), i);
        }
        passedSoFar = passedSoFar && passed[i];
        if (passedSoFar) {
          registry.fill(HIST("Cut/hCutFlowCumulative"), i);
          ++nSurvivors[i];
        }
      }
    }

    for (int i = 0; i < NCuts; ++i) {
      registry.fill(HIST("Cut/h2NTracksPerCollVsCut"), i, nSurvivors[i]);
    }
  }

  // Basic kinematic / TPC properties of a candidate prong
  template <typename TTrack>
  void checkTpcTrackProperties(TTrack const& track)
  {
    registry.fill(HIST("Trk/hPt"), track.pt());
    registry.fill(HIST("Trk/hEta"), track.eta());
    registry.fill(HIST("Trk/hChi2NCl"), track.tpcChi2NCl());
    registry.fill(HIST("Trk/hTpcSignalVsP"), track.p(), track.tpcSignal());
    registry.fill(HIST("Trk/hNSigmaPiVsP"), track.p(), track.tpcNSigmaPi());
    registry.fill(HIST("Trk/hHasIts"), static_cast<int>(track.hasITS()));
    registry.fill(HIST("Trk/hIsPvContrib"), static_cast<int>(track.isPVContributor()));
    registry.fill(HIST("Trk/hTpcNClsFound"), track.tpcNClsFound());
    registry.fill(HIST("Trk/hItsChi2NCl"), track.itsChi2NCl());
    registry.fill(HIST("Trk/hItsNCls"), track.itsNCls());
    registry.fill(HIST("Trk/hItsNClsInnerBarrel"), track.itsNClsInnerBarrel());
    registry.fill(HIST("Trk/hDcaXY"), track.dcaXY());
    registry.fill(HIST("Trk/hDcaZ"), track.dcaZ());

    // TODO: ambiguity (aod::AmbiguousTracks), ...
  }

  // Basic single-track selection used to build the rho candidate
  template <typename TTrack>
  bool isGoodTrack(TTrack const& track)
  {
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsFound() < nMinTpcClusters) {
      return false;
    }
    if (track.pt() < ptTrackMin || std::abs(track.eta()) > etaTrackMax) {
      return false;
    }
    return std::abs(track.tpcNSigmaPi()) <= nSigmaTpcMax;
  }

  // Tracks are grouped by collision automatically
  void processRhoCand(aod::Collision const& collision, TracksExtraWPidPi const& tracks)
  {
    registry.fill(HIST("Coll/hNContrib"), collision.numContrib());
    registry.fill(HIST("Coll/hBCid"), collision.bcId());
    registry.fill(HIST("Coll/hVtxZ"), collision.posZ());
    registry.fill(HIST("Coll/hVtxX"), collision.posX());
    registry.fill(HIST("Coll/hVtxY"), collision.posY());
    registry.fill(HIST("Coll/hVtxChi2"), collision.chi2());

    // Sequential ITS/TPC cuts on ALL tracks of the collision
    fillCutFlow(tracks);

    // Select tracks for the rho candidate
    std::vector<decltype(tracks.begin())> goodTracks;
    for (auto const& track : tracks) {
      if (isGoodTrack(track)) {
        goodTracks.push_back(track);
      }
    }

    // Exactly two tracks with opposite charge -> rho candidate
    if (goodTracks.size() != 2) {
      return;
    }
    auto const& track0 = goodTracks[0];
    auto const& track1 = goodTracks[1];
    if (track0.sign() * track1.sign() >= 0) {
      return;
    }

    ROOT::Math::PxPyPzMVector p0(track0.px(), track0.py(), track0.pz(), MassPion);
    ROOT::Math::PxPyPzMVector p1(track1.px(), track1.py(), track1.pz(), MassPion);
    auto rho = p0 + p1;

    if (rho.M() < massMin || rho.M() > massMax) {
      return;
    }
    if (rho.Pt() < ptCandMin || std::abs(rho.Rapidity()) > yCandMax) {
      return;
    }

    registry.fill(HIST("Cand/hMass"), rho.M());
    registry.fill(HIST("Cand/hPt"), rho.Pt());
    registry.fill(HIST("Cand/hRapidity"), rho.Rapidity());

    checkTpcTrackProperties(track0);
    checkTpcTrackProperties(track1);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoCand, "Rho -> pi pi candidates and track QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcTrackVertexingQA>(cfgc)};
}
