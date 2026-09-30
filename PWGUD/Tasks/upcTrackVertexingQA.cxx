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
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum class CandSpecies {
  kRho = 0,
  kJpsi = 1
};

struct UpcTrackVertexingQA {
  // Configurables
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 0., "min. cand. pT (GeV/c)"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track pT (GeV/c)"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.9, "max. |eta| of tracks"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 3.f, "max. TPC N_sigma (pion)"};
  Configurable<int> nMinTpcClusters{"nMinTpcClusters", 60, "min. number of TPC clusters"};
  Configurable<float> massMin{"massMin", 0.5, "min. inv. mass (GeV/c^2)"};
  Configurable<float> massMax{"massMax", 1.3, "max. inv. mass (GeV/c^2)"};

  // Name shortenings
  // passed* columns are in TrackSelectionExtension; isGlobalTrack* and trackCutFlag in TrackSelection
  using TracksExtraSels = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::TrackSelection, aod::TrackSelectionExtension>;
  using TracksExtraWPidPi = soa::Join<TracksExtraSels, aod::pidTPCFullPi>;
  using TracksExtraWPidMu = soa::Join<TracksExtraSels, aod::pidTPCFullMu>;
  using TracksExtraWPidPiMc = soa::Join<TracksExtraWPidPi, aod::McTrackLabels>;
  using TracksExtraWPidMuMc = soa::Join<TracksExtraWPidMu, aod::McTrackLabels>;
  using CollisionsWithMc = soa::Join<aod::Collisions, aod::McCollisionLabels>;

  static constexpr float MassPion = o2::constants::physics::MassPionCharged;
  static constexpr float MassMuon = o2::constants::physics::MassMuon;

  // Ordered list of cuts. The cumulative cut flow follows this order,
  // so reorder CutLabels and getCutResults together to change what "tighter" means.
  static constexpr int NCuts = 8;
  static constexpr std::array<const char*, NCuts> CutLabels = {
    "all tracks", "hasTPC", "passedTPCNCls",
    "passedEtaRange", "passedTPCChi2NDF", "passedITSNCls",
    "passedITSChi2NDF", "hasITS"};

  // Axes (the mass axis matches the default rho window; override it for the J/psi)
  ConfigurableAxis axisMass{"axisMass", {200, 0.5, 1.3}, "m_{#pi#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPt{"axisPt", {200, 0., 2.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisY{"axisY", {100, -1., 1.}, "Rapidity"};

  HistogramConfigSpec hVtxZ{HistType::kTH1F, {{200, -20., 20.}}};
  HistogramConfigSpec hVtxXY{HistType::kTH1F, {{200, -0.05, 0.05}}};
  HistogramConfigSpec hPtTrk{HistType::kTH1F, {{200, 0.0, 2.0}}};
  HistogramConfigSpec hEtaTrk{HistType::kTH1F, {{100, -1.0, 1.0}}};

  HistogramRegistry registry{
    "registry",
    {
      // collision level
      {"Coll/hNContrib", ";N_{PV contributors};entries", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
      {"Coll/hBCid", ";BCid;entries", {HistType::kTH1F, {{1000, 0., 1000.}}}},
      {"Coll/hVtxZ", ";#it{z}_{vtx} (cm);entries", hVtxZ},
      {"Coll/hVtxX", ";#it{x}_{vtx} (cm);entries", hVtxXY},
      {"Coll/hVtxY", ";#it{y}_{vtx} (cm);entries", hVtxXY},
      {"Coll/hVtxChi2", ";#chi^{2} vtx;entries", {HistType::kTH1F, {{100, 0., 10.}}}},

      // track level (all good tracks, before grouping by collision)
      {"Trk/hPt", ";#it{p}_{T} (GeV/#it{c});entries", hPtTrk},
      {"Trk/hEta", ";#eta;entries", hEtaTrk},
      {"Trk/hChi2NCl", ";#chi^{2}/N_{cls} TPC;entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"Trk/hTpcSignalVsP", ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x}", {HistType::kTH2F, {{200, 0., 2.}, {300, 0., 300.}}}},
      {"Trk/hNSigmaVsP", ";#it{p} (GeV/#it{c});n#sigma^{TPC}", {HistType::kTH2F, {{200, 0., 2.}, {100, -10., 10.}}}},
      {"Trk/hHasIts", ";has ITS;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"Trk/hIsPvContrib", ";is PV contributor;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"Trk/hTpcNClsFound", ";N_{cls} TPC;entries", {HistType::kTH1F, {{160, 0., 160.}}}},
      {"Trk/hItsChi2NCl", ";#chi^{2}/N_{cls} ITS;entries", {HistType::kTH1F, {{100, 0., 40.}}}},
      {"Trk/hItsNCls", ";N_{cls} ITS;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"Trk/hItsNClsInnerBarrel", ";N_{cls} ITS Inner Barrel;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"Trk/hDcaXY", ";DCA_{xy} (cm);entries", {HistType::kTH1F, {{200, -0.1, 0.1}}}},
      {"Trk/hDcaZ", ";DCA_{z} (cm);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}}},

      // track level (prongs of the candidate, after collision matching)
      {"TrkColl/hPt", ";#it{p}_{T} (GeV/#it{c});entries", hPtTrk},
      {"TrkColl/hEta", ";#eta;entries", hEtaTrk},
      {"TrkColl/hChi2NCl", ";#chi^{2}/N_{cls} TPC;entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"TrkColl/hTpcSignalVsP", ";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x}", {HistType::kTH2F, {{200, 0., 2.}, {300, 0., 300.}}}},
      {"TrkColl/hNSigmaVsP", ";#it{p} (GeV/#it{c});n#sigma^{TPC}", {HistType::kTH2F, {{200, 0., 2.}, {100, -10., 10.}}}},
      {"TrkColl/hHasIts", ";has ITS;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"TrkColl/hIsPvContrib", ";is PV contributor;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"TrkColl/hTpcNClsFound", ";N_{cls} TPC;entries", {HistType::kTH1F, {{160, 0., 160.}}}},
      {"TrkColl/hItsChi2NCl", ";#chi^{2}/N_{cls} ITS;entries", {HistType::kTH1F, {{100, 0., 40.}}}},
      {"TrkColl/hItsNCls", ";N_{cls} ITS;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"TrkColl/hItsNClsInnerBarrel", ";N_{cls} ITS Inner Barrel;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
      {"TrkColl/hDcaXY", ";DCA_{xy} (cm);entries", {HistType::kTH1F, {{200, -0.1, 0.1}}}},
      {"TrkColl/hDcaZ", ";DCA_{z} (cm);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}}}
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
  
    // candidate level plots
    registry.add("Cand/hMass", ";m_{#pi#pi} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
    registry.add("Cand/hPt", ";#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Cand/hRapidity", ";#it{y};entries", {HistType::kTH1F, {axisY}});
    registry.add("Cand/hPtVsRapidity", ";#it{p}_{T} (GeV/#it{c});#it{y}", {HistType::kTH2F, {axisPt, axisY}});

    // MC histograms are only booked if at least one MC process function is enabled
    const bool needMc = doprocessRhoTracksBeforeGroupingMcInfo || doprocessJpsiTracksBeforeGroupingMcInfo ||
                        doprocessRhoCandMcInfo || doprocessJpsiCandMcInfo || doprocessMcGenRho || doprocessMcGenJpsi;
    if (!needMc) {
      return;
    }

    // MC reco matched plots
    registry.add("Matched/Cand/hMass", ";m_{#pi#pi}^{MC} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
    registry.add("Matched/Cand/hMcPt", ";#it{p}_{T}^{MC} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Matched/Cand/hMcRapidity", ";#it{y}^{MC};entries", {HistType::kTH1F, {axisY}});
    registry.add("Matched/Cand/hPtVsRapidity", ";#it{p}_{T}^{MC matched} (GeV/#it{c});#it{y}^{MC}", {HistType::kTH2F, {axisPt, axisY}});

    registry.add("Matched/Coll/hNContrib", ";N_{PV contributors}^{MC matched};entries", {HistType::kTH1F, {{10, -0.5, 9.5}}});
    registry.add("Matched/Coll/hVtxX", ";#it{x}_{vtx}^{MC matched} (cm);entries", hVtxXY);
    registry.add("Matched/Coll/hVtxY", ";#it{y}_{vtx}^{MC matched} (cm);entries", hVtxXY);
    registry.add("Matched/Coll/hVtxZ", ";#it{z}_{vtx}^{MC matched} (cm);entries", hVtxZ);

    registry.add("Matched/Trk/hPt", ";#it{p}_{T}^{MC matched} (GeV/#it{c});entries", hPtTrk);
    registry.add("Matched/Trk/hEta", ";#eta^{MC};entries", hEtaTrk);
    registry.add("Matched/Trk/hMcIsPhysicalPrimary", ";is physical primary;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("Matched/Trk/hMotherPdg", ";mother PDG code;entries", {HistType::kTH1F, {{1001, 0.5, 1000.5}}});

    registry.add("Matched/TrkColl/hPt", ";#it{p}_{T}^{MC matched} (GeV/#it{c});entries", hPtTrk);
    registry.add("Matched/TrkColl/hEta", ";#eta^{MC};entries", hEtaTrk);
    registry.add("Matched/TrkColl/hMcIsPhysicalPrimary", ";is physical primary;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("Matched/TrkColl/hMotherPdg", ";mother PDG code;entries", {HistType::kTH1F, {{1001, 0.5, 1000.5}}});

    // MC generated plots
    registry.add("McGen/Coll/hVtxX", ";#it{x}_{vtx}^{MC} (cm);entries", hVtxXY);
    registry.add("McGen/Coll/hVtxY", ";#it{y}_{vtx}^{MC} (cm);entries", hVtxXY);
    registry.add("McGen/Coll/hVtxZ", ";#it{z}_{vtx}^{MC} (cm);entries", hVtxZ);

    registry.add("McGen/Part/hPt", ";#it{p}_{T}^{MC} (GeV/#it{c});entries", hPtTrk);
    registry.add("McGen/Part/hEta", ";#eta^{MC};entries", hEtaTrk);
    registry.add("McGen/Part/hMcIsPhysicalPrimary", ";is physical primary;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});

    registry.add("McGen/MotherPart/hMotherPdg", ";mother PDG code;entries", {HistType::kTH1F, {{1001, 0.5, 1000.5}}});
    registry.add("McGen/MotherPart/hPt", ";#it{p}_{T}^{MC} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("McGen/MotherPart/hRapidity", ";#it{y}^{MC};entries", {HistType::kTH1F, {axisY}});
    registry.add("McGen/MotherPart/hPtVsRapidity", ";#it{p}_{T}^{MC} (GeV/#it{c});#it{y}^{MC}", {HistType::kTH2F, {axisPt, axisY}});
  }

  template <CandSpecies species, typename TTrack>
  float getNSigma(TTrack const& track)
  {
    if constexpr (species == CandSpecies::kRho) {
      return track.tpcNSigmaPi();
    } else {
      return track.tpcNSigmaMu();
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

  // Basic kinematic / TPC / ITS properties of a candidate prong
  template <CandSpecies species, bool isMc, typename TTrack>
  void checkTpcTrackProperties(TTrack const& track)
  {
    registry.fill(HIST("TrkColl/hPt"), track.pt());
    registry.fill(HIST("TrkColl/hEta"), track.eta());
    registry.fill(HIST("TrkColl/hChi2NCl"), track.tpcChi2NCl());
    registry.fill(HIST("TrkColl/hTpcSignalVsP"), track.p(), track.tpcSignal());
    registry.fill(HIST("TrkColl/hNSigmaVsP"), track.p(), getNSigma<species>(track));
    registry.fill(HIST("TrkColl/hHasIts"), static_cast<int>(track.hasITS()));
    registry.fill(HIST("TrkColl/hIsPvContrib"), static_cast<int>(track.isPVContributor()));
    registry.fill(HIST("TrkColl/hTpcNClsFound"), track.tpcNClsFound());
    registry.fill(HIST("TrkColl/hItsChi2NCl"), track.itsChi2NCl());
    registry.fill(HIST("TrkColl/hItsNCls"), track.itsNCls());
    registry.fill(HIST("TrkColl/hItsNClsInnerBarrel"), track.itsNClsInnerBarrel());
    registry.fill(HIST("TrkColl/hDcaXY"), track.dcaXY());
    registry.fill(HIST("TrkColl/hDcaZ"), track.dcaZ());

    if constexpr (isMc) {
      if (!track.has_mcParticle()) {
        return;
      }
      auto mcPart = track.mcParticle();
      registry.fill(HIST("Matched/TrkColl/hPt"), mcPart.pt());
      registry.fill(HIST("Matched/TrkColl/hEta"), mcPart.eta());
      registry.fill(HIST("Matched/TrkColl/hMcIsPhysicalPrimary"), static_cast<int>(mcPart.isPhysicalPrimary()));
      if (mcPart.has_mothers()) {
        for (auto const& mcMother : mcPart.template mothers_as<aod::McParticles>()) {
          registry.fill(HIST("Matched/TrkColl/hMotherPdg"), std::abs(mcMother.pdgCode()));
        }
      }
    }

    // TODO: ambiguity (aod::AmbiguousTracks), ...
  }

  // Basic single-track selection used to build the candidate
  template <CandSpecies species, typename TTrack>
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
    return std::abs(getNSigma<species>(track)) <= nSigmaTpcMax;
  }

  // Loop on tracks before grouping by collision
  template <CandSpecies species, bool isMc, typename TTrack>
  void fillTrackPlotsBeforeGrouping(TTrack const& tracks)
  {
    for (auto const& track : tracks) {
      // select good tracks for the candidate
      if (!isGoodTrack<species>(track)) {
        continue;
      }

      registry.fill(HIST("Trk/hPt"), track.pt());
      registry.fill(HIST("Trk/hEta"), track.eta());
      registry.fill(HIST("Trk/hChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("Trk/hTpcSignalVsP"), track.p(), track.tpcSignal());
      registry.fill(HIST("Trk/hNSigmaVsP"), track.p(), getNSigma<species>(track));
      registry.fill(HIST("Trk/hHasIts"), static_cast<int>(track.hasITS()));
      registry.fill(HIST("Trk/hIsPvContrib"), static_cast<int>(track.isPVContributor()));
      registry.fill(HIST("Trk/hTpcNClsFound"), track.tpcNClsFound());
      registry.fill(HIST("Trk/hItsChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("Trk/hItsNCls"), track.itsNCls());
      registry.fill(HIST("Trk/hItsNClsInnerBarrel"), track.itsNClsInnerBarrel());
      registry.fill(HIST("Trk/hDcaXY"), track.dcaXY());
      registry.fill(HIST("Trk/hDcaZ"), track.dcaZ());

      if constexpr (isMc) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mcPart = track.mcParticle();
        registry.fill(HIST("Matched/Trk/hPt"), mcPart.pt());
        registry.fill(HIST("Matched/Trk/hEta"), mcPart.eta());
        registry.fill(HIST("Matched/Trk/hMcIsPhysicalPrimary"), static_cast<int>(mcPart.isPhysicalPrimary()));
        if (mcPart.has_mothers()) {
          for (auto const& mcMother : mcPart.template mothers_as<aod::McParticles>()) {
            registry.fill(HIST("Matched/Trk/hMotherPdg"), std::abs(mcMother.pdgCode()));
          }
        }
      }
    }
  }

  // Tracks are grouped by collision automatically
  template <CandSpecies species, bool isMc, typename TTrack, typename TCollision>
  void checkCandidateTracks(TCollision const& collision, TTrack const& tracks, float candMass)
  {
    registry.fill(HIST("Coll/hNContrib"), collision.numContrib());
    registry.fill(HIST("Coll/hBCid"), collision.bcId());
    registry.fill(HIST("Coll/hVtxZ"), collision.posZ());
    registry.fill(HIST("Coll/hVtxX"), collision.posX());
    registry.fill(HIST("Coll/hVtxY"), collision.posY());
    registry.fill(HIST("Coll/hVtxChi2"), collision.chi2());
    if constexpr (isMc) {
      if (collision.has_mcCollision()) {
        auto mcColl = collision.mcCollision();
        registry.fill(HIST("Matched/Coll/hNContrib"), collision.numContrib());
        registry.fill(HIST("Matched/Coll/hVtxX"), mcColl.posX());
        registry.fill(HIST("Matched/Coll/hVtxY"), mcColl.posY());
        registry.fill(HIST("Matched/Coll/hVtxZ"), mcColl.posZ());
      }
    }

    // Sequential ITS/TPC cuts on ALL tracks of the collision
    fillCutFlow(tracks);

    // Select tracks for the candidate
    std::vector<decltype(tracks.begin())> goodTracks;
    for (auto const& track : tracks) {
      if (isGoodTrack<species>(track)) {
        goodTracks.push_back(track);
      }
    }

    // Exactly two tracks with opposite charge -> candidate
    if (goodTracks.size() != 2) {
      return;
    }
    auto const& track0 = goodTracks[0];
    auto const& track1 = goodTracks[1];
    if (track0.sign() * track1.sign() >= 0) {
      return;
    }

    ROOT::Math::PxPyPzMVector p0(track0.px(), track0.py(), track0.pz(), candMass);
    ROOT::Math::PxPyPzMVector p1(track1.px(), track1.py(), track1.pz(), candMass);
    auto candidate = p0 + p1;

    if (candidate.M() < massMin || candidate.M() > massMax) {
      return;
    }
    if (candidate.Pt() < ptCandMin || std::abs(candidate.Rapidity()) > yCandMax) {
      return;
    }

    registry.fill(HIST("Cand/hMass"), candidate.M());
    registry.fill(HIST("Cand/hPt"), candidate.Pt());
    registry.fill(HIST("Cand/hRapidity"), candidate.Rapidity());
    registry.fill(HIST("Cand/hPtVsRapidity"), candidate.Pt(), candidate.Rapidity());

    if constexpr (isMc) {
      if (track0.has_mcParticle() && track1.has_mcParticle()) {
        auto mcPart0 = track0.mcParticle();
        auto mcPart1 = track1.mcParticle();
        ROOT::Math::PxPyPzMVector mcP0(mcPart0.px(), mcPart0.py(), mcPart0.pz(), candMass);
        ROOT::Math::PxPyPzMVector mcP1(mcPart1.px(), mcPart1.py(), mcPart1.pz(), candMass);
        auto mcCandidate = mcP0 + mcP1;
        registry.fill(HIST("Matched/Cand/hMass"), mcCandidate.M());
        registry.fill(HIST("Matched/Cand/hMcPt"), mcCandidate.Pt());
        registry.fill(HIST("Matched/Cand/hMcRapidity"), mcCandidate.Rapidity());
        registry.fill(HIST("Matched/Cand/hPtVsRapidity"), mcCandidate.Pt(), mcCandidate.Rapidity());
      }
    }

    checkTpcTrackProperties<species, isMc>(track0);
    checkTpcTrackProperties<species, isMc>(track1);
  }

  // Generated (MC truth) level: daughters and their mothers
  template <CandSpecies species>
  void checkMcGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    using o2::constants::physics::Pdg;
    constexpr int DaughterPdg = (species == CandSpecies::kRho) ? static_cast<int>(PDG_t::kPiPlus)
                                                               : static_cast<int>(PDG_t::kMuonMinus);
    constexpr int MotherPdg = (species == CandSpecies::kRho) ? static_cast<int>(PDG_t::kRho770_0)
                                                             : static_cast<int>(o2::constants::physics::Pdg::kJPsi);

    registry.fill(HIST("McGen/Coll/hVtxX"), mcCollision.posX());
    registry.fill(HIST("McGen/Coll/hVtxY"), mcCollision.posY());
    registry.fill(HIST("McGen/Coll/hVtxZ"), mcCollision.posZ());

    for (auto const& mcPart : mcParticles) {
      if (std::abs(mcPart.pdgCode()) != DaughterPdg) {
        continue;
      }
      if (mcPart.pt() < ptTrackMin || std::abs(mcPart.eta()) > etaTrackMax) {
        continue;
      }

      registry.fill(HIST("McGen/Part/hPt"), mcPart.pt());
      registry.fill(HIST("McGen/Part/hEta"), mcPart.eta());
      registry.fill(HIST("McGen/Part/hMcIsPhysicalPrimary"), static_cast<int>(mcPart.isPhysicalPrimary()));

      if (!mcPart.has_mothers()) {
        continue;
      }
      for (auto const& mcMother : mcPart.template mothers_as<aod::McParticles>()) {
        registry.fill(HIST("McGen/MotherPart/hMotherPdg"), std::abs(mcMother.pdgCode()));

        // fill the mother once, not once per daughter: only via the positive-PDG prong
        if (std::abs(mcMother.pdgCode()) == MotherPdg && mcPart.pdgCode() > 0) {
          registry.fill(HIST("McGen/MotherPart/hPt"), mcMother.pt());
          registry.fill(HIST("McGen/MotherPart/hRapidity"), mcMother.y());
          registry.fill(HIST("McGen/MotherPart/hPtVsRapidity"), mcMother.pt(), mcMother.y());
        }
      }
    }
  }

  // Process functions
  void processRhoTracksBeforeGrouping(TracksExtraWPidPi const& tracks)
  {
    fillTrackPlotsBeforeGrouping<CandSpecies::kRho, false>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoTracksBeforeGrouping, "Process rho tracks before asking for collisions", true);

  void processJpsiTracksBeforeGrouping(TracksExtraWPidMu const& tracks)
  {
    fillTrackPlotsBeforeGrouping<CandSpecies::kJpsi, false>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiTracksBeforeGrouping, "Process J/Psi tracks before asking for collisions", false);

  void processRhoTracksBeforeGroupingMcInfo(TracksExtraWPidPiMc const& tracks, aod::McParticles const&)
  {
    fillTrackPlotsBeforeGrouping<CandSpecies::kRho, true>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoTracksBeforeGroupingMcInfo, "Process rho tracks before asking for collisions, ask for MC info", false);

  void processJpsiTracksBeforeGroupingMcInfo(TracksExtraWPidMuMc const& tracks, aod::McParticles const&)
  {
    fillTrackPlotsBeforeGrouping<CandSpecies::kJpsi, true>(tracks);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiTracksBeforeGroupingMcInfo, "Process J/Psi tracks before asking for collisions, ask for MC info", false);

  void processRhoCand(aod::Collision const& collision, TracksExtraWPidPi const& tracks)
  {
    checkCandidateTracks<CandSpecies::kRho, false>(collision, tracks, MassPion);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoCand, "Rho -> pi pi candidates and track QA", true);

  void processJpsiCand(aod::Collision const& collision, TracksExtraWPidMu const& tracks)
  {
    checkCandidateTracks<CandSpecies::kJpsi, false>(collision, tracks, MassMuon);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiCand, "J/Psi -> mu mu candidates and track QA", false);

  void processRhoCandMcInfo(CollisionsWithMc::iterator const& collision, aod::McCollisions const&, TracksExtraWPidPiMc const& tracks, aod::McParticles const&)
  {
    checkCandidateTracks<CandSpecies::kRho, true>(collision, tracks, MassPion);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processRhoCandMcInfo, "Rho -> pi pi candidates and track QA, ask for MC info", false);

  void processJpsiCandMcInfo(CollisionsWithMc::iterator const& collision, aod::McCollisions const&, TracksExtraWPidMuMc const& tracks, aod::McParticles const&)
  {
    checkCandidateTracks<CandSpecies::kJpsi, true>(collision, tracks, MassMuon);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processJpsiCandMcInfo, "J/Psi -> mu mu candidates and track QA, ask for MC info", false);

  void processMcGenRho(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    checkMcGen<CandSpecies::kRho>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processMcGenRho, "Generated rho -> pi pi", false);

  void processMcGenJpsi(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    checkMcGen<CandSpecies::kJpsi>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(UpcTrackVertexingQA, processMcGenJpsi, "Generated J/Psi -> mu mu", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcTrackVertexingQA>(cfgc)};
}
