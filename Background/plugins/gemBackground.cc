// -*- C++ -*-
//
// Package:    gem-background/gemBackground
// Class:      gemBackground
//
/**\class gemBackground gemBackground.cc gem-background/gemBackground/plugins/gemBackground.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Master student from University of Seoul
//         Created:  Mon, 03 Oct 2022 08:17:53 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMOHStatusCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


typedef std::tuple<int, int> Key2;

class gemBackground : public edm::one::EDAnalyzer<> {
public:
  explicit gemBackground(const edm::ParameterSet&);
  ~gemBackground() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool maskChamberWithError(const GEMDetId& chamber_id, const edm::Handle<GEMOHStatusCollection>);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);// override;
  virtual void beginJob() override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMOHStatusCollection> oh_status_collection_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<OnlineLuminosityRecord> onlineLumiRecord_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_;

  std::map<Key2, TH1D*> bad_chamber_;
  edm::Service<TFileService> fs_;

  TTree *t_rec_hits;
  float b_instLumi;
  int b_event, b_eventTime, b_bunchId, b_orbitNumber;
  int b_region, b_layer, b_chamber, b_ieta;
  int b_first_strip, b_cluster_size;
};

gemBackground::gemBackground(const edm::ParameterSet& iConfig) {
  gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("muonGEMDigis"));
  oh_status_collection_ = consumes<GEMOHStatusCollection>(iConfig.getParameter<edm::InputTag>("OHInputLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  onlineLumiRecord_ = consumes<OnlineLuminosityRecord>(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"));
  t_rec_hits = fs_->make<TTree>("rec_hits", "gem_rec_hits");
  #define BRANCH_(name, suffix) t_rec_hits->Branch(#name, & b_##name, #name "/" #suffix);
  BRANCH_(instLumi, F);
  BRANCH_(event, I);
  BRANCH_(eventTime, I);
  BRANCH_(bunchId, I);
  BRANCH_(orbitNumber, I);
  BRANCH_(region, I);
  BRANCH_(layer, I);
  BRANCH_(chamber, I);
  BRANCH_(ieta, I);
  BRANCH_(first_strip, I);
  BRANCH_(cluster_size, I);
}

gemBackground::~gemBackground() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

bool gemBackground::maskChamberWithError(const GEMDetId& chamber_id,
                                         // const GEMOHStatusCollection* oh_status_collection) {
                                         const edm::Handle<GEMOHStatusCollection> oh_status_collection) {
  const bool mask = true;
  for (auto iter = oh_status_collection->begin(); iter != oh_status_collection->end(); iter++) {
    const auto [oh_id, range] = (*iter);
    if (chamber_id != oh_id) {
      continue;
    }

    for (auto oh_status = range.first; oh_status != range.second; oh_status++) {
      if (oh_status->isBad()) {
        // GEMOHStatus is bad. Mask this chamber.
        // std::cout << "Chamber Id " << oh_id << " isBad " << oh_status->isBad() << std::endl;
        return mask;
      }  // isBad
    }  // range
  }  // collection
  return not mask;
}



// ------------ method called for each event  ------------
void gemBackground::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<GEMOHStatusCollection> oh_status_collection;
  iEvent.getByToken(oh_status_collection_, oh_status_collection);

  edm::Handle<OnlineLuminosityRecord> onlineLumiRecord;
  iEvent.getByToken(onlineLumiRecord_, onlineLumiRecord);

  b_instLumi = onlineLumiRecord->instLumi();
  b_bunchId = iEvent.bunchCrossing();
  b_orbitNumber = iEvent.orbitNumber();
  b_eventTime = iEvent.time().unixTime();
  b_event = iEvent.id().event();

  for (auto rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit) {
    const GEMDetId gem_id = rechit->gemId();
    if (!gem_id.isGE11()) {
      continue;
    }
    GEMDetId chamber_id = gem_id.chamberId();
    if (maskChamberWithError(chamber_id, oh_status_collection)) {
      //  std::cout << "BAD: " << chamber_id << std::endl;
      b_region = gem_id.region();
      b_layer = gem_id.layer();
      b_chamber = gem_id.chamber();
      Key2 key2{b_region, b_layer};
      bad_chamber_[key2]->Fill(b_chamber);
      continue;
    }

    b_region = gem_id.region();
    b_layer = gem_id.layer();
    b_chamber = gem_id.chamber();
    b_ieta = gem_id.ieta();
    b_first_strip = rechit->firstClusterStrip();
    b_cluster_size = rechit->clusterSize();

    t_rec_hits->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void gemBackground::beginJob() {
  int st = 1;
  for (int re = -1; re <= 1; re += 2) {
    for (int la = 1; la <= 2; ++la) {
      Key2 key2{re, la};
      bad_chamber_[key2] = fs_->make<TH1D>(Form("bad_chamber_GE%d/%d-L%d", re, st, la),
                                           Form("Bad Chambers GE%d/%d-L%d", re, st, la),
                                           36, 0.5, 36.5);
    }
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void gemBackground::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void gemBackground::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(gemBackground);
