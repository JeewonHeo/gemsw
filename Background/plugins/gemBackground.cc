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

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"

#include "TFile.h"
#include "TTree.h"

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class gemBackground : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit gemBackground(const edm::ParameterSet&);
  ~gemBackground() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<OnlineLuminosityRecord> onlineLumiRecord_;

  edm::Service<TFileService> fs_;

  TTree *t_rec_hits;
  float b_instLumi;
  int b_event, b_eventTime, b_bunchId, b_orbitNumber;
  int b_region, b_layer, b_chamber, b_ieta;
  int b_first_strip, b_cluster_size;
};

gemBackground::gemBackground(const edm::ParameterSet& iConfig) {
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

// ------------ method called for each event  ------------
void gemBackground::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<OnlineLuminosityRecord> onlineLumiRecord;
  iEvent.getByToken(onlineLumiRecord_, onlineLumiRecord);

  // std::cout << "TIME: " << iEvent.time().unixTime()
  //     << " EVENTID: " << iEvent.id()
  //     << " BX: " << iEvent.bunchCrossing()
  //     << " ORBITNUMBER: " << iEvent.orbitNumber()
  //     << " INST.LUMI.: " << onlineLumiRecord->instLumi()
  //     << std::endl;

  b_instLumi = onlineLumiRecord->instLumi();
  b_bunchId = iEvent.bunchCrossing();
  b_orbitNumber = iEvent.orbitNumber();
  b_eventTime = iEvent.time().unixTime();
  b_event = iEvent.id().event();

  for (auto rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit) {
    GEMDetId gem_id = rechit->gemId();
    if (!gem_id.isGE11()) continue;
    b_region = gem_id.region();
    b_layer = gem_id.layer();
    b_chamber = gem_id.chamber();
    b_ieta = gem_id.ieta();
    b_first_strip = rechit->firstClusterStrip();
    b_cluster_size = rechit->clusterSize();

    t_rec_hits->Fill();

    // std::cout << "region: " << b_region
    //     << " layer: " << b_layer
    //     << " chamber: " << b_chamber
    //     << " ieta: " << b_ieta
    //     << " first_strip: " << b_first_strip
    //     << " cluster_size: " << b_cluster_size
    //     << std::endl;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void gemBackground::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void gemBackground::endJob() {
  // please remove this method if not needed
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
