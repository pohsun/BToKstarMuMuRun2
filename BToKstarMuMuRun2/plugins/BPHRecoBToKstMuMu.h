#ifndef BPHRECOLAMBTOJPSILAM_H
#define BPHRECOLAMBTOJPSILAM_H

#include "BPHAnalysis/RecoDecay/interface/BPHAnalyzerTokenWrapper.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include "BphAna/BToKstarMuMuRun2/interface/CommonBDecayAssistant.h"
#include "BphAna/BToKstarMuMuRun2/interface/EvtBranchPlanter.h"
#include "BphAna/BToKstarMuMuRun2/interface/MuonBranchPlanter.h"

class BPHRecoBToKstMuMu: public BPHAnalyzerWrapper<BPHModuleWrapper::one_analyzer> {
    
public:
    
    explicit BPHRecoBToKstMuMu( const edm::ParameterSet& ps );
    virtual ~BPHRecoBToKstMuMu() {};
    
    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );
    
    virtual void beginJob();
    virtual void beginRun(const edm::Run&, const edm::EventSetup&  );
    virtual void analyze( edm::Event const& , edm::EventSetup const& );
    virtual void endJob();
 
    friend class CommonBDecayAssistant;
private:
    edm::Service<TFileService> fs;

    // CommonTools and specified treePlanters   
    CommonBDecayAssistant   assistant;
    TTree *evtTree;
    EvtBranchPlanter        EvtInfo; 
    MuonBranchPlanter       MuonInfo;

};
#endif
