#include "BphAna/BToKstarMuMuRun2/interface/CommonBDecayAssistant.h"

#include <set>

#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

CommonBDecayAssistant::CommonBDecayAssistant(const edm::ParameterSet& ps){
    // Remember to add description to the ParameterSet used in EDAnalyzer.
    // Labels to retrieve objects
    gpLabel = ps.getParameter<edm::InputTag>("gpLabel");
    ccLabel = ps.getParameter<edm::InputTag>("ccLabel");
    bsLabel = ps.getParameter<edm::InputTag>("bsLabel");
    pvLabel = ps.getParameter<edm::InputTag>("pvLabel");
    hltLabel= ps.getParameter<edm::InputTag>("hltLabel");
    l1mLabel= ps.getParameter<edm::InputTag>("l1mLabel");

    // parameters to be used
    triggerNames    = ps.getParameter< std::vector<std::string> >("triggerNames");
    lastFilterNames = ps.getParameter< std::vector<std::string> >("lastFilterNames");

    // initialization
    hltConfigPtr_ = new HLTConfigProvider();
}

CommonBDecayAssistant::~CommonBDecayAssistant(){
    delete hltConfigPtr_;
    hltConfigPtr_ = 0;
}

void CommonBDecayAssistant::updateHandles(const edm::Event& ev, const edm::EventSetup& es){

    bsToken.get(ev, bsHandle);
    gpToken.get(ev, gpHandle);
    ccToken.get(ev, ccHandle);
    pvToken.get(ev, pvHandle);
    hltToken.get(ev,hltHandle);
    l1mToken.get(ev,l1mHandle);
    //es.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    // Corresponding objects for getters should be cleaned
    // All following updates at the first time the getter is called
    bsPtr_      = 0;
    bsVtxPtr_   = 0;
    pvPtr_      = 0;
    trkPtr_     = 0;
    muPtr_      = 0;
    mu_.clear();
    trk_.clear();
}

const reco::BeamSpot* CommonBDecayAssistant::getBS(){
    if (bsPtr_ ==0){
        if ( !bsHandle.isValid() ){
            printf("DEBUG\t: (CommonBDecayAssistant) No beam spot available from EventSetup\n");
        }else{
            bsPtr_= &*bsHandle;
        }
    }
    return bsPtr_;
}

reco::Vertex* CommonBDecayAssistant::getBSVtx(){
    if (bsVtxPtr_ == 0 && getBS() != 0){
        bsVtx_ = reco::Vertex(bsPtr_->position(),bsPtr_->covariance3D());
        bsVtxPtr_ = &bsVtx_;
    }
    return bsVtxPtr_;
}

reco::Vertex* CommonBDecayAssistant::getPV(){
    if (pvPtr_ == 0){
        for(std::vector<reco::Vertex>::const_iterator iVtx = pvHandle->begin();iVtx != pvHandle->end(); iVtx++ ){
            if (iVtx->isValid()){
                pv_ = *iVtx;
                break;
            }
        }
        if (pvPtr_ == 0) pvPtr_ = getBSVtx();
    }
    return pvPtr_;
}

std::vector<const pat::GenericParticle*>* CommonBDecayAssistant::getTrk(){
    if (trkPtr_ == 0){
        int nCand = gpHandle->size();
        trk_.clear();
        trk_.reserve(nCand);
        
        for(int iCand = 0; iCand < nCand; ++iCand){
            const pat::GenericParticle& cand = gpHandle->at( iCand );
            if ( !trkPresel_(&cand) ) continue;
            trk_.push_back(&cand);
        }
        trkPtr_ = &trk_;
    }
    return trkPtr_;
}

bool CommonBDecayAssistant::trkPresel_(const pat::GenericParticle* trk){
    if ( trk->pt() < 0.7 ) return false;
    if ( fabs(trk->eta()) > 2.5 ) return false;
    return true;
}

std::vector<const pat::Muon*>* CommonBDecayAssistant::getMu(){
    if (muPtr_ == 0){
        // This section taken from BPHSkim example code
        // Scan all composite candidate daughters and fill non-duplicated muon to the vector.
        std::set<const pat::Muon*> muonSet;
    
        int nCCand = ccHandle->size();
        mu_.reserve(nCCand);
        muonSet.clear();
        std::set<const pat::Muon*>::const_iterator iter;
        std::set<const pat::Muon*>::const_iterator iend;
        for(int iCCand = 0; iCCand < nCCand; ++iCCand){
            const pat::CompositeCandidate& cCand = ccHandle->at( iCCand );
            int nCCandDau = cCand.numberOfDaughters();
            for ( int iCCandDau = 0; iCCandDau < nCCandDau; ++iCCandDau ) {
                const reco::Candidate* dp = cCand.daughter( iCCandDau );
                const pat::Muon* mp = dynamic_cast<const pat::Muon*>( dp );
                if ( !muPresel_(mp) ) continue;
                iter = muonSet.begin();
                iend = muonSet.end();
                bool add = ( mp != 0 ) && ( muonSet.find( mp ) == iend );
                while ( add && ( iter != iend ) ) {
                    if ( BPHRecoBuilder::sameTrack( mp, *iter++, 1.0e-5 ) ) add = false;
                }
                if ( add ) muonSet.insert( mp );
            }    
            iter = muonSet.begin();
            iend = muonSet.end();
            while ( iter != iend ) mu_.push_back( *iter++ );
        }
        muPtr_ = &mu_;
    }
    return muPtr_;
}

bool CommonBDecayAssistant::muPresel_(const pat::Muon* mu){
    // BPH soft muon condition.
    return true;
}

void CommonBDecayAssistant::updateHltConfig(edm::Run const& iRun, edm::EventSetup const& iSetup){
    bool dummy(true);
    hltConfigPtr_->init(iRun, iSetup, "CommonBDecayAssistant", dummy);
}

void CommonBDecayAssistant::writePSetDescription(edm::ParameterSetDescription &desc){
    // Ref: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideConfigurationValidationAndHelp

    desc.add<edm::InputTag>("gpLabel" ,edm::InputTag("patSelectedTracks"));
    desc.add<edm::InputTag>("ccLabel" ,edm::InputTag("onia2MuMuPAT::RECO"));
    desc.add<edm::InputTag>("bsLabel" ,edm::InputTag("offlineBeamSpot"));
    desc.add<edm::InputTag>("pvLabel" ,edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("hltLabel",edm::InputTag("TriggerResults"));
    desc.add<edm::InputTag>("l1mLabel",edm::InputTag("gtDigis"));
    desc.add<std::vector<std::string> >("triggerNames");
    desc.add<std::vector<std::string> >("lastFilterNames");

    // In case nested config is used.
    //edm::ParameterSetDescription muonDesc; 
    //muonDesc.addUntracked<int>("maxEta",2.5);
    //desc.addUntracked<edm::ParameterSetDescription>("muonPS",muonDesc);

    return;
}

