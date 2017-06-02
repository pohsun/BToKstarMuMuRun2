#ifndef COMMON_BDECAY_ASSISTANT_H
#define COMMON_BDECAY_ASSISTANT_H

#include <string>


#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "BPHAnalysis/RecoDecay/interface/BPHAnalyzerTokenWrapper.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

/*

   This class handles all routine stuff for an analysis.
   The edanalyzer should inherite BPHAnalyzerWrapper<BPHModuleWrapper::one_analyzer>

*/

// Global constant
const int MUONMINUS_PDG_ID          = 13;
const int PIONPLUS_PDG_ID           = 211;
const int KLONGZERO_PDG_ID          = 130;
const int KSHORTZERO_PDG_ID         = 310;
const int KZERO_PDG_ID              = 311;
const int KSTARPLUS_PDG_ID          = 323;
const int BPLUS_PDG_ID              = 521;
const int JPSI_PDG_ID               = 443;
const int PSI2S_PDG_ID              = 100443;
const int ETA_PDG_ID                = 221;
const int DZERO_PDG_ID              = 421;
const int DSTAR2007_PDG_ID          = 423;
const int DSPLUS_PDG_ID             = 431;
const int DS1PLUS2460_PDG_ID        = 20443;
const int SIGMACSTARPLUSPLUS_PDG_ID = 4224;
const int DELTAPLUS_PDG_ID          = 2214;

const double PI                     = 3.141592653589793;

// 
class CommonBDecayAssistant{
    public:
        CommonBDecayAssistant(const edm::ParameterSet& ps);
        ~CommonBDecayAssistant();

        // SOP for building the ntuple
        // Tokens(Label) -> Handle -> Object
            // Tokens are retrieved at constructor level
            // consumes, mayConsume is inherited from EDM, treat them from outside.
            // Handles are updated in for every event.
            // Handles and getters are used as input for branch planters.
        void updateHandles(const edm::Event&, edm::EventSetup const&);
        void updateHltConfig(edm::Run const&, edm::EventSetup const&);// Put this in "beginRun" for efficiency purpose

        // adders, setters, and getters
        const reco::BeamSpot*                       getBS();
        reco::Vertex*                               getBSVtx();
        reco::Vertex*                               getPV();
        HLTConfigProvider*                          getHltConfig(){return hltConfigPtr_;}
        std::vector<const pat::GenericParticle*>*   getTrk();
        std::vector<const pat::Muon*>*              getMu();
        //std::vector<pat::CompositeCandidate>*       getKshort();
        
        // Labels and PSets retrived in constructor
        edm::InputTag   gpLabel;
        edm::InputTag   ccLabel;
        edm::InputTag   bsLabel;
        edm::InputTag   pvLabel;
        edm::InputTag   hltLabel;
        edm::InputTag   l1mLabel; 
        std::vector<std::string> triggerNames;
        std::vector<std::string> lastFilterNames;

        // Tokens
        BPHTokenWrapper< std::vector<pat::GenericParticle>      > gpToken;
        BPHTokenWrapper< std::vector<pat::CompositeCandidate>   > ccToken;
        BPHTokenWrapper< reco::BeamSpot                         > bsToken;
        BPHTokenWrapper< reco::VertexCollection                 > pvToken;
        BPHTokenWrapper< edm::TriggerResults                    > hltToken;
        BPHTokenWrapper< std::vector<L1MuGMTCand>               > l1mToken;

        // Handles
        edm::ESHandle<MagneticField>                        bFieldHandle;
        edm::Handle< reco::BeamSpot >                       bsHandle;
        edm::Handle< reco::VertexCollection >               pvHandle;
        edm::Handle< edm::TriggerResults >                  hltHandle;
        edm::Handle< std::vector<pat::GenericParticle> >    gpHandle;
        edm::Handle< std::vector<pat::CompositeCandidate> > ccHandle;
        edm::Handle< std::vector<L1MuGMTCand> >             l1mHandle;

        // services
        virtual const std::string getClassName() const {return "CommonBDecayAssistant";}
        virtual const std::string getMotherClassName() const {return "";}
        static void writePSetDescription(edm::ParameterSetDescription&);


        // friend
        friend class EvtBranchPlanter;
        friend class MuonBranchPlanter;

    private:
        // Corresponding objects for getters
        // Make sure getters  before accessing.
        const reco::BeamSpot*   bsPtr_ = 0;
        reco::Vertex*           bsVtxPtr_ = 0;
        reco::Vertex            bsVtx_;
        reco::Vertex*           pvPtr_ = 0;
        reco::Vertex            pv_;
        std::vector<const pat::Muon*>*              muPtr_= 0;
        std::vector<const pat::Muon*>               mu_;
        std::vector<const pat::GenericParticle*>*   trkPtr_ = 0;
        std::vector<const pat::GenericParticle*>    trk_;
        HLTConfigProvider *hltConfigPtr_ = 0;

        // object preslection applied
        // TODO: Taking cut value from PSet
        bool muPresel_(const pat::Muon*);
        bool trkPresel_(const pat::GenericParticle*);
};

#endif
