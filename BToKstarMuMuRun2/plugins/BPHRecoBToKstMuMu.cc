#include "BphAna/BToKstarMuMuRun2/plugins/BPHRecoBToKstMuMu.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

BPHRecoBToKstMuMu::BPHRecoBToKstMuMu( const edm::ParameterSet& ps ):
    assistant(ps)
{
    consume< std::vector<pat::GenericParticle>      >( assistant.gpToken , assistant.gpLabel );
    consume< std::vector<pat::CompositeCandidate>   >( assistant.ccToken , assistant.ccLabel );
    consume< reco::BeamSpot                         >( assistant.bsToken , assistant.bsLabel );
    consume< reco::VertexCollection                 >( assistant.pvToken , assistant.pvLabel );
    consume< edm::TriggerResults                    >( assistant.hltToken, assistant.hltLabel);
    consume< std::vector<L1MuGMTCand>               >( assistant.l1mToken, assistant.l1mLabel);

    EvtInfo .linkAssistant(&assistant);
    MuonInfo.linkAssistant(&assistant);

    LogDebug("BPHRecoBToKstMuMu") << "End of constructor";
}

void BPHRecoBToKstMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions ) {
    edm::LogInfo("BPHRecoBToKstMuMu") << "FillDescriptions";
    edm::ParameterSetDescription desc;
    CommonBDecayAssistant::writePSetDescription(desc);
    descriptions.add( "BPHRecoBToKstMuMu", desc );
    return;
}

void BPHRecoBToKstMuMu::beginJob() {
    LogDebug("BPHRecoBToKstMuMu") << "beginJob";
    evtTree     = fs->make<TTree>("evtTree" ,"evtTree");
    EvtInfo     .cuttage(evtTree,"EvtInfo");
    MuonInfo    .cuttage(evtTree,"MuonInfo");
}

void BPHRecoBToKstMuMu::analyze( const edm::Event& ev, const edm::EventSetup& es ) {
    LogDebug("BPHRecoBToKstMuMu") << "analyze";

    assistant.updateHandles(ev,es);
    EvtInfo     .plant(ev,es);
    MuonInfo    .plant(ev,es);
    //TrkInfo     .plant(ev,es);
    //OniaInfo    .plant(ev.es);
    //GenInfo     .plant(ev,es);

    evtTree->Fill();
    
    // clean up
    //memset(&JpsiBr,0x00,sizeof(JpsiBr));
    //memset(&LamBr,0x00,sizeof(LamBr));
    //memset(&LambBr,0x00,sizeof(LambBr));
    
    //get beamspot information
    //edm::Handle< reco::BeamSpot > bsHandle;
    //bsToken.get( ev, bsHandle );
    
    //reco::Vertex BeamSpotVtx;
    //reco::BeamSpot beamSpot;
    //if (bsHandle.isValid()){
    //    beamSpot = *bsHandle;
    //    BeamSpotVtx = reco::Vertex(beamSpot.position(), beamSpot.covariance3D());
    //}else{
    //    cout<< "No beam spot available from EventSetup" << endl;
    //}
    
    // get object collections
    // collections are got through "BPHTokenWrapper" interface to allow
    // uniform access in different CMSSW versions
    
    // get pat::GenericParticle collection (in skimmed data)
    //edm::Handle< vector<pat::GenericParticle> > gpCands;
    //gpCandsToken.get( ev, gpCands );
    
    // get muons from pat::CompositeCandidate objects describing onia;
    // muons from all composite objects are copied to an unique std::vector
    //vector<const reco::Candidate*> muDaugs;
    //set<const pat::Muon*> muonSet;

    //edm::Handle< vector<pat::CompositeCandidate> > ccCands;
    //ccCandsToken.get( ev, ccCands );
    //int n = ccCands->size();
    //muDaugs.clear();
    //muDaugs.reserve( n );
    //muonSet.clear();
    //set<const pat::Muon*>::const_iterator iter;
    //set<const pat::Muon*>::const_iterator iend;
    //int i;
    //for ( i = 0; i < n; ++i ) {
    //    const pat::CompositeCandidate& cc = ccCands->at( i );
    //    int j;
    //    int m = cc.numberOfDaughters();
    //    for ( j = 0; j < m; ++j ) {
    //        const reco::Candidate* dp = cc.daughter( j );
    //        const pat::Muon* mp = dynamic_cast<const pat::Muon*>( dp );
    //        iter = muonSet.begin();
    //        iend = muonSet.end();
    //        bool add = ( mp != 0 ) && ( muonSet.find( mp ) == iend );
    //        while ( add && ( iter != iend ) ) {
    //            if ( BPHRecoBuilder::sameTrack( mp, *iter++, 1.0e-5 ) ) add = false;
    //        }
    //        if ( add ) muonSet.insert( mp );
    //    }
    //    iter = muonSet.begin();
    //    iend = muonSet.end();
    //    while ( iter != iend ) muDaugs.push_back( *iter++ );
    //}
    
    // -----------------------------------
    // starting objects selection
    // -----------------------------------
    
    // muon selections
    class MuonSelect: public BPHRecoSelect {
        public:
            MuonSelect( double pt, double eta ): ptCut(pt), etaCut(eta) {}
            virtual bool accept( const reco::Candidate& cand ) const {
                const pat::Muon* p = reinterpret_cast<const pat::Muon*>( &cand );
                if ( p == 0 ) return false;
                return ( p->p4().pt()>ptCut && fabs(p->p4().eta())<etaCut );
            }
        private:
            double ptCut, etaCut;
    };
    
    // track selections
    class TrackSelect: public BPHRecoSelect {
        public:
            TrackSelect( double pt, double eta ): ptCut(pt), etaCut(eta) {}
            virtual bool accept( const reco::Candidate& cand ) const {
                return ( cand.p4().pt()>ptCut && fabs(cand.p4().eta())<etaCut );
            }
        private:
            double ptCut, etaCut;
    };
    
    // -----------------------------------
    // reconstructed object selection
    // -----------------------------------
    
    // selection by mass
    class MassSelect: public BPHMomentumSelect {
        public:
            MassSelect( double minMass, double maxMass ): mMin( minMass ), mMax( maxMass ) {}
            virtual bool accept( const BPHDecayMomentum& cand ) const {
                double mass = cand.composite().mass();
                return ( (mass>mMin) && (mass<mMax) );
            }
        private:
            double mMin, mMax;
    };
    
    // selection by vertex prob
    class VtxProbSelect: public BPHVertexSelect {
        public:
            VtxProbSelect( double minProb ): mProb( minProb ) {}
            virtual bool accept( const BPHDecayVertex& cand ) const {
                const reco::Vertex& v = cand.vertex();
                if ( v.isFake() ) return false;
                if ( !v.isValid() ) return false;
                return ( TMath::Prob( v.chi2(), lround( v.ndof() ) ) > mProb );
            }
        private:
            double mProb;
    };
    
    // build J/psi
    //BPHRecoBuilder bJPsi( es );
    //if ( useCC ) {
    //    bJPsi.add( "Muon1", BPHRecoBuilder::createCollection( muDaugs, "cfmig" ), 0.105658 );
    //    bJPsi.add( "Muon2", BPHRecoBuilder::createCollection( muDaugs, "cfmig" ), 0.105658 );
    //}
    
    //bJPsi.filter( "Muon1", MuonSelect(4.0, 2.1) ); // pt>4.0, |eta|<2.1
    //bJPsi.filter( "Muon2", MuonSelect(4.0, 2.1) ); // pt>4.0, |eta|<2.1
    //bJPsi.filter( MassSelect(2.5, 3.9) ); // 2.5<M(mumu)<3.9
    //bJPsi.filter( VtxProbSelect(0.0) ); // vertex prob > 0.0
    
    //vector<BPHPlusMinusConstCandPtr> lJPsi = BPHPlusMinusCandidate::build( bJPsi, "Muon1", "Muon2" );
    
    //// full the J/psi tree
    //for ( unsigned int iJPsi = 0; iJPsi < lJPsi.size(); iJPsi++ ) {
    //    const BPHRecoCandidate* cand = lJPsi[iJPsi].get();
        
    //    const reco::Vertex& vx = cand->vertex();
    //    TVector3 bs_to_vtx(vx.position().X()-BeamSpotVtx.position().X(),
    //                       vx.position().Y()-BeamSpotVtx.position().Y(),
    //                       vx.position().Z()-BeamSpotVtx.position().Z());
        
    //    JpsiBr.mass    = cand->composite().mass();
    //    JpsiBr.pt      = cand->composite().pt();
    //    JpsiBr.eta     = cand->composite().eta();
    //    JpsiBr.phi     = cand->composite().phi();
        
    //    JpsiBr.vtxprob = TMath::Prob( vx.chi2(), vx.ndof() );
    //    JpsiBr.fl2d    = bs_to_vtx.Perp();
    //    JpsiBr.cosa2d  = (cand->composite().px()*bs_to_vtx.X()+cand->composite().py()*bs_to_vtx.Y())/bs_to_vtx.Perp()/cand->composite().pt();
        
    //    JpsiTree->Fill();
    //}

    // build Lambda
    
    //BPHRecoBuilder bLam( es );
    //if ( useGP ) {
    //    bLam.add( "Proton", BPHRecoBuilder::createCollection( gpCands ), 0.938272 );
    //    bLam.add( "Pion",   BPHRecoBuilder::createCollection( gpCands ), 0.139570 );
    //}
    
    //bLam.filter( "Proton", TrackSelect(0.7, 2.5) ); // pt>0.7, |eta|<2.5
    //bLam.filter( "Pion", TrackSelect(0.7, 2.5) ); // pt>0.7, |eta|<2.5
    //bLam.filter( MassSelect(1.06, 1.16) ); // 1.06<M(p,pi)<1.16
    
    //vector<BPHPlusMinusConstCandPtr> lLam = BPHPlusMinusCandidate::build( bLam, "Proton", "Pion" );
    
    //// full the Lambda tree
    //for ( unsigned int iLam = 0; iLam < lLam.size(); iLam++ ) {
    //    const BPHRecoCandidate* cand = lLam[iLam].get();
        
    //    // perform vertex fit
    //    vector<RefCountedKinematicParticle> kComp = cand->kinParticles();
        
    //    KinematicParticleVertexFitter vtxFitter;
    //    RefCountedKinematicTree compTree = vtxFitter.fit( kComp );
    //    if (compTree->isEmpty()) continue;
    //    compTree->movePointerToTheTop();
    //    const RefCountedKinematicParticle kPart = compTree->currentParticle();
    //    const RefCountedKinematicVertex kVtx = compTree->currentDecayVertex();
    //    const KinematicState& kState = kPart->currentState();
    //    if (!kState.isValid()) continue;
        
    //    double vtxProb = TMath::Prob( kVtx->chiSquared(), kVtx->degreesOfFreedom() );
    //    if (vtxProb <=0.) continue; // fit prob > 0.0
        
    //    TVector3 v3_fit(kState.kinematicParameters().momentum().x(),
    //                    kState.kinematicParameters().momentum().y(),
    //                    kState.kinematicParameters().momentum().z());

    //    TVector3 bs_to_vtx(kVtx->position().x()-BeamSpotVtx.position().X(),
    //                       kVtx->position().y()-BeamSpotVtx.position().Y(),
    //                       kVtx->position().z()-BeamSpotVtx.position().Z());
        
    //    const reco::Candidate* cand_p  = cand->getDaug( "Proton" );
    //    const reco::Candidate* cand_pi = cand->getDaug( "Pion" );
        
    //    LamBr.mass    = cand->composite().mass();
    //    LamBr.pt      = cand->composite().pt();
    //    LamBr.eta     = cand->composite().eta();
    //    LamBr.phi     = cand->composite().phi();
        
    //    LamBr.vtxmass = kState.mass();
    //    LamBr.vtxprob = vtxProb;
    //    LamBr.fl2d    = bs_to_vtx.Perp();
    //    LamBr.cosa2d  = v3_fit.XYvector()*bs_to_vtx.XYvector()/v3_fit.Perp()/bs_to_vtx.Perp();
        
    //    LamBr.ppt     = cand_p->pt();
    //    LamBr.pipt    = cand_pi->pt();
        
    //    LamTree->Fill();
    //}
    
    //// build Lambda_b
    
    //BPHRecoBuilder bLamb( es );
    //bLamb.setMinPDiffererence( 1.0e-5 );
    //bLamb.add( "JPsi", lJPsi );
    //bLamb.add( "Lambda",  lLam );
    //bLamb.filter( "JPsi", MassSelect( 3.0969-0.15, 3.0969+0.15 ) );
    //bLamb.filter( "Lambda",  MassSelect( 1.115683-0.03, 1.115683+0.03 ) );
    //bLamb.filter( MassSelect(4.6, 6.7) );
    //vector<BPHRecoConstCandPtr> lLamb = BPHRecoCandidate::build( bLamb );

    //// apply kinematic fit, fill the tree
    //for ( unsigned int iLamb = 0; iLamb < lLamb.size(); ++iLamb ) {
    //    const BPHRecoCandidate* cand = lLamb[iLamb].get();
    //    const BPHRecoCandidate* cand_jpsi = cand->getComp( "JPsi" ).get();
    //    const BPHRecoCandidate* cand_lam = cand->getComp( "Lambda" ).get();
        
    //    // --------------------------------------------------
    //    // perform vertex fit for Lambda
    //    vector<RefCountedKinematicParticle> kComp_lam = cand_lam->kinParticles();
        
    //    KinematicParticleVertexFitter vtxFitter;
    //    RefCountedKinematicTree compTree_lam = vtxFitter.fit( kComp_lam );
    //    if (compTree_lam->isEmpty()) continue;
    //    compTree_lam->movePointerToTheTop();
    //    const RefCountedKinematicParticle kPart_lam = compTree_lam->currentParticle();
    //    const RefCountedKinematicVertex kVtx_lam = compTree_lam->currentDecayVertex();
    //    const KinematicState& kState_lam = kPart_lam->currentState();
    //    if (!kState_lam.isValid()) continue;
        
    //    double vtxProb_lam = TMath::Prob( kVtx_lam->chiSquared(), kVtx_lam->degreesOfFreedom() );
    //    if (vtxProb_lam <=0.) continue; // Lambda fit prob > 0.0

    //    // --------------------------------------------------
    //    // perform J/psi+Lambda vertex mass constrained fit
    //    vector<RefCountedKinematicParticle> kComp = cand_jpsi->kinParticles(); // mu+, mu-
        
    //    VirtualKinematicParticleFactory vFactory;
    //    float lam_chi = kVtx_lam->chiSquared();
    //    float lam_ndf = kVtx_lam->degreesOfFreedom();
    //    kComp.push_back(vFactory.particle(kState_lam,lam_chi,lam_ndf,kPart_lam)); // Lambda

    //    ParticleMass jpsi_mass = 3.096916;
    //    TwoTrackMassKinematicConstraint *jpsi_const = new TwoTrackMassKinematicConstraint(jpsi_mass);
    //    KinematicConstrainedVertexFitter kcvFitter;
    //    RefCountedKinematicTree compTree = kcvFitter.fit(kComp, jpsi_const);
    //    if (compTree->isEmpty()) continue;
    //    compTree->movePointerToTheTop();
    //    const RefCountedKinematicParticle kPart = compTree->currentParticle();
    //    const RefCountedKinematicVertex kVtx = compTree->currentDecayVertex();
    //    const KinematicState& kState = kPart->currentState();
    //    if (!kState.isValid()) continue;
        
    //    double vtxProb = TMath::Prob( kVtx->chiSquared(), kVtx->degreesOfFreedom() );
    //    if (vtxProb <=0.) continue; // Lambda_b fit prob > 0.0
        
    //    // --------------------------------------------------
    //    // fill tree
    //    TVector3 v3_fit(kState.kinematicParameters().momentum().x(),
    //                    kState.kinematicParameters().momentum().y(),
    //                    kState.kinematicParameters().momentum().z());
        
    //    TVector3 v3_fit_lam(kState_lam.kinematicParameters().momentum().x(),
    //                        kState_lam.kinematicParameters().momentum().y(),
    //                        kState_lam.kinematicParameters().momentum().z());
        
    //    TVector3 bs_to_vtx(kVtx->position().x()-BeamSpotVtx.position().X(),
    //                       kVtx->position().y()-BeamSpotVtx.position().Y(),
    //                       kVtx->position().z()-BeamSpotVtx.position().Z());
        
    //    TVector3 bs_to_lamvtx(kVtx_lam->position().x()-BeamSpotVtx.position().X(),
    //                          kVtx_lam->position().y()-BeamSpotVtx.position().Y(),
    //                          kVtx_lam->position().z()-BeamSpotVtx.position().Z());
        
    //    TVector3 vtx_to_lamvtx(kVtx_lam->position().x()-kVtx->position().x(),
    //                           kVtx_lam->position().y()-kVtx->position().y(),
    //                           kVtx_lam->position().z()-kVtx->position().z());

    //    LambBr.mass          = cand->composite().mass();
    //    LambBr.pt            = cand->composite().pt();
    //    LambBr.eta           = cand->composite().eta();
    //    LambBr.phi           = cand->composite().phi();
        
    //    LambBr.vtxmass       = kState.mass();
    //    LambBr.vtxprob       = vtxProb;
    //    LambBr.fl2d          = bs_to_vtx.Perp();
    //    LambBr.cosa2d        = v3_fit.XYvector()*bs_to_vtx.XYvector()/v3_fit.Perp()/bs_to_vtx.Perp();
        
    //    LambBr.jpsimass      = cand_jpsi->composite().mass();
    //    LambBr.jpsipt        = cand_jpsi->composite().pt();
    //    LambBr.jpsivtxprob   = TMath::Prob( cand_jpsi->vertex().chi2(), cand_jpsi->vertex().ndof() );
        
    //    LambBr.lammass       = cand_lam->composite().mass();
    //    LambBr.lampt         = cand_lam->composite().pt();
        
    //    LambBr.lamvtxmass    = kState_lam.mass();
    //    LambBr.lamvtxprob    = vtxProb_lam;
    //    LambBr.lamfl2d       = bs_to_lamvtx.Perp();
    //    LambBr.lamcosa2d     = v3_fit_lam.XYvector()*bs_to_lamvtx.XYvector()/v3_fit_lam.Perp()/bs_to_lamvtx.Perp();
    //    LambBr.lamfl2dbvtx   = vtx_to_lamvtx.Perp();
    //    LambBr.lamcosa2dbvtx = v3_fit_lam.XYvector()*vtx_to_lamvtx.XYvector()/v3_fit_lam.Perp()/vtx_to_lamvtx.Perp();
        
    //    LambTree->Fill();
    //}
}

void BPHRecoBToKstMuMu::beginRun(edm::Run const &iRun, edm::EventSetup const &iEvent){
    assistant.updateHltConfig(iRun, iEvent);
}

void BPHRecoBToKstMuMu::endJob() {
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( BPHRecoBToKstMuMu );
