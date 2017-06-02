#include "BphAna/BToKstarMuMuRun2/interface/MuonBranchPlanter.h"
#include "TString.h"

HistArgs muonHists[]={
    {"nMuons","",1,0,1},
    {"muonPt","",1,0,1},
    {"muonEta","",1,0,1},
    {"muonPhi","",1,0,1}
};

MuonBranchPlanter::MuonBranchPlanter(){
}

MuonBranchPlanter::~MuonBranchPlanter(){
}


bool MuonBranchPlanter::setBranchAddress(TTree* tree, const char prefix[]){
    if (AbstractBranchPlanter::setBranchAddress(tree,prefix)){
        tree_->SetBranchAddress(TString::Format("%s.nMu"       , prefix).Data() , &nMu);
        tree_->SetBranchAddress(TString::Format("%s.charge"    , prefix).Data() , &charge);
        tree_->SetBranchAddress(TString::Format("%s.pt"        , prefix).Data() , &pt);
        tree_->SetBranchAddress(TString::Format("%s.eta"       , prefix).Data() , &eta);
        tree_->SetBranchAddress(TString::Format("%s.phi"       , prefix).Data() , &phi);
        tree_->SetBranchAddress(TString::Format("%s.nPixHit"   , prefix).Data() , &nPixHit);
        tree_->SetBranchAddress(TString::Format("%s.nTrkHit"   , prefix).Data() , &nTrkHit);
        tree_->SetBranchAddress(TString::Format("%s.nPixLayer" , prefix).Data() , &nPixLayer);
        tree_->SetBranchAddress(TString::Format("%s.nTrkLayer" , prefix).Data() , &nTrkLayer);
        tree_->SetBranchAddress(TString::Format("%s.dxy"       , prefix).Data() , &dxy);
        tree_->SetBranchAddress(TString::Format("%s.dz"        , prefix).Data() , &dz);
        tree_->SetBranchAddress(TString::Format("%s.normChi2"  , prefix).Data() , &normChi2);
        tree_->SetBranchAddress(TString::Format("%s.dcaBS"     , prefix).Data() , &dcaBS);
        tree_->SetBranchAddress(TString::Format("%s.dcaBSErr"  , prefix).Data() , &dcaBSErr);
        //pt,eta,phi,normChi2,isL1Matched,isMCMatched
        return true;
    }else{
        return false;
    }
}

void MuonBranchPlanter::print(){
    return;
}

bool MuonBranchPlanter::isHandlesOk(){
    return true;
}

void MuonBranchPlanter::regBranches(const char prefix[]){
    AbstractBranchPlanter::regBranches(prefix);
    tree_->Branch(TString::Format("%s.nMu"       , prefix).Data() , &nMu);
    tree_->Branch(TString::Format("%s.charge"    , prefix).Data() , &charge);
    tree_->Branch(TString::Format("%s.pt"        , prefix).Data() , &pt);
    tree_->Branch(TString::Format("%s.eta"       , prefix).Data() , &eta);
    tree_->Branch(TString::Format("%s.phi"       , prefix).Data() , &phi);
    tree_->Branch(TString::Format("%s.nPixHit"   , prefix).Data() , &nPixHit);
    tree_->Branch(TString::Format("%s.nTrkHit"   , prefix).Data() , &nTrkHit);
    tree_->Branch(TString::Format("%s.nPixLayer" , prefix).Data() , &nPixLayer);
    tree_->Branch(TString::Format("%s.nTrkLayer" , prefix).Data() , &nTrkLayer);
    tree_->Branch(TString::Format("%s.dxy"       , prefix).Data() , &dxy);
    tree_->Branch(TString::Format("%s.dz"        , prefix).Data() , &dz);
    tree_->Branch(TString::Format("%s.normChi2"  , prefix).Data() , &normChi2);
    tree_->Branch(TString::Format("%s.dcaBS"     , prefix).Data() , &dcaBS);
    tree_->Branch(TString::Format("%s.dcaBSErr"  , prefix).Data() , &dcaBSErr);
    return;
}

void MuonBranchPlanter::clearVariables(){
    nMu=0;
    charge   .clear();
    pt       .clear();
    eta      .clear();
    phi      .clear();
    nPixHit  .clear();
    nTrkHit  .clear();
    nPixLayer.clear();
    nTrkLayer.clear();
    dxy      .clear();
    dz       .clear();
    normChi2 .clear();
    dcaBS    .clear();
    dcaBSErr .clear();
    return;
}

void MuonBranchPlanter::buildBranches(const edm::Event& ev, const edm::EventSetup& es){
    if ( !isHandlesOk() ){
        return;
    }
    std::vector<const pat::Muon*>* cands = assistant_->getMu();
    reco::TrackRef innerTrack;
    reco::Vertex*  primaryVertex = assistant_->getPV();
    double l_dcaBS, l_dcaBSErr;
    for(unsigned int iMu = 0; iMu < cands->size(); iMu++){

        nMu += 1;
        charge      .push_back(cands->at(iMu)->charge());
        pt          .push_back(cands->at(iMu)->pt());
        eta         .push_back(cands->at(iMu)->eta());
        phi         .push_back(cands->at(iMu)->phi());

        innerTrack = cands->at(iMu)->innerTrack();
        nPixHit     .push_back(innerTrack->hitPattern().numberOfValidPixelHits());
        nPixLayer   .push_back(innerTrack->hitPattern().pixelLayersWithMeasurement());
        nTrkHit     .push_back(innerTrack->hitPattern().numberOfValidTrackerHits());
        nTrkLayer   .push_back(innerTrack->hitPattern().trackerLayersWithMeasurement());
        dxy         .push_back(innerTrack->dxy(primaryVertex->position()));
        dz          .push_back(innerTrack->dz(primaryVertex->position()));
        normChi2    .push_back(innerTrack->normalizedChi2());
        
        calcDcaBS(innerTrack,l_dcaBS,l_dcaBSErr);
        dcaBS       .push_back(l_dcaBS);
        dcaBSErr    .push_back(l_dcaBSErr);

    }
    return;
}

void MuonBranchPlanter::calcDcaBS(reco::TrackRef trk, double& l_dcaBS, double& l_dcaBSErr){
    l_dcaBS     = -1;
    l_dcaBSErr  = -1;
    assistant_->getBS();

    const reco::TransientTrack trkT(trk, &*(assistant_->bFieldHandle));
    TrajectoryStateClosestToPoint theDcaBS = trkT.trajectoryStateClosestToPoint(
            GlobalPoint(assistant_->bsPtr_->position().x(),
                        assistant_->bsPtr_->position().y(),
                        assistant_->bsPtr_->position().z()));
    if ( ! theDcaBS.isValid() ){
        return;
    }

    l_dcaBS = theDcaBS.perigeeParameters().transverseImpactParameter();
    l_dcaBSErr = theDcaBS.perigeeError().transverseImpactParameterError();
    return;
}
