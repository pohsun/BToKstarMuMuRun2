#include "BphAna/BToKstarMuMuRun2/interface/BdToKstMuMuBranchPlanter.h"
#include "TString.h"

HistArgs kstmmHists[]={
    //{"nKst","",1,0,1},
    //{"kstPt" ,"",1,0,1},
    //{"kstEta","",1,0,1},
    //{"kstPhi","",1,0,1}
};

BdToKstMuMuBranchPlanter::BdToKstMuMuBranchPlanter(){
}

BdToKstMuMuBranchPlanter::~BdToKstMuMuBranchPlanter(){
}


bool BdToKstMuMuBranchPlanter::setBranchAddress(TTree* tree, const char prefix[]){
    if (AbstractBranchPlanter::setBranchAddress(tree,prefix)){
        //tree_->SetBranchAddress(TString::Format("%s.nMu"       , prefix).Data() , &nMu);
        //tree_->SetBranchAddress(TString::Format("%s.charge"    , prefix).Data() , &charge);
        return true;
    }else{
        return false;
    }
}

void BdToKstMuMuBranchPlanter::print(){
    return;
}

bool BdToKstMuMuBranchPlanter::isHandlesOk(){
    return true;
}

void BdToKstMuMuBranchPlanter::regBranches(const char prefix[]){
    AbstractBranchPlanter::regBranches(prefix);
    tree_->Branch(TString::Format("%s.nKst"       , prefix).Data() , &nKst);
    tree_->Branch(TString::Format("%s.nBd"        , prefix).Data() , &nBd);
    return;
}

void BdToKstMuMuBranchPlanter::clearVariables(){
    //nMu=0;
    //charge.clear();
    return;
}

void BdToKstMuMuBranchPlanter::buildBranches(const edm::Event& ev){
    if ( !isHandlesOk() ){
        return;
    }
    //BPHRecoBuilder bdToKstMuMuBuilder(assistant_)
    return;
}

void BdToKstMuMuBranchPlanter::calLS(
    double Vx      , double Vy           , double Vz     ,
    double Wx      , double Wy           , double Wz     ,
    double VxErr2  , double VyErr2       , double VzErr2 ,
    double VxyCov  , double VxzCov       , double VyzCov ,
    double WxErr2  , double WyErr2       , double WzErr2 ,
    double WxyCov  , double WxzCov       , double WyzCov ,
    double* deltaD , double* deltaDErr){
    *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
    if (*deltaD > 0.){
        *deltaDErr = sqrt(  (Vx-Wx) * (Vx-Wx) * VxErr2 +
                            (Vy-Wy) * (Vy-Wy) * VyErr2 +
                            (Vz-Wz) * (Vz-Wz) * VzErr2 +
        
                            (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
                            (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
                            (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
        
                            (Vx-Wx) * (Vx-Wx) * WxErr2 +
                            (Vy-Wy) * (Vy-Wy) * WyErr2 +
                            (Vz-Wz) * (Vz-Wz) * WzErr2 +
        
                            (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
                            (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
                            (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
    }else{
        *deltaDErr = 0.;
    }
}

void BdToKstMuMuBranchPlanter::calCosAlpha(
    double Vx, double Vy, double Vz,
    double Wx, double Wy, double Wz,
    double VxErr2, double VyErr2, double VzErr2,
    double VxyCov, double VxzCov, double VyzCov,
    double WxErr2, double WyErr2, double WzErr2,
    double WxyCov, double WxzCov, double WyzCov,
    double *cosAlpha, double *cosAlphaErr){

    double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
    double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
    double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

    if ((Vnorm > 0.) && (Wnorm > 0.)) {
        *cosAlpha = VdotW / (Vnorm * Wnorm);
        *cosAlphaErr = sqrt(
            ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
             (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
             (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
             (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
             (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
             (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov )/
            (Wnorm*Wnorm*Wnorm*Wnorm) +
            
            ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
             (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
             (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
             (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
             (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
             (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
            (Vnorm*Vnorm*Vnorm*Vnorm)) / 
            (Wnorm*Vnorm);
    }else{
        *cosAlpha = 0.;
        *cosAlphaErr = 0.;
    }
}

void BdToKstMuMuBranchPlanter::calCtau(
    RefCountedKinematicTree vertexFitTree,
    double &bctau, double &bctauerr){
    //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)

    vertexFitTree->movePointerToTheTop();
    RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
    RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();

    double betagamma = (b_KP->currentState().globalMomentum().mag()/BuMass_);

    // calculate ctau error. Momentum error is negligible compared to
    // the vertex errors, so don't worry about it

    GlobalPoint BVP = GlobalPoint( b_KV->position() );
    GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
        primaryVertex_.position().y(),
        primaryVertex_.position().z() );
    GlobalVector sep3D = BVP-PVP;
    GlobalVector pBV = b_KP->currentState().globalMomentum(); 
    bctau = (BuMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));

    GlobalError BVE = b_KV->error();
    GlobalError PVE = GlobalError( primaryVertex_.error() );
    VertexDistance3D theVertexDistance3D;
    Measurement1D TheMeasurement = theVertexDistance3D.distance(
        VertexState(BVP, BVE), VertexState(PVP, PVE) );
    double myError = TheMeasurement.error();  

    //  ctau is defined by the portion of the flight distance along
    //  the compoenent of the B momementum, so only consider the error
    //  of that component, too, which is accomplished by scaling by
    //  ((VB-VP)(dot)PB)/|VB-VP|*|PB| 

    double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );       
    bctauerr =  (myError*scale)/betagamma;
}

bool BdToKstMuMuBranchPlanter::hasGoodDcaBs (const reco::TransientTrack trackTT,
                double &dcaBs, double &dcaBsErr, double &maxDcaBs)
{
    TrajectoryStateClosestToPoint theDCAXBS =
        trackTT.trajectoryStateClosestToPoint(
            GlobalPoint(beamSpot_.position().x(),beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() ){  return false; }

  dcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  dcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(dcaBs) > maxDcaBs ){ return false; }
  return true;
}

bool BdToKstMuMuBranchPlanter::calClosestApproachTracks (const reco::TransientTrack trackpTT,
    const reco::TransientTrack trackmTT,
    double &trk_R,
    double &trk_Z,
    double &trk_DCA)
{
    ClosestApproachInRPhi ClosestApp;
    ClosestApp.calculate(trackpTT.initialFreeState(), trackmTT.initialFreeState());
    if (! ClosestApp.status() ){return false ;}

    GlobalPoint XingPoint = ClosestApp.crossingPoint();

    trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
    trk_Z = fabs(XingPoint.z());

    // if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) >
    //      TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_))  return false;

    trk_DCA = ClosestApp.distance();
    // if (DCAmumu > MuMuMaxDca_) return false;

    return true;
}

