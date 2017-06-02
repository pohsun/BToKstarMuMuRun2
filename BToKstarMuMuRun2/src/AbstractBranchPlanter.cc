#include "BphAna/BToKstarMuMuRun2/interface/AbstractBranchPlanter.h"
#include "TString.h"

AbstractBranchPlanter::AbstractBranchPlanter(){
}

AbstractBranchPlanter::~AbstractBranchPlanter(){
    // Delete histogram
    for(std::map<std::string,TH1D*>::iterator it=histMap_.begin(); it!= histMap_.end(); ++it){
        delete it->second;
        it->second = 0;
        histMap_.erase(it);
    }
}

bool AbstractBranchPlanter::isHandlesOk(){
    // Check if all needed materials are ready.
    // If not, raise ERROR.
    return true;
}

void AbstractBranchPlanter::cuttage(TTree* tree, const char prefix[]){
    // Do this in "beginJob" to register branches to a given tree.
    // Remark: It's not always necessary to keep a tree.
    //         Sometimes, we use only 'plant' funciton to fill content.
    tree_ = tree;
    regBranches(prefix);
    return;
}

void AbstractBranchPlanter::plant(const edm::Event& ev, const edm::EventSetup& es){
    // The value of all contents is decided here!
    // Don't put TTree::Fill in case that branch sets in a shared tree.
    // This class handles branches, NOT the tree itself
    runNo = ev.id().run();
    evtNo = ev.id().event();
    clearVariables();
    buildBranches(ev,es);
    return;
}

void AbstractBranchPlanter::reap(const char prefix[]){
    // Do this in  "endJob".
    if ( tree_ != 0 ){
        tree_->BuildIndex(TString::Format("%s.runNo",prefix).Data(),TString::Format("%s.evtNo",prefix).Data());
    }
    for (std::map<std::string,TH1D*>::iterator it=histMap_.begin(); it!=histMap_.end(); ++it){
        it->second->Write();
    }
    return;
}

void AbstractBranchPlanter::regBranches(const char prefix[]){
    // Remark: The runNo and evtNo must be kept in all trees for indexing.
    //         See TTree::BuildIndex for more information.
    tree_->Branch(TString::Format("%s.runNo",prefix).Data(),&runNo);
    tree_->Branch(TString::Format("%s.evtNo",prefix).Data(),&evtNo);
    return;
}

void AbstractBranchPlanter::linkAssistant(CommonBDecayAssistant *ptr){
    assistant_ = ptr;
}

void AbstractBranchPlanter::print() const{
    //getClassName();
    //getMotherClassName();
    return;
}

bool AbstractBranchPlanter::setHistFile(TFile *file){
    if (file != 0){
        file_ = file;
        return true;
    }
    return false;
}

bool AbstractBranchPlanter::setBranchAddress(TTree *tree, const char prefix[]){
    if (tree != 0){
        tree->SetBranchAddress(TString::Format("%s.runNo",prefix).Data(),&runNo);
        tree->SetBranchAddress(TString::Format("%s.evtNo",prefix).Data(),&evtNo);
        return true;
    }
    return false;
}

TH1D* AbstractBranchPlanter::getHist(std::string name) const{
    return histMap_.find(name)->second;
}

std::vector<std::string> AbstractBranchPlanter::getHistList() const{
    std::vector<std::string> output;
    for (std::map<std::string,TH1D*>::const_iterator it=histMap_.begin(); it!=histMap_.end(); ++it){
        output.push_back(it->first);
    }
    return output;
}

void AbstractBranchPlanter::getEntryWithIndex(unsigned int run, unsigned int evt){
    // Ref: https://root.cern.ch/root/html/tutorials/tree/treefriend.C.html
    if (run != runNo || evt != evtNo){
        tree_->GetEntryWithIndex(run,evt);
    }else{
    }
    return;
}

bool AbstractBranchPlanter::addHist(const TH1D &h1d){
    std::string name(h1d.GetName());
    if (histMap_.find(name) == histMap_.end()){
        // TODO: Check duplication.
        histMap_[name] = new TH1D(h1d);
        return true;
    }
    return false;
}

bool AbstractBranchPlanter::addHist(const char* name, const char* title, unsigned int nBins, double x_min, double x_max){
    if (histMap_.find(name) == histMap_.end()){
        histMap_[name] = new TH1D(name, title, nBins, x_min, x_max);
        return true;
    }
    return false;
}

bool AbstractBranchPlanter::addHist(HistArgs &args){
    return addHist(args.name, args.title, args.nBins, args.x_min, args.x_max);
}

