#include "BphAna/BToKstarMuMuRun2/interface/EvtBranchPlanter.h"
#include "TString.h"

HistArgs evtHists[]={
    {"nEvents","",1,0,1}
};

EvtBranchPlanter::EvtBranchPlanter(){
    triggeredNames = new std::vector<std::string>;
    triggeredPrescales = new std::vector<int>;
}

EvtBranchPlanter::~EvtBranchPlanter(){
    delete triggeredNames;
    delete triggeredPrescales;
}

bool EvtBranchPlanter::setBranchAddress(TTree* tree, const char prefix[]){
    if (AbstractBranchPlanter::setBranchAddress(tree,prefix)){
        tree->SetBranchAddress(TString::Format("%s.lumiBlock"         ,prefix).Data() , &lumiBlock);
        tree->SetBranchAddress(TString::Format("%s.nPV"               ,prefix).Data() , &nPV);
        tree->SetBranchAddress(TString::Format("%s.triggeredNames"    ,prefix).Data() , &triggeredNames);
        tree->SetBranchAddress(TString::Format("%s.triggeredPrescales",prefix).Data() , &triggeredPrescales);
        return true;
    }else{
        return false;
    }
}

void EvtBranchPlanter::print(){
    return;
}

bool EvtBranchPlanter::isHandlesOk(){
    if ( !AbstractBranchPlanter::isHandlesOk()) {return false;}
    if ( !assistant_->hltHandle.isValid()) {return false;}
    if ( !assistant_->pvHandle.isValid() ) {return false;}
    return true;
}

void EvtBranchPlanter::regBranches(const char prefix[]){
    AbstractBranchPlanter::regBranches(prefix);
    tree_->Branch(TString::Format("%s.lumiBlock"       ,prefix).Data() , &lumiBlock);
    tree_->Branch(TString::Format("%s.nPV"             ,prefix).Data() , &nPV);
    tree_->Branch(TString::Format("%s.triggeredNames"    ,prefix).Data() , &triggeredNames);
    tree_->Branch(TString::Format("%s.triggeredPrescales",prefix).Data() , &triggeredPrescales);
    return;
}

void EvtBranchPlanter::clearVariables(){
    lumiBlock  = 0;
    nPV        = 0;
    triggeredNames       ->clear();
    triggeredPrescales   ->clear();
    return;
}

void EvtBranchPlanter::hltReport(const edm::Event& ev){
    if (assistant_->hltHandle.isValid()) {
        const edm::TriggerNames& evTriggerNames_ = ev.triggerNames(*(assistant_->hltHandle));

        for (unsigned int iTrg = 0; iTrg < assistant_->hltHandle->size(); iTrg++){

            // Only consider the triggered case.
            if (assistant_->hltHandle->at(iTrg).accept() == 1){

                std::string triggeredname = evTriggerNames_.triggerName(iTrg);
                int triggeredprescale = assistant_->hltConfigPtr_->prescaleValue(iTrg, triggeredname);

                // Loop over our interested HLT trigger names to find if this event contains.
                for(unsigned int it=0; it < assistant_->triggerNames.size(); it++){
                    if (triggeredname.find(assistant_->triggerNames[it]) != std::string::npos){
                        // save the name without version
                        triggeredNames      ->push_back(assistant_->triggerNames[it]);
                        triggeredPrescales  ->push_back(triggeredprescale);
                    }
                }
            }
        }
    }
    return;
}

void EvtBranchPlanter::buildBranches(const edm::Event& ev, const edm::EventSetup& es){
    lumiBlock = ev.luminosityBlock();
    nPV = assistant_->pvHandle->size();
    hltReport(ev);
    return;
}

