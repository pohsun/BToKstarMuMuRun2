#include "BphAna/BToKstarMuMuRun2/interface/pctools.h"

using namespace pctools;

PMessanger::PMessanger(std::string fname)
{
    ofName_ = fname;
    ofPtr_ = fopen(ofName_.c_str(), "w");
}

void PMessanger::out(std::string msg, int verb, int mode){
    if ( mode != SIO && verb > this->verbosityLevel_ ) {
        return;
    }

    std::string newMsg(verbosityTitle[verb]+msg+'\n');

    switch(mode){
        case SIO:
            if ( isRedirectStdio_ ){
                switchRedirectStdio(false);
                fprintf(stdout, newMsg.c_str());
                switchRedirectStdio(true);
            }else{
                fprintf(stdout, newMsg.c_str());
            }
            break;
        case FIO:
            fprintf(ofPtr_, newMsg.c_str());
            break;
        default: // STDIO
            fprintf(stdout, newMsg.c_str());
    }
}

void PMessanger::out(std::string msg, int verb){
    out(msg, verb, this->currentMode_);
}

void PMessanger::out(std::string msg){
    out(msg, this->verbosityLevel_, this->currentMode_);
}

void PMessanger::sprint(std::string msg, int verb){
    out(msg, verb, SIO);
}

void PMessanger::sprint(std::string msg){
    sprint(msg, this->verbosityLevel_);
}

void PMessanger::switchRedirectStdio(bool redirectStdio ){
    // This function works ONLY on unix-like system.
    // Or the stdout cannot be restored after running freopen.

    if ( redirectStdio == isRedirectStdio_ ) { return; }

    struct stat fiBuff;
    if (stat("/dev/tty",&fiBuff) != 0){
        sprint("\"/dev/tty\" is NOT available.",WARNING);
        return;
    }

    if ( redirectStdio ){
        sprint("Direct stdout/stderr to screen.",HINT);
        freopen("/dev/tty","a",stdout);
        freopen("/dev/tty","a",stderr);
        isRedirectStdio_ = false;
    }else{
        sprint(std::string("Redirect stdout/stderr to "+ofName_+"\".").c_str(),HINT);
        freopen(ofName_.c_str(),"a",stdout);
        freopen(ofName_.c_str(),"a",stderr);
        isRedirectStdio_ = true;
    }

    return;
}

void PMessanger::switchRedirectStdio(){
    switchRedirectStdio( isRedirectStdio_ ? false : true );
}

int PMessanger::setMode(int mode = -1){
    ( mode < 0 || mode >= nMODES ) ? this->currentMode_ = (this->currentMode_+1)%nMODES : this->currentMode_ = mode ;
    return this->currentMode_;
}

int PMessanger::getMode() const{
    return this->currentMode_;
}

int PMessanger::setVerbosity(int verb){
    if ( verb >= 0 && verb < nLEVELS) this->verbosityLevel_=verb;
    return this->verbosityLevel_;
}

int PMessanger::getVerbosity() const{
    return this->verbosityLevel_;
}

bool PMessanger::isRedirectStdio() const{
    return this->isRedirectStdio_;
}

FILE* PMessanger::getOfile() const{
    return this->ofPtr_;
}
