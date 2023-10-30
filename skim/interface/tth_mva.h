#ifndef TTH_MVA_H_
#define TTH_MVA_H_

class Electron;
class Jet;

#include "TMVA/Reader.h"

class TTH_MVA {
public:
    TTH_MVA() {};
    void setup(std::string xml_file, Electron* ele);
    float get_mva(size_t i, Jet& jet);

private:
    TMVA::Reader* reader;

    Electron* ele;

    Float_t pt;
    Float_t eta;
    Float_t pfRelIso03_all;
    Float_t miniPFRelIso_chg;
    Float_t miniRelIsoNeutral;
    Float_t jetNDauCharged;
    Float_t jetPtRelv2;
    Float_t jetPtRatio;
    Float_t jetBTagDeepFlavB;
    Float_t sip3d;
    Float_t log_dxy;
    Float_t log_dz;
    Float_t mvaFall17V2noIso;

    float min_val = 1e-10;
};

#endif // TTH_MVA_H_
