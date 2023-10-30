#include "analysis_suite/skim/interface/tth_mva.h"

#include "analysis_suite/skim/interface/Electron.h"
#include "analysis_suite/skim/interface/Jet.h"

void TTH_MVA::setup(std::string xml_file, Electron* electron) {
    reader = new TMVA::Reader("Silent");

    ele = electron;
    reader->AddVariable("Electron_pt", &pt);
    reader->AddVariable("Electron_eta", &eta);
    reader->AddVariable("Electron_pfRelIso03_all", &pfRelIso03_all);
    reader->AddVariable("Electron_miniPFRelIso_chg", &miniPFRelIso_chg);
    reader->AddVariable("Electron_miniRelIsoNeutral := Electron_miniPFRelIso_all - Electron_miniPFRelIso_chg", &miniRelIsoNeutral);
    reader->AddVariable("Electron_jetNDauCharged", &jetNDauCharged);
    reader->AddVariable("Electron_jetPtRelv2", &jetPtRelv2);
    reader->AddVariable("Electron_jetPtRatio := min(1 / (1 + Electron_jetRelIso), 1.5)", &jetPtRatio);
    reader->AddVariable("Electron_jetBTagDeepFlavB := Electron_jetIdx > -1 ? Jet_btagDeepFlavB[Electron_jetIdx] : 0", &jetBTagDeepFlavB);
    reader->AddVariable("Electron_sip3d", &sip3d);
    reader->AddVariable("Electron_dxy := log(abs(Electron_dxy))", &log_dxy);
    reader->AddVariable("Electron_dz  := log(abs(Electron_dz))", &log_dz);
    reader->AddVariable("Electron_mvaFall17V2noIso", &mvaFall17V2noIso);

    reader->BookMVA("BDTG", xml_file);
}

float TTH_MVA::get_mva(size_t i, Jet& jet) {
    pt = ele->pt(i);
    eta = ele->eta(i);
    pfRelIso03_all = ele->ptRel.at(i);
    miniPFRelIso_chg = ele->iso_chg.at(i);
    miniRelIsoNeutral = ele->iso.at(i) - ele->iso_chg.at(i);
    jetNDauCharged = ele->jetNDauCharged.at(i);
    jetPtRelv2 = ele->ptRel.at(i);
    jetPtRatio = std::min(ele->ptRatio(i), (float)1.5);
    jetBTagDeepFlavB = ele->jetIdx.at(i) > -1 ? jet.btag.at(ele->jetIdx.at(i)) : 0;
    sip3d = ele->sip3d.at(i);
    log_dxy = std::log(fabs(std::max(min_val, ele->dxy.at(i))));
    log_dz = std::log(fabs(std::max(min_val, ele->dz.at(i))));
    mvaFall17V2noIso = ele->mva_vals.at(i);

    return reader->EvaluateMVA("BDTG");
}
