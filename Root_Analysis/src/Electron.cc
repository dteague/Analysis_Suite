#include "analysis_suite/Root_Analysis/interface/Electron.h"
#include "analysis_suite/Root_Analysis/interface/Jet.h"

void Electron::setup(TTreeReader& fReader, int year) {
    Lepton::setup("Electron", fReader, year);
    eCorr = new TTRArray<Float_t>(fReader, "Electron_eCorr");
    lostHits = new  TTRArray<UChar_t>(fReader, "Electron_lostHits");
    convVeto = new  TTRArray<Bool_t>(fReader, "Electron_convVeto");
    iso = new  TTRArray<Float_t>(fReader, "Electron_miniPFRelIso_all");
    dz = new  TTRArray<Float_t>(fReader, "Electron_dz");
    dxy = new  TTRArray<Float_t>(fReader, "Electron_dxy");
    sieie = new  TTRArray<Float_t>(fReader, "Electron_sieie");
    hoe = new  TTRArray<Float_t>(fReader, "Electron_hoe");
    eInvMinusPInv = new  TTRArray<Float_t>(fReader, "Electron_eInvMinusPInv");
    sip3d = new  TTRArray<Float_t>(fReader, "Electron_sip3d");
    tightCharge = new  TTRArray<Int_t>(fReader, "Electron_tightCharge");
    ecalSumEt = new TTRArray<Float_t>(fReader, "Electron_dr03EcalRecHitSumEt");
    hcalSumEt = new  TTRArray<Float_t>(fReader, "Electron_dr03HcalDepth1TowerSumEt");
    tkSumPt = new  TTRArray<Float_t>(fReader, "Electron_dr03TkSumPt");

    ptRatioCut = 0.78;
    ptRelCut = pow(8.0, 2);
        
    if(year_ == yr2016) {
        ptRatioCut = 0.8;
        ptRelCut = pow(7.2, 2);
        mva = new  TTRArray<Float_t>(fReader, "Electron_mvaSpring16GP");
        mvaLoose = {{-0.46, -0.48, 0, -0.85},
                    {-0.03, -0.67, 0, -0.91},
                    {0.06, -0.49, 0, -0.83}};
        mvaTight = {{10, 0.77, 0, 0.52},
                    {10, 0.56, 0, 0.11},
                    {10, 0.48, 0, -0.01}};
        /// Setup interpolation used for 25>pt>15
        for(size_t i = 0; i < mvaLoose.size(); i++) {
            mvaLoose[i][2] = (mvaLoose[i][3] - mvaLoose[i][1])/10;
            mvaTight[i][2] = (mvaTight[i][3] - mvaTight[i][1])/10;
        }
    } else if(year_ == yr2017) {
        mva = new  TTRArray<Float_t>(fReader, "Electron_mvaFall17V1noIso");
        mvaLoose = {{0.488, -0.738667, 0.00986667, -0.64},
                    {-0.045, -0.825, 0.005, -0.775},
                    {0.176, -0.784333, 0.00513333, -0.733}};
        mvaTight = {{10, 0.36, 0.032, 0.68},
                    {10, 0.225, 0.025, 0.475},
                    {10, 0.04, 0.028, 0.32}};
    } else if(year_ == yr2018) {
        mva = new  TTRArray<Float_t>(fReader, "Electron_mvaFall17V2noIso");
        mvaLoose = {{1.32, 0.544, 0.066, 1.204},
                    {0.192, -0.246, 0.033, 0.084},
                    {0.362, -0.653, 0.053, -0.123}};
        mvaTight = {{100, 3.157, 0.112, 4.277},
                    {100, 2.552, 0.06, 3.152},
                    {100, 1.489, 0.087, 2.359}};
    }

}

void Electron::createLooseList() {
    for (size_t i = 0; i < pt->GetSize(); i++) {
        if (pt->At(i)/eCorr->At(i) > 7 &&
            abs(eta->At(i)) < 2.5 &&
            convVeto->At(i) &&
            lostHits->At(i) <= 1 &&
            abs(dz->At(i)) < 0.1 &&
            abs(dxy->At(i)) < 0.05 &&
            passMVACut(i, false) &&
            iso->At(i) < 0.4 && 
            abs(eInvMinusPInv->At(i) < 0.01) &&
            hoe->At(i) < 0.08 &&
            ((abs(eta->At(i)) < BARREL_ETA && sieie->At(i) < 0.011)  ||
             (abs(eta->At(i)) >= BARREL_ETA && sieie->At(i) < 0.031))
        )
            looseList.push_back(i);
    }
}

/// closejet w iso
void Electron::createFakeList(Particle& jets) {
    std::vector<size_t> pre_list;
    for (auto i: looseList) {
        if (pt->At(i)/eCorr->At(i) >= 10 &&
            sip3d->At(i) < 4 &&
            lostHits->At(i) == 0 &&
            tightCharge->At(i) == 2
        )
            pre_list.push_back(i);
    }
    for(auto i: pre_list) {
        auto closeJetInfo = getCloseJet(i, jets);
        if (passJetIsolation(i, jets))
            fakeList.push_back(i);
        dynamic_cast<Jet&>(jets).closeJetDr_by_index.insert(closeJetInfo);
    }
}

// need iso
void Electron::createTightList() {
    for (auto i: fakeList) {
        float pt_cor = pt->At(i)/eCorr->At(i);
        if (pt_cor > 15 &&
            iso->At(i) < 0.12 &&
            ecalSumEt->At(i) / pt_cor < 0.45 &&
            hcalSumEt->At(i) / pt_cor < 0.25 &&
            tkSumPt->At(i) / pt_cor < 0.2 &&
            passMVACut(i, true) 
        )
            tightList.push_back(i);
    }
}

bool Electron::passMVACut(size_t idx, bool isTight) {
    int caseIdx = 0;
    float pt_cor = pt->At(idx)/eCorr->At(idx);

    //// ETA Splitting
    if(pt_cor < 5)       return false;
    else if(abs(eta->At(idx)) < 0.8)        caseIdx = 0;
    else if(abs(eta->At(idx)) < BARREL_ETA) caseIdx = 1;
    else  caseIdx = 2;
    auto& mvaCuts = (isTight) ? mvaTight[caseIdx] : mvaLoose[caseIdx];

    //// PT Splitting
    if(pt_cor < 10) caseIdx = 0;
    else if(pt_cor < 15 && year_ == yr2016) caseIdx = 1;
    else if(pt_cor < 25) caseIdx = 2;
    else caseIdx = 3;

    double mvaValue = mva->At(idx);
    if(year_ == yr2018) mvaValue = atanh(mva->At(idx));

    if(caseIdx % 4 != 2) return mvaValue > mvaCuts[caseIdx%4];
    else                 return mvaValue > mvaCuts[1] + mvaCuts[2]*(pt_cor-15);
}
