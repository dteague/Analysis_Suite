#include "analysis_suite/skim/interface/TagProbe.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

namespace Channel {
    enum {
        OS,
        None,
    };
}

enum class Subchannel {
    MM,
    EE,
    None,
};

void TagProbe::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    // Charge Mis-id Fake Rate
    createTree("OS", Channel::OS);
    final_tree = trees.at(Channel::OS).tree;
    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    trigobj.setup(fReader);

    if (dataset_ == Dataset::Single_M)
        file_chan = Subchannel::MM;
    else if (dataset_ == Dataset::Single_E || dataset_ == Dataset::DoubleEG)
        file_chan = Subchannel::EE;
    elec.mvaCut = -2;
    muon.mvaCut = -2;

    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) {
        if (data_run != "H") {
            trig_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");
        } else {
            trig_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
        }
        trig_ee.setup(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    } else if (year_ == Year::yr2017) {
        if (data_run == "B") {
            trig_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
        } else {
            trig_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
        }
        trig_ee.setup(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
    } else if (year_ == Year::yr2018) {
        trig_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
        trig_ee.setup(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    }
    LOG_FUNC << "End of Init";
}

void TagProbe::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    tree->Branch("tag_pt", &tag_pt);
    tree->Branch("tag_abseta", &tag_abseta);
    tree->Branch("tag_iso", &tag_iso);
    tree->Branch("tag_tight", &tag_tight);

    tree->Branch("probe_pt", &probe_pt);
    tree->Branch("probe_fakept", &probe_fakept);
    tree->Branch("probe_abseta", &probe_abseta);
    tree->Branch("probe_iso", &probe_iso);
    tree->Branch("probe_tight", &probe_tight);
    tree->Branch("probe_mva", &probe_mva);

    tree->Branch("mass", &mass);
    tree->Branch("mass_fake", &mass_fake);
    tree->Branch("njets", &njets);

    tree->Branch("isMuon", &isMuon);
    tree->Branch("passEE", &pass_ee);
    tree->Branch("passMM", &pass_mm);

    tree->Branch("wgt", &wgt);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void TagProbe::setOtherGoodParticles(size_t syst)
{
    if (muon.size(Level::Loose) == 2 || elec.size(Level::Loose) == 2) {
        trigobj.setupGoodLists(elec, muon);

        //// setup probe pair
        for (auto lidx : muon.list(Level::Loose)) {
            for (auto fidx : muon.list(Level::Fake)) {
                if (lidx != fidx) muon_pair[lidx] = fidx;
            }
        }
        for (auto lidx : elec.list(Level::Loose)) {
            for (auto fidx : elec.list(Level::Fake)) {
                if (lidx != fidx) elec_pair[lidx] = fidx;
            }
        }
    }

}

void TagProbe::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    trigobj.clear();
    muon_pair.clear();
    elec_pair.clear();
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void TagProbe::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    // tag_pt.clear();
    // tag_abseta.clear();
    // tag_iso.clear();
    // probe_pt.clear();
    // probe_fakept.clear();
    // probe_abseta.clear();
    // probe_iso.clear();
    // probe_mva.clear();
    // mass_fake.clear();
    // mass.clear();
    // njets.clear();
    // isMuon.clear();
    // wgt.clear();
    LOG_FUNC << "End of clearOutputs";
}

void TagProbe::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= jet.getScaleFactor();
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

bool TagProbe::isOppositeSign(Level level)
{
    int q_total = 0;
    for (size_t idx : muon.list(level)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(level)) {
        q_total += elec.charge(idx);
    }
    return q_total == 0;
}


void TagProbe::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Loose) == 2) {
        if (muon.size(Level::Loose) == 2) {
            subChannel_ = Subchannel::MM;
        } else if (elec.size(Level::Loose) == 2) {
            subChannel_ = Subchannel::EE;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}



bool TagProbe::getCutFlow()
{
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (closure_cuts()) {
        (*currentChannel_) = Channel::OS;
        FillValues(0);
    }

    return false;
}

bool TagProbe::closure_cuts() {
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2Lepton", nLeps(Level::Loose) == 2);
    passCuts &= cuts.setCut("passSameFlavor", muon.size(Level::Loose)*elec.size(Level::Loose) == 0);
    passCuts &= cuts.setCut("passOppositeSign", isOppositeSign(Level::Loose));
    passCuts &= cuts.setCut("passMatchedTrigObj", trigobj.is_match());
    passCuts &= cuts.setCut("passDataset", isMC_ || file_chan == subChannel_);

    fillCutFlow(Channel::OS, cuts);

    return passCuts;
}

float TagProbe::get_mass(size_t t, size_t p, bool useRaw) {
    if (subChannel_ == Subchannel::MM) {
        auto tag = muon.p4(t);
        tag.SetPt(muon.rawpt(t));
        auto probe = muon.p4(p);
        if (useRaw) probe.SetPt(muon.rawpt(p));
        return (tag+probe).M();
    } else if (subChannel_ == Subchannel::EE) {
        auto tag = elec.p4(t);
        tag.SetPt(elec.rawpt(t));
        auto probe = elec.p4(p);
        if (useRaw) probe.SetPt(elec.rawpt(p));
        return (tag+probe).M();
    }
    return -1;
}

void TagProbe::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    setupSyst(0);

    if (subChannel_ == Subchannel::MM) {
        for (auto tag_i : trigobj.muon_match) {
            if (muon_pair.find(tag_i) == muon_pair.end()) continue;
            size_t probe_i = muon_pair[tag_i];
            isMuon = true;
            tag_pt = muon.rawpt(tag_i);
            tag_abseta = fabs(muon.eta(tag_i));
            tag_iso = muon.iso.at(tag_i);
            tag_tight = muon.tid.at(tag_i);
            probe_pt = muon.rawpt(probe_i);
            probe_fakept = muon.pt(probe_i);
            probe_abseta = fabs(muon.eta(probe_i));
            probe_iso = muon.iso.at(probe_i);
            probe_mva = muon.mvaTTH.at(probe_i);
            probe_tight = muon.tid.at(probe_i);
            mass_fake = get_mass(tag_i, probe_i, false);
            mass = get_mass(tag_i, probe_i, true);
            njets = static_cast<int>(jet.size(Level::Loose));
            wgt = *weight;
            pass_mm = *trig_mm;
            pass_ee = *trig_ee;
            final_tree->Fill();
        }
    } else if (subChannel_ == Subchannel::EE) {
        for (auto tag_i : trigobj.elec_match) {
            if (elec_pair.find(tag_i) == elec_pair.end()) continue;
            size_t probe_i = elec_pair[tag_i];
            isMuon = false;
            tag_pt = elec.rawpt(tag_i);
            tag_abseta = fabs(elec.eta(tag_i));
            tag_iso = elec.iso.at(tag_i);
            tag_tight = elec.mva_80.at(tag_i);
            probe_pt = elec.rawpt(probe_i);
            probe_fakept = elec.pt(probe_i);
            probe_abseta = fabs(elec.eta(probe_i));
            probe_iso = elec.iso.at(probe_i);
            probe_mva = elec.mvaTTH.at(probe_i);
            probe_tight = elec.mva_80.at(probe_i);
            mass_fake = get_mass(tag_i, probe_i, false);
            mass = get_mass(tag_i, probe_i, true);
            njets = static_cast<int>(jet.size(Level::Loose));
            wgt = *weight;
            pass_mm = *trig_mm;
            pass_ee = *trig_ee;
            final_tree->Fill();
        }
    }

    LOG_FUNC << "End of FillValues";
}

