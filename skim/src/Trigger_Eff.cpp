#include "analysis_suite/skim/interface/Trigger_Eff.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

namespace Channel {
    enum {
        Signal,
        None,
    };
}

void Trigger_Eff::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);
    met_type = MET_Type::PF;
    createTree("Signal", Channel::Signal);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) {
        trig_cuts.setup_channel(Subchannel::MET, Dataset::MET, fReader,
                                {"HLT_PFMET300",
                                 "HLT_MET200",
                                 "HLT_PFMET170_HBHECleaned",
                                 "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                 "HLT_PFMET120_PFMHT120_IDTight"});
        // trig_cuts.setup_channel(Subchannel::MHT, Dataset::MHT, fReader,
        //                         {"HLT_PFHT300_PFMET110"});
    } else {
        trig_cuts.setup_channel(Subchannel::MET, Dataset::MET, fReader,
                                {"HLT_PFMET200_HBHECleaned",
                                 "HLT_PFMET200_HBHE_BeamHaloCleaned",
                                 "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned",
                                 "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                 "HLT_PFMET120_PFMHT120_IDTight",
                                 "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                                 "HLT_PFMET120_PFMHT120_IDTight_PFHT60", });
        // trig_cuts.setup_channel(Subchannel::MHT, Dataset::MHT, fReader,
        //                         {"HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        //                          "HLT_PFHT700_PFMET85_PFMHT85_IDTight",
        //                          "HLT_PFHT800_PFMET75_PFMHT75_IDTight",});
    }
    LOG_FUNC << "End of Init";
}

void Trigger_Eff::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "ParticleOut", &o_tightMuons);
    tree->Branch("TightElectron", "ParticleOut", &o_tightElectrons);

    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("HT", &o_ht);
    tree->Branch("MHT", &o_mht);
    tree->Branch("Dilepton_trigger", &o_dilepton_trigger);
    // tree->Branch("Met_trigger", &o_met_trigger);

    LOG_FUNC << "End of SetupOutTreeBranches";
}

void Trigger_Eff::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_met.clear();
    o_metphi.clear();
    o_ht.clear();
    o_mht.clear();
    o_dilepton_trigger.clear();
    // o_met_trigger.clear();
    LOG_FUNC << "End of clearOutputs";
}

void Trigger_Eff::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    // (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getScaleFactor();
    // (*weight) *= jet.getTotalBTagWeight("M");
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    // (*weight) *= rTop.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void Trigger_Eff::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Tight) >= 2) {
        if (muon.size(Level::Tight) * elec.size(Level::Tight) > 0) {
            subChannel_ = Subchannel::EM;
        } else if (elec.size(Level::Tight) >= 2) {
            subChannel_ = Subchannel::EE;
        }  else {
            subChannel_ = Subchannel::MM;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}

bool Trigger_Eff::passMetTrigger()
{
    return (trig_cuts.pass_cut(Subchannel::MET)
            // || trig_cuts.pass_cut(Subchannel::MHT)
            );
}

bool Trigger_Eff::passDileptonTrigger()
{
    if (subChannel_ == Subchannel::EE) {
        return (trig_cuts.pass_cut_any(Subchannel::EE)
                || trig_cuts.pass_cut_any(Subchannel::Single_E));
    } else if (subChannel_ == Subchannel::MM) {
        return (trig_cuts.pass_cut_any(Subchannel::MM)
                || trig_cuts.pass_cut_any(Subchannel::Single_M));
    } else if (subChannel_ == Subchannel::EM) {
        return (trig_cuts.pass_cut_any(Subchannel::EM)
                || trig_cuts.pass_cut_any(Subchannel::Single_M)
                || trig_cuts.pass_cut_any(Subchannel::Single_E));
    } else {
        return false;
    }
}


bool Trigger_Eff::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (signal_cuts())
        (*currentChannel_) = Channel::Signal;

    if (trees.find(*currentChannel_) == trees.end()) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool Trigger_Eff::signal_cuts()
{
    LOG_FUNC << "Start of baseline_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());

    passCuts &= cuts.setCut("passLeadPt", getLeadPt() > 25);
    passCuts &= cuts.setCut("passLeptonNum", nLeps(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passLowMassCut", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));
    passCuts &= cuts.setCut("passLeadPt", !isSameSign(Level::Tight));

    passCuts &= cuts.setCut("passMetTrigger", passMetTrigger());
    passCuts &= cuts.setCut("passMetCut", met.pt() > 130);
    // passCuts &= cuts.setCut("passMhtCut", jet.getMHT() > 130);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Loose) > 100);

    fillCutFlow(Channel::Signal, cuts);

    LOG_FUNC << "End of baseline_cuts";
    return passCuts;
}

void Trigger_Eff::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillOutput(*o_tightMuons, Level::Tight, event_bitmap);
    elec.fillOutput(*o_tightElectrons, Level::Tight, event_bitmap);

    for (size_t systNum = 0; systNum < numSystematics(); ++systNum) {
        size_t syst = syst_to_index.at(systNum);
        if (syst == 0 && systNum != 0) {
            continue;
        }
        setupSyst(systNum);

        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_ht.push_back(jet.getHT(Level::Loose));
        o_mht.push_back(jet.getMHT());
        o_dilepton_trigger.push_back(passDileptonTrigger());
        // o_met_trigger.push_back(passMetTrigger());
    }
    LOG_FUNC << "End of FillValues";
}
