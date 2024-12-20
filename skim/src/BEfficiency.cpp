#include "analysis_suite/skim/interface/BEfficiency.h"

#include "analysis_suite/skim/interface/logging.h"

namespace Channel {
    enum  {
        Signal,
        None,
    };
}

void BEfficiency::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    DileptonBase::Init(tree);

    createTree("Signal", Channel::Signal);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    hlt_mu.setup(fReader, "HLT_IsoMu24");
    hlt_el.setup(fReader, "HLT_Ele32_WPTight_Gsf");
    hlt_mm.setup(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
    hlt_em.setup(fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");
    hlt_me.setup(fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
    hlt_ee.setup(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    hlt_ee_high.setup(fReader, "HLT_DoubleEle25_CaloIdL_MW");

    LOG_FUNC << "End of Init";
}

void BEfficiency::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    // tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("Muons", "ParticleOut", &o_muon);
    tree->Branch("Electrons", "ElectronOut_Endcap", &o_elec);
    // tree->Branch("BJets", "BEffOut", &o_beff);
    // tree->Branch("HLT_dilepton_HT", &o_hlt_dilepton_ht);
    // tree->Branch("HLT_dilepton", &o_hlt_dilepton);
    // tree->Branch("HLT_trilepton", &o_hlt_trilepton);
    // tree->Branch("passZVeto_loose", &o_pass_zveto_loose);
    // tree->Branch("passZVeto_fake", &o_pass_zveto_fake);
    tree->Branch("hlt_mu", &o_hlt_mu);
    tree->Branch("hlt_el", &o_hlt_el);
    tree->Branch("hlt_mm", &o_hlt_mm);
    tree->Branch("hlt_em", &o_hlt_em);
    tree->Branch("hlt_me", &o_hlt_me);
    tree->Branch("hlt_ee", &o_hlt_ee);
    tree->Branch("hlt_ee_high", &o_hlt_ee_high);
    
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void BEfficiency::clearOutputs()
{
    o_hlt_mu.clear();
    o_hlt_el.clear();
    o_hlt_mm.clear();
    o_hlt_em.clear();
    o_hlt_me.clear();
    o_hlt_ee.clear();
    o_hlt_ee_high.clear();
}

void BEfficiency::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";

    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getScaleFactor();
    (*weight) *= jet.getTotalBTagWeight("M");
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

bool BEfficiency::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    setSubChannel();

    if (signal_cuts()) {
        (*currentChannel_) = Channel::Signal;
        return true;
    } else {
        (*currentChannel_) = Channel::None;
        return false;
    }
}


bool BEfficiency::signal_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());

    if (nLeps(Level::Tight) > 0) {
        passCuts &= cuts.setCut("passLeadPt", getLeadPt() > 25);
    }
    passCuts &= cuts.setCut("passFakeLeptonNum", nLeps(Level::Tight) >= 2);

    // passCuts &= cuts.setCut("passTrigger", getTriggerValue());
    passCuts &= cuts.setCut("passSSLeptons", isSameSign(Level::Tight));
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);
    passCuts &= cuts.setCut("passLowMassCut", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));

    return passCuts;
}

void BEfficiency::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    elec.fillElectron_Endcap(*o_elec, Level::Tight, event_bitmap);
    // jet.fillJet(*o_jets, Level::Loose, event_bitmap);
    muon.fillOutput(*o_muon, Level::Tight, event_bitmap);
    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        o_hlt_mu.push_back(*hlt_mu);
        o_hlt_el.push_back(*hlt_el);
        o_hlt_mm.push_back(*hlt_mm);
        o_hlt_me.push_back(*hlt_me);
        o_hlt_em.push_back(*hlt_em);
        o_hlt_ee.push_back(*hlt_ee);
        o_hlt_ee_high.push_back(*hlt_ee_high);
        // o_hlt_dilepton_ht.push_back(trig_cuts.pass_cut_any(Subchannel::All_Dilep_HT));
        // o_hlt_dilepton.push_back(trig_cuts.pass_cut_any(Subchannel::All_Dilep));
        // o_hlt_trilepton.push_back(trig_cuts.pass_cut_any(Subchannel::All_Multi));
        // o_pass_zveto_loose.push_back(!muon.isInMassRange(Level::Loose) && !elec.isInMassRange(Level::Loose));
        // o_pass_zveto_fake.push_back(!muon.isInMassRange(Level::Fake) && !elec.isInMassRange(Level::Fake));
    }
    LOG_FUNC << "End of FillValues";
}
