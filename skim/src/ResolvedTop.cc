#include "analysis_suite/skim/interface/ResolvedTop.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void ResolvedTop::setup(TTreeReader& fReader)
{
    GenericParticle::setup("FatJet", fReader);
    mass_softdrop.setup(fReader, "FatJet_msoftdrop");
    tau2.setup(fReader, "FatJet_tau2");
    tau3.setup(fReader, "FatJet_tau3");

    setup_map(Level::Loose);

    wp = wp_by_name.at("VL");
}

void ResolvedTop::createLooseList()
{
    for (size_t i = 0; i < size(); i++) {
        if (mass_softdrop.at(i) > 105
            && mass_softdrop.at(i) < 210
            && tau3.at(i)/tau2.at(i) < wp)
            m_partList[Level::Loose]->push_back(i);
    }
}

float ResolvedTop::getScaleFactor()
{
    float weight = 1.;
    return weight;
}

void ResolvedTop::fillTop(TopOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            output.disc.push_back(tau3.at(idx)/tau2.at(idx));
        }
    }
}
