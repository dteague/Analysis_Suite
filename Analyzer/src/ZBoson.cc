#include "analysis_suite/Analyzer/interface/ZBoson.h"

void ZBoson::setup()
{
    lep_pairs[Level::Muon] = std::vector<std::pair<size_t, size_t>>();
    lep_pairs[Level::Electron] = std::vector<std::pair<size_t, size_t>>();
}

void ZBoson::createTightList(Lepton& leptons, Level level)
{
    for (auto it1 = leptons.list(Level::Tight).begin(); it1 != leptons.list(Level::Tight).end(); ++it1) {
        for (auto it2 = it1+1; it2 != leptons.list(Level::Tight).end(); ++it2) {
            size_t lep1 = *it1;
            size_t lep2 = *it2;
            auto z = leptons.p4(lep1) + leptons.p4(lep2);
            if (leptons.charge(lep1) * leptons.charge(lep2) < 0
                && z.M() > 60 && z.M() < 120)  {
                zList.push_back(z);
                lep_pairs[level].push_back({lep1, lep2});
            }
        }
    }
}
