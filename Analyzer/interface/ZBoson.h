#ifndef ZBOSON_H_
#define ZBOSON_H_

#include "analysis_suite/Analyzer/interface/Particle.h"
#include "analysis_suite/Analyzer/interface/Lepton.h"

class ZBoson : public Particle {
public:

    size_t _size() const override { return zList.size(); }
    float _pt(size_t idx) const override { return zList.at(idx).Pt(); }
    float _eta(size_t idx) const override { return zList.at(idx).Eta(); }
    float _phi(size_t idx) const override { return zList.at(idx).Phi(); }
    float _mass(size_t idx) const override { return zList.at(idx).M(); }

    void setup();

    void createTightList(Lepton& leptons, Level level);

    void setupGoodLists(Lepton& muons, Lepton& elecs)
    {
        createTightList(muons, Level::Muon);
        createTightList(elecs, Level::Electron);
    }

    void clear() override
    {
        Particle::clear();
        zList.clear();
        for (auto& [level, vec] : lep_pairs) {
            vec.clear();
        }
    }

protected:
    std::vector<LorentzVector> zList;
    std::unordered_map<Level, std::vector<std::pair<size_t, size_t>>> lep_pairs;

    const float ZMASS = 91.188;
    const float ZWINDOW = 15;
    const float LOW_ENERGY_CUT = 12;

};

#endif // ZBOSON_H_
