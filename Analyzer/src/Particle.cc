#include "analysis_suite/Analyzer/interface/Particle.h"

void GenericParticle::setup(std::string name, TTreeReader& fReader)
{
    m_pt.setup(fReader, (name + "_pt").c_str());
    m_eta.setup(fReader, (name + "_eta").c_str());
    m_phi.setup(fReader, (name + "_phi").c_str());
    m_mass.setup(fReader, (name + "_mass").c_str());
}

void GenericParticle::clear() {
    for (auto& [key, plist] : m_partList) {
        plist->clear();
    }
}

void GenericParticle::setup_map(Level level) {
    m_partList[level] = new std::vector<size_t>(); // maybe try for virtual soon?
}

size_t GenericParticle::idx(Level level, size_t i) const
{
    if (i >= size(level)) {
        throw std::out_of_range("Particle of level {} has size ("
                                + std::to_string(size(level))
                                + ") and not large enough to look at index ("
                                + std::to_string(i) + ")" );
    }
    return list(level).at(i);
}

void Particle::clear()
{
    for (auto& [key, plist] : m_partArray) {
        m_bitArray[key].clear();
        for (size_t i = 0; i < nSyst; ++i) {
            plist[i].clear();
        }
    }
}

void Particle::setup_map(Level level)
{
    m_partArray[level] = PartList(nSyst);
    m_partList[level] = nullptr;
    m_bitArray[level] = std::vector<size_t>();
}

void Particle::moveLevel(Level level_start, Level level_end)
{
    m_partList[level_end] = m_partList[level_start];
    m_bitArray[level_end] = m_bitArray[level_start];
}

void Particle::xorLevel(Level big, Level small, Level target)
{
    auto small_list = list(small);
    for (auto idx: list(big)) {
        if (std::find(small_list.begin(), small_list.end(), idx) != small_list.end()) {
            continue;
        }
        m_partList[target]->push_back(idx);
        m_bitArray[target][idx] = m_bitArray[big][idx];
    }
}
