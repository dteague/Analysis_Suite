#ifndef COMMONSTRUCTS_H_
#define COMMONSTRUCTS_H_

#include "analysis_suite/skim/interface/Variable.h"
#include "analysis_suite/skim/interface/CommonEnums.h"

#include "TSelectorList.h"
#include "TH1.h"
#include "TH2.h"

#include <set>

enum class Channel;
enum class Subchannel;

struct TriggerInfo {
    std::unordered_map<Subchannel, std::vector<std::string>> trigger_names;
    std::unordered_map<Subchannel, std::vector<TRVariable<Bool_t>>> trigs;
    std::unordered_map<Subchannel, Dataset> dataset_by_chan;
    std::unordered_map<std::string, std::string> l1_by_trig;
    std::vector<Int_t> trig_output;
    Dataset dataset;

    void setup_channel(Subchannel chan, Dataset dataset_, TTreeReader& fReader, std::vector<std::string> trigger_names_ = {}) {
        trigger_names[chan] = trigger_names_;
        dataset_by_chan[chan] = dataset_;
        trigs[chan] = std::vector<TRVariable<Bool_t>>();
        for (auto trig_name: trigger_names_) {
            trigs[chan].push_back(TRVariable<Bool_t>(fReader, trig_name, false));
        }
    }

    void set_dataset(Dataset dataset_) {
        dataset = dataset_;
    }

    std::string trig_name(Subchannel chan, size_t i) {
        return trigger_names[chan].at(i);
    }

    bool pass_cut(Subchannel chan) {
        if (!pass_dataset(chan)) { // Dataset correct for the channel
            return false;
         }

        for (auto trig : trigs[chan]) {
            if (*trig) return true;
        }
        return false;
    }

    bool pass_cut_any(Subchannel chan) {
        for (auto trig : trigs[chan]) {
            if (*trig) return true;
        }
        return false;
    }

    bool dataset_or_trig(Subchannel chan) {
        if (pass_dataset(chan)) { // Dataset correct for the channel
            return true;
        }

        for (auto trig : trigs[chan]) {
            if (*trig) return true;
        }
        return false;
    }

    bool pass_dataset(Subchannel chan) {
        return dataset == Dataset::None || dataset_by_chan[chan] == dataset;
    }

    bool pass_cut(Subchannel chan, size_t i) {
        return pass_dataset(chan) && *trigs[chan].at(i);
    }
};

struct CutInfo {
    std::vector<std::pair<std::string, bool>> cuts;

    bool passCutFlow() {
        for (auto& [_, pass] : cuts) {
            if (!pass) return false;
        }
        return true;
    }

    bool setCut(std::string name, bool pass) {
        cuts.push_back(std::make_pair(name, pass));
        return pass;
    }

    size_t size() { return cuts.size(); }

    void clear() { cuts.clear(); }

    std::string name(size_t i) { return cuts.at(i).first; }

};



struct TreeInfo {
    TTree* tree;
    TH1F *cutflow, *cutflow_ind;
    bool initialize_axis = false;

    TreeInfo(std::string name,  TDirectory* outdir, TSelectorList* fOutput);

    void fillCutFlow(CutInfo& cuts, float weight);
};

struct TrigEff {
    TH2F* passEvents;
    TH2F* allEvents;

    void setup(TSelectorList* fOutput, std::string name, int nChans, int bins, float low, float high) {
        passEvents = new TH2F((name+"_pass").c_str(), (name+"_pass").c_str(), nChans, 0, nChans, bins, low, high);
        fOutput->Add(passEvents);
        allEvents = new TH2F((name+"_all").c_str(), (name+"_all").c_str(), nChans, 0, nChans, bins, low, high);
        fOutput->Add(allEvents);
    }

    void fill(float val, bool passTrig, Subchannel chan, float weight) {
        allEvents->Fill(static_cast<int>(chan), val, weight);
        if (passTrig) passEvents->Fill(static_cast<int>(chan), val, weight);
    }
};

template <typename T>
struct VarCollection {
    std::vector<TRVariable<T>> variables;

    void setup(TTreeReader& fReader, std::vector<std::string> names) {
        for (auto name : names) {
            variables.push_back(TRVariable<T>(fReader, name));
        }
    }

    bool pass() {
        for (auto var: variables) {
            if (!*var) return false;
        }
        return true;
    }
};



#endif // COMMONSTRUCTS_H_
