#ifndef VARIABLE_H_
#define VARIABLE_H_

#include <bitset>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <stdexcept>

#include "analysis_suite/skim/interface/logging.h"

#define BITMAP_SIZE 64
typedef std::bitset<BITMAP_SIZE> Bitmap;

template <class T>
class TRVariable {
public:
    TRVariable() {}

    TRVariable(TTreeReader& fReader, std::string branch, bool needed=true) {
        setup(fReader, branch, needed);
    }

    void setup(TTreeReader& fReader, std::string branch, bool needed=true) {
        branch_ = branch;
        if (fReader.GetTree()->GetBranchStatus(branch.c_str())) {
            val = new TTreeReaderValue<T>(fReader, branch.c_str());

        } else if (!needed) {
            std::cout << "Branch (" << branch << ") not found and not needed, will continue to run" << std::endl;
        } else {
            throw std::out_of_range("Branch (" + branch + ") not found and needed, will stop");
        }
    }
    T operator*() const {
        if (val) {
            return **val;
        } else if (std::is_same<T, bool>::value) {
            return false;
        } else {
            return 0.;
        }
    }
protected:
    TTreeReaderValue<T>* val = nullptr;
    std::string branch_;
};


template <class T>
class OutValue : public TRVariable<T> {
public:
    void setupBranch(TTree* tree) {
        tree->Branch(this->branch_.c_str(), &var_raw);
    }

    void fill() {
        var_raw = this->operator*();
    }

private:
    T var_raw;
};



template <class T>
class TRArray {

public:
    TRArray() {}

    TRArray(TTreeReader& fReader, std::string branch) {
        setup(fReader, branch);
    }

    void setup(TTreeReader& fReader, std::string branch, bool needed=true) {
        if (fReader.GetTree()->GetBranchStatus(branch.c_str())) {
            array = new TTreeReaderArray<T>(fReader, branch.c_str());
        } else if (!needed) {
            std::cout << "Branch (" << branch << ") not found and not needed, will continue to run" << std::endl;
        } else {
            throw std::out_of_range("Branch (" + branch + ") not found and needed, will stop");
        }
    }
    T at(size_t idx) const {
        if (array) {
            return array->At(idx);
        } else if (std::is_same<T, bool>::value) {
            return false;
        } else {
            return 0.;
        }
    }
    size_t size() const {
        if (array) {
            return array->GetSize();
        } else {
            return 0;
        }
    }
    auto begin() {
        return array->begin();
    }
    auto end() {
        return array->end();
    }
private:
    TTreeReaderArray<T>* array = nullptr;
};




#endif // VARIABLE_H_
