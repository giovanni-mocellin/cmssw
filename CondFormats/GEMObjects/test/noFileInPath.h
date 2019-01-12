#ifndef noFileInPath_H
#define noFileInPath_H

#include <string>
#include <iostream>

// remove dependency on FWCore ... FileInPath
#define FWCore_ParameterSet_FileInPath_h
#define COND_SERIALIZABLE

namespace edm {
  struct FileInPath {

    FileInPath(std::string fname) : fullFName() {
      const char *t="CondFormats/GEMObjects/data/";
      fname.replace(fname.find(t),strlen(t),"");
      fullFName= "../data/" + fname;
    }

    const std::string& fullPath() const { return fullFName; }
  public:
    std::string fullFName;
  };
};

#endif
