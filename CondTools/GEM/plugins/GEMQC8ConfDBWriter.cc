#include "CondCore/PopCon/interface/PopConAnalyzer.h"
#include "CondTools/GEM/interface/GEMQC8ConfSourceHandler.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef popcon::PopConAnalyzer<popcon::GEMQC8ConfSourceHandler> GEMQC8ConfDBWriter;
//define this as a plug-in
DEFINE_FWK_MODULE(GEMQC8ConfDBWriter);
