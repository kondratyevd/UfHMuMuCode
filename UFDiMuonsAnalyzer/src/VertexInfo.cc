
#import "UfHMuMuCode/UFDiMuonsAnalyzer/interface/VertexInfo.h"

void VertexInfo::init() {

  isValid  = -999;
  x        = -999;
  y        = -999;
  z        = -999;
  xErr     = -999;
  yErr     = -999;
  zErr     = -999;
  chi2     = -999;
  ndof     = -999;
  normChi2 = -999;

}

void VertexInfos::init() {

  nVertices = 0;
  vertices.clear();

}
