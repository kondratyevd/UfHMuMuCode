
#import"UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _VertexInfo::init(){

  nVertices = 0;
  for (unsigned int i = 0; i < arraySize; i++) {
    isValid[i]  = -999;
    x[i]        = -999;
    y[i]        = -999;
    z[i]        = -999;
    xErr[i]     = -999;
    yErr[i]     = -999;
    zErr[i]     = -999;
    chi2[i]     = -999;
    ndof[i]     = -999;
    normChi2[i] = -999;
  }

}

void _VertexInfo::fill( const reco::VertexCollection verticesSelected ) {

  nVertices = verticesSelected.size();

  for (int i = 0; i < nVertices; i++) {

    if ( i >= int(arraySize) ) {
      // Suppress printout because we only usually save the primary vertex
      // std::cout << "Found " << i+1 << "th vertex; only " << arraySize << " allowed in array" << std::endl;
      return;
    }

    reco::Vertex vertex = verticesSelected.at(i);

    isValid[i]  = vertex.isValid();
    x[i]        = vertex.position().X();
    y[i]        = vertex.position().Y();
    z[i]        = vertex.position().Z();
    xErr[i]     = vertex.xError();
    yErr[i]     = vertex.yError();
    zErr[i]     = vertex.zError();
    chi2[i]     = vertex.chi2();
    ndof[i]     = vertex.ndof();
    normChi2[i] = vertex.normalizedChi2();
  } // End loop: for (int i = 0; i < nVertices; i++)

}

reco::VertexCollection _VertexInfo::select( const edm::Handle<reco::VertexCollection>& vertices, const double _vertex_ndof_min,
					   const double _vertex_rho_max, const double _vertex_z_max ) {

  // Following no official recipe: should find? - AWB 12.11.16
  // Modeled after https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc#L284

  reco::VertexCollection verticesSelected;
  verticesSelected.clear();

  if ( !vertices.isValid() ) {
    std::cout << "No valid vertex collection" << std::endl;
    return verticesSelected;
  }

  for (reco::VertexCollection::const_iterator vertex = vertices->begin(), verticesEnd = vertices->end(); vertex !=verticesEnd; ++vertex) {

    if ( vertex->isFake()                                    ) continue;
    if ( vertex->ndof()                   < _vertex_ndof_min ) continue;
    if ( vertex->position().Rho()         > _vertex_rho_max  ) continue;
    if ( fabs( vertex->position().Z() )   > _vertex_z_max    ) continue;

    verticesSelected.push_back(*vertex);
  }
  
  return verticesSelected;
}
