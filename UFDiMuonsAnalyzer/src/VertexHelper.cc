
#import "UfHMuMuCode/UFDiMuonsAnalyzer/interface/VertexHelper.h"

void FillVertexInfos( VertexInfos& _vertexInfos,
		      const reco::VertexCollection verticesSelected ) {

  _vertexInfos.init();
  int nVertices = verticesSelected.size();

  for (int i = 0; i < nVertices; i++) {

    reco::Vertex vertex = verticesSelected.at(i);
    VertexInfo _vertexInfo;
    _vertexInfo.init();
    
    _vertexInfo.isValid  = vertex.isValid();
    _vertexInfo.x        = vertex.position().X();
    _vertexInfo.y        = vertex.position().Y();
    _vertexInfo.z        = vertex.position().Z();
    _vertexInfo.xErr     = vertex.xError();
    _vertexInfo.yErr     = vertex.yError();
    _vertexInfo.zErr     = vertex.zError();
    _vertexInfo.chi2     = vertex.chi2();
    _vertexInfo.ndof     = vertex.ndof();
    _vertexInfo.normChi2 = vertex.normalizedChi2();

    _vertexInfos.vertices.push_back( _vertexInfo );
    _vertexInfos.nVertices += 1;
  } // End loop: for (int i = 0; i < nVertices; i++)                                                                                            

  if ( _vertexInfos.nVertices != int(_vertexInfos.vertices.size()) )
    std::cout << "Bizzare error: _vertexInfos.nVertices = " << _vertexInfos.nVertices
	      << ", _vertexInfos.vertices.size() = " << _vertexInfos.vertices.size() << std::endl;
  
} // End function: void FillVertexInfo     
  

reco::VertexCollection SelectVertices( const edm::Handle<reco::VertexCollection>& vertices, const double _vertex_ndof_min,
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
} // End function: reco::VertexCollection select
