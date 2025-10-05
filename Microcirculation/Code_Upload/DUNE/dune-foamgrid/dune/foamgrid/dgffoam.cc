// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include "dgffoam.hh"

namespace Dune
{

  /*namespace dgf
  {

    // Implementation of FoamGridParameterBlock
    // --------------------------------------
    struct FoamGridParameterBlock
      : public GridParameterBlock
    {
      FoamGridParameterBlock ( std::istream &input )
      { }
    }

  } // namespace dgf

  */

  // Implementation of DGFGridFactory< FoamGrid >
  // -------------------------------------------

  template< int dim , int dimworld>
  void DGFGridFactory< FoamGrid< dim, dimworld  > >::generate ( std::istream &input )
  {
    dgf_.element = DuneGridFormatParser::General;

    if( !dgf_.readDuneGrid( input, dim, dimworld  ) )
      DUNE_THROW( DGFException, "Error: Failed to build grid");

    dgf_.setOrientation( 0, 1 );

    // get grid parameter block
    //dgf::FoamGridParameterBlock gridParam( input );

    // create grid here to set heap size
    // create grid factory (passed grid is returned by createGrid method)
    //if( gridParam.heapSize() > 0 )
    //  FoamGrid< dim, dimworld  >::setDefaultHeapSize( gridParam.heapSize() );

    for( int n = 0; n < dgf_.nofvtx; n++ )
    {
      FieldVector< double, dimworld > v;
      for( int j = 0; j < dimworld; j++ )
        v[ j ] = dgf_.vtx[ n ][ j ];
      factory_.insertVertex( v );
    }

    std::vector< unsigned int > el;
    for( int n = 0; n < dgf_.nofelements; n++ )
    {
      el.clear();
      for( size_t j = 0; j < dgf_.elements[ n ].size(); ++j )
        el.push_back( ( dgf_.elements[ n ][ j ] ) );

      // simplices
      if( el.size() == std::size_t( dim+1 ) )
        factory_.insertElement( GeometryType( GeometryType::simplex, dim ), el );
      else
        DUNE_THROW( DGFException, "Invalid number of element vertices: " << el.size() );
    }

    // create grid
    grid_ = factory_.createGrid();
  }

} // namespace Dune
