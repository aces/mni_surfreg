/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

// -*- C++ -*-


namespace MNI {


template <class _Surface>
template <class ControlIndexMap, class TargetIndexMap>
void 
SurfaceMap<_Surface>::write( const std::string& filename,
			     const ControlIndexMap& control_index,
			     const TargetIndexMap& target_index ) const
{
    std::ofstream os(filename.c_str());
    write( os, control_index, target_index );
}


template <class _Surface>
void 
SurfaceMap<_Surface>::write( const std::string& filename ) const
{
    std::ofstream os(filename.c_str());
    MNI::write_surface_map( os, *this );
}



template <class _Surface>
template <class ControlIndexMap, class TargetIndexMap>
void
SurfaceMap<_Surface>::write( std::ostream& os,
                             const ControlIndexMap& control_index,
                             const TargetIndexMap& target_index ) const
{
    MNI::write_surface_map( os, *this, control_index, target_index );
}


}

