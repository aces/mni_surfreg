/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef DISJOINT_VERTEX_SETS_HPP
#define DISJOINT_VERTEX_SETS_HPP

#include <CGAL/Unique_hash_map.h>


namespace MNI {


/*! \brief Collection of disjoint vertex sets.
 *
 * Implemented using a parent-pointer forest.
 */
template <class Surface_>
class DisjointVertexSets
{
public:
    typedef Surface_ Surface;
    typedef typename Surface::Vertex_handle Vertex_handle;


    DisjointVertexSets(Surface& s) 
    {
	typename Surface::Vertex_iterator v = s.vertices_begin();
	CGAL_For_all(v,s.vertices_end()) {
	    parent[v] = v;
	    rank[v] = 0;
	    size[v] = 1;
	}
    }

    void make_union( Vertex_handle u, Vertex_handle v )
    {
	u = find_set(u);
	v = find_set(v);
	if ( u != v )
	    link( u,v );
    }

    Vertex_handle find_set( Vertex_handle v )
    {
	if ( parent[v] != v )
	    parent[v] = find_set(parent[v]);
	return parent[v];
    }

    bool same_component( Vertex_handle u, Vertex_handle v )
    {
	return find_set(u) == find_set(v);
    }

    template <class Out>
    void components( Surface& s, Out it )
    {
	typename Surface::Vertex_iterator v = s.vertices_begin();
	CGAL_For_all(v,s.vertices_end()) {
	    if (parent[v] == v)
		*it++ = v;
	}
    }

    int component_size( Vertex_handle v )
    {
	return size[find_set(v)];
    }

    
private:
    // Link together two set representatives.
    void link( Vertex_handle u, Vertex_handle v )
    {
	CGAL_assertion( parent[u] == u );
	CGAL_assertion( parent[v] == v );
	CGAL_assertion( u != v );

	if (rank[u] > rank[v]) {
	    parent[v] = u;
	    size[u] += size[v];
	} else {
	    parent[u] = v;
	    size[v] += size[u];
	    if (rank[u] == rank[v])
		++rank[v];
	}
    }

private:
    CGAL::Unique_hash_map<Vertex_handle,Vertex_handle> parent;
    CGAL::Unique_hash_map<Vertex_handle,int> rank;
    CGAL::Unique_hash_map<Vertex_handle,int> size;
};


}
#endif
