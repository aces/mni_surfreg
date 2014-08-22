/*
   SURFACE DEPTH POTENTIAL

   Compute the depth potential for a cortical surface. This is the dataterm
   to be used for surface registration.

   Can also compute:
     - surface area (geometric or Voronoi)
     - mean curvature
     - Gaussian curvature
     - normal vectors
     - blurred surface field

   depth_potential [-SOR val] [-alpha val] [-normals] [-gaussian_curvature] 
                   [-mean_curvature] [-area_voronoi] [-area_simple] [-depth_potential] 
                   [-smooth <fwhm> <in_file.txt>] file.obj out_file.txt 

   Values: SOR = successive over-relaxation parameter for Gauss-Seidel solver 
                 (1 < val < 2, default=1.90)
           alpha = damping factor for depth potential (val > 0, default=0.0015)
           normals = compute the surface normals
           gaussian_curvature = compute the Gaussian curvature
           mean_curvature = compute the mean curvature
           area_voronoi = compute the surface area using Voronoi's method (default)
           area_simple = compute the surface area using equal weights
           depth_potential = compute the depth potential
           fwhm = full width half maximum for surface smoothing
           in_file.txt = input field to smooth

   References:

   M. Boucher, S. Whitesides, and A. Evans, "Depth potential function for 
     folding pattern representation, registration and analysis", Medical
     Image Analysis, Volume 13, Issue 2, pp 203-214, April 2009.
     
   AUTHOR:    Maxime Boucher

   HISTORY:   VERSION 1.0  February 12, 2007 (conversion from Matlab to C,
                           by Claude Lepage)
              VERSION 1.1  July 24, 2007 (addition of surface smoothing,
                           by Claude Lepage)

   COPYRIGHT: Copyright Alan C. Evans
              Professor of Neurology
              McGill University

              Corresponding Address: boucher@bic.mni.mcgill.ca
                                     claude@bic.mni.mcgill.ca
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <volume_io.h>
#include <bicpl.h>

#define vec_sub( a, b, c ) { \
  c.coords[0] = a.coords[0] - b.coords[0]; \
  c.coords[1] = a.coords[1] - b.coords[1]; \
  c.coords[2] = a.coords[2] - b.coords[2]; \
}

#define vec_dot_product( a, b ) \
  ( a.coords[0] * b.coords[0] + a.coords[1] * b.coords[1] + a.coords[2] * b.coords[2] )

#define vec_normalize( a ) { \
  Real aa = sqrt( a.coords[0] * a.coords[0] + a.coords[1] * a.coords[1] + \
                  a.coords[2] * a.coords[2] ); \
  a.coords[0] /= aa; \
  a.coords[1] /= aa; \
  a.coords[2] /= aa; \
}

#define vec_cross( a, b, c ) { \
  c.coords[0] = a.coords[1] * b.coords[2] - a.coords[2] * b.coords[1]; \
  c.coords[1] = a.coords[2] * b.coords[0] - a.coords[0] * b.coords[2]; \
  c.coords[2] = a.coords[0] * b.coords[1] - a.coords[1] * b.coords[0]; \
}

struct csr_matrix {
  int   n;        // size of matrix
  int   nnz;      // number of non-zero coeffs
  int * ia;       // row pointers
  int * ja;       // column pointers
  Real * A;       // coefficients
};

enum { SN, GC, MC, AV, AS, DP, SM };

// Prototypes of local functions.

private void usage( char * );
private Status read_surface_obj( STRING, int *, Point *[],
                                 Vector *[], int *[], int **[] );
private Status get_surface_neighbours( polygons_struct *, int *[],
                                       int ** [] );
private void save_vector( int n_points, Vector vec[], char * out_file );
private Real * read_scalar( int n_points, char * in_file );
private void save_scalar( int n_points, Real scalar[], char * out_file );

private void init_csr_matrix( int n_points, int * n_ngh, int ** ngh, 
                              struct csr_matrix * mat );
private void free_csr_matrix( struct csr_matrix * mat );
private void assemble( Real val, int row, int col, struct csr_matrix * mat );
private void print_csr_matrix( struct csr_matrix * mat );
private void gauss_seidel( int nnode, int * ia, int * ja, Real * mat, Real * rhs, 
                           Real * soln, int max_iter, double tol, double SOR,
                           int verbose );

private void stable_normals( int n_points, Point coords[], Vector normals[], 
                             int * n_ngh, int ** ngh );
private Real * compute_areas( int n_points, Point coords[], int * n_ngh, 
                              int ** ngh, int lambda );
private Real * compute_gaussian_curvature( int n_points, Point coords[], int * n_ngh, 
                                           int ** ngh, Real * areas );
private Real * compute_mean_curvature( int n_points, Point coords[], Real * areas,
                                       Vector normals[] , struct csr_matrix * mat );
private void cot_laplacian_operator( int n_points, Point coords[], struct csr_matrix * mat,
                                     int * n_ngh, int ** ngh );
private Real * compute_depth_potential( int n_points, Point coords[], Real * areas,
                                        struct csr_matrix * mat, Real * mc, 
                                        double alpha, double SOR );
private void smooth( int n_points, Point coords[], Real * areas,
                     struct csr_matrix * mat, Real * scalar, 
                     double fwhm );

// -------------------------------------------------------------------
// Help message on how to use this module.
//
private void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s [-SOR <val>] [-alpha <val>] [-normals] [-gaussian_curvature] [-mean_curvature] \n\
       [-area_voronoi] [-area_simple] [-depth_potential] [-smooth <fwhm> <in_file.txt>] \n\
       <file.obj> <out_file.txt> \n\n\
Values: -SOR <val> : successive over-relaxation parameter for Gauss-Seidel solver (1 < val < 2, default=1.90)\n\
        -alpha <val> : damping factor for depth potential (val > 0, default=0.0015)\n\
        -normals : compute the surface normals \n\
        -gaussian_curvature : compute the Gaussian curvature \n\
        -mean_curvature : compute the mean curvature \n\
        -area_voronoi : compute the surface area using Voronoi's method (default) \n\
        -area_simple : compute the surface area using equal weights \n\
        -depth_potential : compute the depth potential \n\
        -smooth <fwhm> <file> : full width half maximum for surface smoothing and \n\
                                file name for input field to smooth \n\n";

  print_error( usage_format, executable_name );
}

// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
private Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_neighbours[],
                                 int ** neighbours[] ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  File_formats      format;
  STRING            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( ERROR );
  }

  if( n_objects != 1 ||
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  // Make a copy of the coordinates and the normals, since
  // delete_object_list will destroy them.

  *n_points = surface->n_points;
  ALLOC( *points, surface->n_points );
  ALLOC( *normals, surface->n_points );
  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  get_surface_neighbours( surface, n_neighbours, neighbours );

  delete_object_list( n_objects, object_list );

  return( OK );
}

// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
private Status get_surface_neighbours( polygons_struct * surface,
                                       int * n_neighbours_return[],
                                       int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  ALLOC( ngh, surface->n_points );
  ALLOC( ngh_array, 3*surface->n_items );

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

  return OK;

}

// -------------------------------------------------------------------
// Save a vector to a file.
//
private void save_vector( int n_points, Vector vec[], char * out_file ) {

  int i;

  FILE * fp = fopen( out_file, "w" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp, "%g %g %g\n", vec[i].coords[0], vec[i].coords[1], vec[i].coords[2] );
  }
  fclose( fp );
}

// -------------------------------------------------------------------
// Read a scalar from a file.
//
private Real * read_scalar( int n_points, char * in_file ) {

  int i;

  Real * scalar = (Real*)malloc( n_points * sizeof( Real ) );
  FILE * fp = fopen( in_file, "r" );
  for( i = 0; i < n_points; i++ ) {
    fscanf( fp, "%lg\n", &scalar[i] );
  }
  fclose( fp );
  return( scalar );
}

// -------------------------------------------------------------------
// Save a scalar to a file.
//
private void save_scalar( int n_points, Real scalar[], char * out_file ) {

  int i;

  FILE * fp = fopen( out_file, "w" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp, "%g\n", scalar[i] );
  }
  fclose( fp );
}

// -------------------------------------------------------------------
// Initialize the data structures for the sparse matrix in CSR format
// based on the connectivity (assumed to be triangles). The matrix is
// initialized to zero.
//
private void init_csr_matrix( int n_points, int * n_ngh, int ** ngh, 
                              struct csr_matrix * mat ) {

  int i, j, nnz;

  mat->nnz = n_points;
  for( i = 0; i < n_points; i++ ) {
    mat->nnz += n_ngh[i];
  }

  mat->n = n_points;
  mat->ia = (int *)malloc( ( n_points + 1 ) * sizeof( int ) );
  mat->ja = (int *)malloc( mat->nnz * sizeof( int ) );
  mat->A = (Real *)malloc( mat->nnz * sizeof( Real ) );

  for( i = 0; i < mat->nnz; i++ ) {
    mat->A[i] = 0.0;
  }

  nnz = 0;
  mat->ia[0] = 0;
  for( i = 0; i < n_points; i++ ) {
    mat->ia[i+1] = mat->ia[i] + n_ngh[i] + 1;
    mat->ja[nnz] = i;
    nnz++;
    for( j = 0; j < n_ngh[i]; j++ ) {
      mat->ja[nnz] = ngh[i][j];
      nnz++;
    }
  }
}

// -------------------------------------------------------------------
// Free the memory for the data structures of a sparse matrix in CSR 
// format.
//
private void free_csr_matrix( struct csr_matrix * mat ) {

  free( mat->ia );
  free( mat->ja );
  free( mat->A );

  mat->nnz = 0;
  mat->n = 0;
  mat->ia = NULL;
  mat->ja = NULL;
  mat->A = NULL;
}

// -------------------------------------------------------------------
// Add a coefficient to a symmetric sparse matrix in CSR format, at
// (row,col) and (col,row).
//
private void assemble( Real val, int row, int col, struct csr_matrix * mat ) {

  int j;

  for( j = mat->ia[row]; j < mat->ia[row+1]; j++ ) {
    if( col == mat->ja[j] ) mat->A[j] += val;
  }
  for( j = mat->ia[col]; j < mat->ia[col+1]; j++ ) {
    if( row == mat->ja[j] ) mat->A[j] += val;
  }
}

// -------------------------------------------------------------------
// Print a sparse matrix in CSR format, for debugging purposes only.
// Output is stdout.
//
private void print_csr_matrix( struct csr_matrix * mat ) {

  int  i, j;

  printf( "%d %d\n", mat->n, mat->n );
  for( i = 0; i < mat->n; i++ ) {
    for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
      printf( "%d %d %g\n", mat->ja[j]+1, i+1, mat->A[j] );
    }
  }
}

// -------------------------------------------------------------------
// Solve a linear system with SOR Gauss-Seidel. The matrix is in CSR
// format. The convergence threshold is a relative tolerance from the
// initial residual, up to a maximum number of iterations, whichever 
// comes first.
//
private void gauss_seidel( int nnode, int * ia, int * ja, Real * mat, 
                           Real * rhs, Real * soln, int max_iter, 
                           double tol, double SOR, int verbose ) {

  int i, j, iter;
  double res, res0, val, diag;

  for( i = 0; i < nnode; i++ ) {
    soln[i] = 0;
  }
 
  for( iter = 1; iter <= max_iter; iter++ ) { 
    res = 0.0;
    for( i = 0; i < nnode; i++ ) {
      val = rhs[i];
      for( j = ia[i]; j < ia[i+1]; j++ ) {
        if( i != ja[j] ) {
          val -= mat[j] * soln[ja[j]];
        } else {
          diag = mat[j];
        }
      }
      val /= diag;
      val = soln[i] + SOR * ( val - soln[i] );
      res += fabs( val - soln[i] );
      soln[i] = val;
    }
    res /= (float)nnode;
    if( iter == 1 ) res0 = res;
    res /= res0;
    if( verbose && iter%10 == 1 ) printf( "Iteration %d  Res %g\n", iter, res );
    if( res < tol ) break;
  }
  if( verbose ) printf( "Last iteration %d  Res %g\n", iter, res );
}

// -------------------------------------------------------------------
// Compute the normal vectors.
//
// function res = stable_normals(coords,tri)
// 
// nb_points = size(coords,2);
// nb_tri = size(tri,2);
// v1 = normalize(coords(:,tri(2,:)) - coords(:,tri(1,:)));
// v2 = normalize(coords(:,tri(3,:)) - coords(:,tri(2,:)));
// v3 = normalize(coords(:,tri(1,:)) - coords(:,tri(3,:)));
// 
// n = normalize(cross(v1,v2));
// 
// w(1,:) = (1-dot_product(v3,v1))./norm_vec(cross(v3,v1));
// w(2,:) = (1-dot_product(v1,v2))./norm_vec(cross(v1,v2));
// w(3,:) = (1-dot_product(v2,v3))./norm_vec(cross(v2,v3));
// 
// res(1,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
//                       [tri(1,:),tri(2,:),tri(3,:)],
//                       [n(1,:).*w(1,:),n(1,:).*w(2,:),n(1,:).*w(3,:)]));
// res(2,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
//                       [tri(1,:),tri(2,:),tri(3,:)],
//                       [n(2,:).*w(1,:),n(2,:).*w(2,:),n(2,:).*w(3,:)]));
// res(3,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
//                       [tri(1,:),tri(2,:),tri(3,:)],
//                       [n(3,:).*w(1,:),n(3,:).*w(2,:),n(3,:).*w(3,:)]));
// res = normalize(full(res));
//
private void stable_normals( int n_points, Point coords[], Vector normals[], 
                             int * n_ngh, int ** ngh ) {

  int     i, j, n1, n2, nn;
  Real    mag1, mag2, mag3, mag3x1, mag1x2, mag2x3, w1, w2, w3;
  Vector  v1, v2, v3, cross;

  for( i = 0; i < n_points; i++ ) {
    normals[i].coords[0] = 0.0;
    normals[i].coords[1] = 0.0;
    normals[i].coords[2] = 0.0;
  }

  for( i = 0; i < n_points; i++ ) {
    for( nn = 0; nn < n_ngh[i]; nn++ ) {

      n1 = ngh[i][nn];
      n2 = ngh[i][(nn+1)%n_ngh[i]];
      if( i < n1 && i < n2 ) {

        vec_sub( coords[n1], coords[i], v1 );
        vec_sub( coords[n2], coords[n1], v2 );
        vec_sub( coords[i], coords[n2], v3 );
        vec_normalize( v1 );
        vec_normalize( v2 );
        vec_normalize( v3 );

        // weights
        vec_cross( v3, v1, cross );
        mag3x1 = sqrt( vec_dot_product( cross, cross ) );
        vec_cross( v2, v3, cross );
        mag2x3 = sqrt( vec_dot_product( cross, cross ) );
        vec_cross( v1, v2, cross );
        mag1x2 = sqrt( vec_dot_product( cross, cross ) );

        w1 = ( 1.0 - vec_dot_product( v3, v1 ) ) / mag3x1;
        w2 = ( 1.0 - vec_dot_product( v1, v2 ) ) / mag1x2;
        w3 = ( 1.0 - vec_dot_product( v2, v3 ) ) / mag2x3;

        // contribution of tri_norm to each node of triangle
        for( j = 0; j < 3; j++ ) {
          normals[i].coords[j] += w1 * cross.coords[j] / mag1x2;
          normals[n1].coords[j] += w2 * cross.coords[j] / mag1x2;
          normals[n2].coords[j] += w3 * cross.coords[j] / mag1x2;
        }
      }
    }
  }

  // normalization of normals
  for( i = 0; i < n_points; i++ ) {
    vec_normalize( normals[i] );
  }
}

// -------------------------------------------------------------------
// Compute Voronoi area around nodes.
//
//     if lambda ==1
//         area_mat = sparse([1:nb_tri,1:nb_tri,1:nb_tri],[tri(1,:),tri(2,:),tri(3,:)],
//                           (1/6)*[area_tri,area_tri,area_tri]);
//     else
//         %mixed area
//         good = (dot_product(v1,v2)<=0).*(dot_product(v2,v3)<=0).*(dot_product(v3,v1)<=0);
//         lgood = find(good);
//         lbad = find(1-good);
//         area_mat = sparse([lgood,lgood,lgood],[tri(1,lgood),tri(2,lgood),tri(3,lgood)],
//                   -0.125*[(norms(3,lgood).*dot_product(v1(:,lgood),v2(:,lgood))+
//                            norms(1,lgood).*dot_product(v2(:,lgood),v3(:,lgood)))./
//                           area_tri(lgood),
//                           (norms(2,lgood).*dot_product(v1(:,lgood),v3(:,lgood))+
//                            norms(1,lgood).*dot_product(v2(:,lgood),v3(:,lgood)))./
//                           area_tri(lgood),
//                           (norms(2,lgood).*dot_product(v1(:,lgood),v3(:,lgood))+
//                            norms(3,lgood).*dot_product(v1(:,lgood),v2(:,lgood)))
//                            ./area_tri(lgood)],nb_tri,nb_points);
//         area_mixed = sparse([lbad,lbad,lbad],[tri(1,lbad),tri(2,lbad),tri(3,lbad)],
//                      0.125*[area_tri(lbad),area_tri(lbad),area_tri(lbad)],nb_tri,nb_points);
//         area_mat = area_mat + area_mixed;
//         bad1 = find(dot_product(v3,v1)>0);
//         bad2 = find(dot_product(v1,v2)>0);
//         bad3 = find(dot_product(v2,v3)>0);
//         area_mixed = sparse([bad1,bad2,bad3],[tri(1,bad1),tri(2,bad2),tri(3,bad3)],
//                      0.125*[area_tri(bad1),area_tri(bad2),area_tri(bad3)],nb_tri,nb_points);
//         area_mat = area_mat + area_mixed;
//     end
//     areas = full(sum(area_mat));
// 
private Real * compute_areas( int n_points, Point coords[], int * n_ngh, 
                              int ** ngh, int lambda ) {

  int     i, n1, n2, nn;
  Real    area_tri;
  Vector  v1, v2, v3, cross;

  Real * areas = (Real*)malloc( n_points * sizeof( Real ) );
  for( i = 0; i < n_points; i++ ) {
    areas[i] = 0.0;
  }

  // Loop over the triangles to compute the area around nodes.

  for( i = 0; i < n_points; i++ ) {
    for( nn = 0; nn < n_ngh[i]; nn++ ) {
      n1 = ngh[i][nn];
      n2 = ngh[i][(nn+1)%n_ngh[i]];
      if( i < n1 && i < n2 ) {
        vec_sub( coords[n1], coords[i], v1 );
        vec_sub( coords[n2], coords[n1], v2 );
        vec_sub( coords[i], coords[n2], v3 );
        vec_cross( v1, v2, cross );
        area_tri = sqrt( vec_dot_product( cross, cross ) );
        if( lambda == 1 ) {
          // equal weights
          areas[i] += area_tri / 6.0;
          areas[n1] += area_tri / 6.0;
          areas[n2] += area_tri / 6.0;
        } else {
          // Voronoi
          int bad1 = vec_dot_product( v3, v1 ) > 0;
          int bad2 = vec_dot_product( v1, v2 ) > 0;
          int bad3 = vec_dot_product( v2, v3 ) > 0;
          int bad = bad1 || bad2 || bad3;
          if( !bad ) {
            Real w1 = vec_dot_product( v3, v3 ) * vec_dot_product( v1, v2 ) +
                      vec_dot_product( v1, v1 ) * vec_dot_product( v2, v3 );
            Real w2 = vec_dot_product( v2, v2 ) * vec_dot_product( v1, v3 ) +
                      vec_dot_product( v1, v1 ) * vec_dot_product( v2, v3 );
            Real w3 = vec_dot_product( v2, v2 ) * vec_dot_product( v1, v3 ) +
                      vec_dot_product( v3, v3 ) * vec_dot_product( v1, v2 );
            areas[i] -= 0.125 * w1 / area_tri;
            areas[n1] -= 0.125 * w2 / area_tri;
            areas[n2] -= 0.125 * w3 / area_tri;
          } else {
            areas[i] += 0.125 * area_tri;
            areas[n1] += 0.125 * area_tri;
            areas[n2] += 0.125 * area_tri;
            if( bad1 ) areas[i] += 0.125 * area_tri;
            if( bad2 ) areas[n1] += 0.125 * area_tri;
            if( bad3 ) areas[n2] += 0.125 * area_tri;
          }
        }
      }
    }
  }
  return( areas );
}

// -------------------------------------------------------------------
// Compute Gaussian curvature.
//
//     v1 = normalize(v1);
//     v2 = normalize(v2);
//     v3 = normalize(v3);
// 
//     gc = (2*pi - sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
//                             [tri(1,:),tri(2,:),tri(3,:)],
//                             [acos(dot_product(-v1,v3)),
//                              acos(dot_product(-v1,v2)),
//                              acos(dot_product(-v2,v3))],
//                             nb_tri,nb_points)))./areas;
//
private Real * compute_gaussian_curvature( int n_points, Point coords[], int * n_ngh, 
                                           int ** ngh, Real * areas ) {

  int     i, n1, n2, nn;
  Vector  v1, v2, v3;
  Real    two_pi;

  two_pi = 4.0 * atan2( 1.0, 0.0 );

  Real * gc = (Real*)malloc( n_points * sizeof( Real ) );
  for( i = 0; i < n_points; i++ ) {
    gc[i] = two_pi;
  }

  // Loop over the triangles to compute the Gaussian curvature.

  for( i = 0; i < n_points; i++ ) {
    for( nn = 0; nn < n_ngh[i]; nn++ ) {
      n1 = ngh[i][nn];
      n2 = ngh[i][(nn+1)%n_ngh[i]];
      if( i < n1 && i < n2 ) {
        vec_sub( coords[n1], coords[i], v1 );
        vec_sub( coords[n2], coords[n1], v2 );
        vec_sub( coords[i], coords[n2], v3 );
        vec_normalize( v1 );
        vec_normalize( v2 );
        vec_normalize( v3 );
        gc[i] -= acos( -vec_dot_product( v1, v3 ) );
        gc[n1] -= acos( -vec_dot_product( v1, v2 ) );
        gc[n2] -= acos( -vec_dot_product( v2, v3 ) );
      }
    }
  }
  for( i = 0; i < n_points; i++ ) {
    gc[i] /= areas[i];
  }
  return( gc );
}

// -------------------------------------------------------------------
// Compute mean curvature.
//
// mc obtained using: mean_curvature white.obj mc.dat
// which is: [mat,area] = cot_laplacian_operator(coords,tri);
//           n = stable_normals(coords,tri);
//           output_vv(output_file,0.5*dot_product(coords*mat,n)./area);
//
private Real * compute_mean_curvature( int n_points, Point coords[], Real * areas,
                                       Vector normals[], struct csr_matrix * mat ) {

  int     i, j;
  Vector  vv;

  Real * mc = (Real*)malloc( n_points * sizeof( Real ) );
  for( i = 0; i < n_points; i++ ) {
    mc[i] = 0.0;
  }

  // Loop over the triangles to compute the mean curvature.

  for( i = 0; i < n_points; i++ ) {
    vv.coords[0] = 0.0;
    vv.coords[1] = 0.0;
    vv.coords[2] = 0.0;
    for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
      vv.coords[0] += mat->A[j] * coords[mat->ja[j]].coords[0];
      vv.coords[1] += mat->A[j] * coords[mat->ja[j]].coords[1];
      vv.coords[2] += mat->A[j] * coords[mat->ja[j]].coords[2];
    }
    mc[i] = 0.5 * vec_dot_product( vv, normals[i] ) / areas[i];
  }
  return( mc );
}

// -------------------------------------------------------------------
// Creation of the matrix for the Laplacian operator.
//
// function [mat,areas,mc,area_mat,gc] = cot_laplacian_operator(coords,tri,lambda)
// 
// nb_points = size(coords,2);
// nb_tri = size(tri,2);
// v1 = coords(:,tri(2,:)) - coords(:,tri(1,:));
// v2 = coords(:,tri(3,:)) - coords(:,tri(2,:));
// v3 = coords(:,tri(1,:)) - coords(:,tri(3,:));
// 
// norms(1,:) = norm_vec(v1).^2;
// norms(2,:) = norm_vec(v2).^2;
// norms(3,:) = norm_vec(v3).^2;
// 
// area_tri = norm_vec(cross(v1,v2));
// mat = sparse([1:nb_points,tri(2,:),tri(3,:),tri(1,:)],
//              [1:nb_points,tri(1,:),tri(2,:),tri(3,:)],
//              [ones(1,nb_points),dot_product(v2,v3)./area_tri,
//               dot_product(v1,v3)./area_tri,dot_product(v1,v2)./area_tri]);
// mat = 0.5*(mat + mat');
// mat = mat - sparse(1:nb_points,1:nb_points,sum(mat));
//
private void cot_laplacian_operator( int n_points, Point coords[], struct csr_matrix * mat,
                                     int * n_ngh, int ** ngh ) {

  int     i, n1, n2, nn;
  Real    area_tri;
  Vector  v1, v2, v3, cross;

  // Loop over the triangles to compute the coefficients of the matrix.

  for( i = 0; i < n_points; i++ ) {
    for( nn = 0; nn < n_ngh[i]; nn++ ) {
      n1 = ngh[i][nn];
      n2 = ngh[i][(nn+1)%n_ngh[i]];
      if( i < n1 && i < n2 ) {
        vec_sub( coords[n1], coords[i], v1 );
        vec_sub( coords[n2], coords[n1], v2 );
        vec_sub( coords[i], coords[n2], v3 );

        vec_cross( v1, v2, cross );
        area_tri = sqrt( vec_dot_product( cross, cross ) );

        Real dot23 = vec_dot_product( v2, v3 ) / area_tri;
        Real dot13 = vec_dot_product( v1, v3 ) / area_tri;
        Real dot12 = vec_dot_product( v1, v2 ) / area_tri;

        // Assemble dot23 at (i,n1) and (n1,i); dot13 at (n1,n2) and (n2,n1);
        // dot12 at (i,n2) and (n2,i). assemble() is symmetric.
        assemble( dot23, i, n1, mat );
        assemble( dot13, n1, n2, mat );
        assemble( dot12, i, n2, mat );

        assemble( -0.5*(dot23+dot12), i, i, mat );
        assemble( -0.5*(dot13+dot23), n1, n1, mat );
        assemble( -0.5*(dot12+dot13), n2, n2, mat );
      }
    }   
  }
}

// -------------------------------------------------------------------
// Compute depth potential.
//
// Linear system is:
// [mat,area] = cot_laplacian_operator(coords,triangles);
// mat = 0.5*mat;
// mc = mean_curvature(coords,triangles);
// c = (mc*area-area*sum(mc*area)/sum(area));
// mat = mat + alpha*sparse(1:nb_coords,1:nb_coords,area);
// 
// [R]=cholinc(mat,1e-3);
// res = pcg(mat,c',1e-12,10000,R',R)';
//
private Real * compute_depth_potential( int n_points, Point coords[], Real * areas,
                                        struct csr_matrix * mat, Real * mc, 
                                        double alpha, double SOR ) {

  int i, j;

  Real sum_area = 0.0;
  Real sum_area_mc = 0.0;
  for( i = 0; i < n_points; i++ ) {
    sum_area += areas[i];
    sum_area_mc += mc[i] * areas[i];
  }

  Real * rhs = (Real*)malloc( n_points * sizeof( Real ) );
  Real factor = sum_area_mc / sum_area;
  for( i = 0; i < n_points; i++ ) {
    rhs[i] = ( mc[i] - factor ) * areas[i];
  }

  Real * dp_mat_array = (Real*)malloc( mat->nnz * sizeof( Real ) );

  for( i = 0; i < n_points; i++ ) {
    for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
      dp_mat_array[j] = 0.50 * mat->A[j];
      if( mat->ja[j] == i ) {
        dp_mat_array[j] += alpha * areas[i];
      }
    }
  }

  // Solution by SOR Gauss Seidel

  Real * dp = (Real*)malloc( n_points * sizeof( Real ) );

  gauss_seidel( n_points, mat->ia, mat->ja, dp_mat_array, rhs, dp, 1000, 1.0e-10, SOR, 1 );

  free( dp_mat_array );
  free( rhs );

  return( dp );
}


// -------------------------------------------------------------------
// Smoothing of a field on a surface
//
// Linear system is:
// [mat,area] = cot_laplacian_operator(coords,triangles);
// mat = 0.5*mat;
// mat = mat + alpha*sparse(1:nb_coords,1:nb_coords,area);
//
private void smooth( int n_points, Point coords[], Real * areas,
                     struct csr_matrix * mat, Real * scalar, 
                     double fwhm ) {

  int i, j;

  int  Niter = 100;
  Real Tfinal = fwhm * fwhm / ( 16.0 * log( 2.0 ) );
  Real dT = Tfinal / (float)Niter;
  double SOR = 1.0;   // best value for implicit smoothing

  Real * sm_mat_array = (Real*)malloc( mat->nnz * sizeof( Real ) );
  Real * sm = (Real*)malloc( n_points * sizeof( Real ) );

  for( i = 0; i < n_points; i++ ) {
    for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
      sm_mat_array[j] = dT * 0.50 * mat->A[j] / areas[i];
      if( mat->ja[j] == i ) {
        sm_mat_array[j] += 1.0;
      }
    }
    sm[i] = scalar[i];
  }

  for( i = 0; i < Niter; i++ ) {
    for( j = 0; j < n_points; j++ ) {
      scalar[j] = sm[j];
    }
    gauss_seidel( n_points, mat->ia, mat->ja, sm_mat_array, scalar, sm, 1000, 1.0e-8, SOR, 0 );
  }

  for( j = 0; j < n_points; j++ ) {
    scalar[j] = sm[j];
  }

  free( sm );
  free( sm_mat_array );

}


// Main program.

int main(int argc, char* argv[]) {

  char   * cortex_file_name;   // file name for cortex file
  char   * out_file_name;      // file name for output
  char   * smooth_file_name;   // file name for field to smooth
  double   SOR;                // successive over-relaxation coefficient
  double   alpha;              // damping factor for depth potential
  double   fwhm;               // full width half maximum for blurring
  int      lambda;             // parameter for areas: 1 or 2
  int      option;             // solver option
  int      n_points;           // number of grid points
  Point  * coords = NULL;      // coordinates
  Vector * normals = NULL;     // normal vectors
  Real   * areas = NULL;       // Voronoi area at nodes
  Real   * gc = NULL;          // Gaussian curvature
  Real   * mc = NULL;          // mean curvature
  Real   * dp = NULL;          // depth potential
  Real   * sm = NULL;          // input field for smoothing
  int    * n_ngh = NULL;       // node neighbours (inverse connectivity)
  int   ** ngh = NULL;

  struct csr_matrix laplacian; // the laplacian matrix operator
  STRING   arg;

  // Defaults
  cortex_file_name = NULL;
  out_file_name = NULL;
  SOR = 1.90;       // this value works best with alpha=0.0015
  alpha = 0.0015;
  lambda = 2;

  // Parse the command line arguments for the options.

  initialize_argument_processing( argc, argv );

  while( get_string_argument( NULL, &arg ) ) {

    if( equal_strings( arg, "-SOR" ) ) {
      if( !get_real_argument( 1.90, &SOR ) ) {
        print_error( "Error in -SOR arguments.\n" );
        usage( argv[0] );
        return( ERROR );
      }
    } else if( equal_strings( arg, "-alpha" ) ) {
      if( !get_real_argument( 0.0015, &alpha ) ) {
        print_error( "Error in -alpha arguments.\n" );
        usage( argv[0] );
        return( ERROR );
      }
    } else if( equal_strings( arg, "-normals" ) ) {
      option = SN;
    } else if( equal_strings( arg, "-gaussian_curvature" ) ) {
      option = GC;
    } else if( equal_strings( arg, "-mean_curvature" ) ) {
      option = MC;
    } else if( equal_strings( arg, "-area_voronoi" ) ) {
      option = AV;
      lambda = 2;
    } else if( equal_strings( arg, "-area_simple" ) ) {
      option = AS;
      lambda = 1;
    } else if( equal_strings( arg, "-depth_potential" ) ) {
      option = DP;
    } else if( equal_strings( arg, "-smooth" ) ) {
      if( !get_real_argument( 10.0, &fwhm ) ) {
        print_error( "Error in -smooth arguments.\n" );
        usage( argv[0] );
        return( ERROR );
      }
      if( !get_string_argument( NULL, &smooth_file_name ) ) {
        print_error( "Error in -smooth arguments.\n" );
        usage( argv[0] );
        return( ERROR );
      }
      option = SM;
    } else {
      if( cortex_file_name == NULL ) {
        cortex_file_name = arg;
      } else if( out_file_name == NULL ) {
        out_file_name = arg;
      } else {
        usage( argv[0] );
        return( ERROR );
      }
    }
  }

  // Check validity of the arguments.

  if( cortex_file_name == NULL ) {
    print_error( "Must supply a name for the surface file.\n" );
    usage( argv[0] );
    return( ERROR );
  }
  if( out_file_name == NULL ) {
    print_error( "Must supply a name for the output file.\n" );
    usage( argv[0] );
    return( ERROR );
  }

  // Read the surface file.

  if( read_surface_obj( cortex_file_name, &n_points,
                        &coords, &normals, &n_ngh, &ngh ) != OK ) {
    return( ERROR );
  }

  // Compute normal vectors.
  if( option == SN ) {
    stable_normals( n_points, coords, normals, n_ngh, ngh );
    save_vector( n_points, normals, out_file_name );
  }

  // Compute Voronoi areas.
  if( option == AV || option == AS ) {
    areas = compute_areas( n_points, coords, n_ngh, ngh, lambda );
    save_scalar( n_points, areas, out_file_name );
  }

  // Gaussian curvature.
  if( option == GC ) {
    areas = compute_areas( n_points, coords, n_ngh, ngh, lambda );
    gc = compute_gaussian_curvature( n_points, coords, n_ngh, ngh, areas );
    save_scalar( n_points, gc, out_file_name );
  }

  // Mean curvature.
  if( option == MC ) {
    stable_normals( n_points, coords, normals, n_ngh, ngh );
    areas = compute_areas( n_points, coords, n_ngh, ngh, lambda );
    init_csr_matrix( n_points, n_ngh, ngh, &laplacian );
    cot_laplacian_operator( n_points, coords, &laplacian, n_ngh, ngh );
    mc = compute_mean_curvature( n_points, coords, areas, normals, &laplacian );
    save_scalar( n_points, mc, out_file_name );
  }

  // Depth potential.
  if( option == DP ) {
    stable_normals( n_points, coords, normals, n_ngh, ngh );
    areas = compute_areas( n_points, coords, n_ngh, ngh, lambda );
    init_csr_matrix( n_points, n_ngh, ngh, &laplacian );
    cot_laplacian_operator( n_points, coords, &laplacian, n_ngh, ngh );
    mc = compute_mean_curvature( n_points, coords, areas, normals, &laplacian );
    dp = compute_depth_potential( n_points, coords, areas, &laplacian, mc, 
                                  alpha, SOR );
    save_scalar( n_points, dp, out_file_name );
  }

  // Surface smoothing.
  if( option == SM ) {
    sm = read_scalar( n_points, smooth_file_name );
    areas = compute_areas( n_points, coords, n_ngh, ngh, lambda );
    init_csr_matrix( n_points, n_ngh, ngh, &laplacian );
    cot_laplacian_operator( n_points, coords, &laplacian, n_ngh, ngh );
    smooth( n_points, coords, areas, &laplacian, sm, fwhm );
    save_scalar( n_points, sm, out_file_name );
  }

  if( option == MC || option == DP ) {
    free_csr_matrix( &laplacian );
  }

  if( coords ) FREE( coords );
  if( normals ) FREE( normals );
  if( n_ngh ) FREE( n_ngh );
  if( ngh ) {
    FREE( ngh[0] );   // this is ngh_array
    FREE( ngh );
  }
  if( areas ) free( areas );
  if( gc ) free( gc );
  if( mc ) free( mc );
  if( dp ) free( dp );
  if( sm ) free( sm );

  return( OK );
}

