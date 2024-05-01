// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "GeomUtilities.hpp"
#include "ContactPlane.hpp"
#include "tribol/utils/Math.hpp"

#include "axom/core.hpp" 
#include "axom/slic.hpp" 

#include <float.h> 
#include <cmath>
#include <iostream> 

namespace tribol
{

void ProjectPointToPlane( const real x, const real y, const real z, 
                          const real nx, const real ny, const real nz,
                          const real ox, const real oy, const real oz, 
                          real& px, real& py, real& pz )
{
   // compute the vector from input point to be projected to 
   // the origin point on the plane
   real vx = x - ox;
   real vy = y - oy;
   real vz = z - oz;

   // compute the projection onto the plane normal
   real dist = vx * nx + vy * ny + vz * nz;

   // compute the projected coordinates of the input point
   px = x - dist*nx;
   py = y - dist*ny;
   pz = z - dist*nz;

   return;

} // end ProjectPointToPlane()

//------------------------------------------------------------------------------
void ProjectPointToSegment( const real x, const real y,
                            const real nx, const real ny,
                            const real ox, const real oy,
                            real& px, real& py )
{
   // compute the vector from input point to be projected to 
   // the origin point on the plane
   real vx = x - ox;
   real vy = y - oy;

   // compute the projection onto the plane normal
   real dist = vx * nx + vy * ny;

   // compute the projected coordinates of the input point
   px = x - dist * nx;
   py = y - dist * ny;

   return;

} // end ProjectPointToSegment()

//------------------------------------------------------------------------------
void PolyInterYCentroid( const int namax,
                         const real* const RESTRICT xa,
                         const real* const RESTRICT ya,
                         const int nbmax,
                         const real* const RESTRICT xb,
                         const real* const RESTRICT yb,
                         const int isym,
                         real & area,
                         real & ycent )
{

   real vol;

   // calculate origin shift to avoid roundoff errors
   real xorg = FLT_MAX;
   real yorg = FLT_MAX;
   real xa_min = FLT_MAX;
   real xa_max = -FLT_MAX;
   real ya_min = FLT_MAX;
   real ya_max = -FLT_MAX;
   real xb_min = FLT_MAX;
   real xb_max = -FLT_MAX;
   real yb_min = FLT_MAX;
   real yb_max = -FLT_MAX;

   real qy = 0.0;

   if (nbmax < 1 || namax < 1) {
      area = 0.0;
      vol = 0.0;
      ycent = 0.0;
      return;
   }

   for (int na = 0 ; na < namax ; ++na) {
      if (xa[na]<xa_min) {
         xa_min=xa[na];
      }
      if (ya[na]<ya_min) {
         ya_min=ya[na];
      }
      if (xa[na]>xa_max) {
         xa_max=xa[na];
      }
      if (ya[na]>ya_max) {
         ya_max=ya[na];
      }
      xorg = axom::utilities::min(xorg, xa[na]);
      yorg = axom::utilities::min(yorg, ya[na]);
   }
   for (int nb = 0 ; nb < nbmax ; ++nb) {
      if (xb[nb]<xb_min) {
         xb_min=xb[nb];
      }
      if (yb[nb]<yb_min) {
         yb_min=yb[nb];
      }
      if (xb[nb]>xb_max) {
         xb_max=xb[nb];
      }
      if (yb[nb]>yb_max) {
         yb_max=yb[nb];
      }
      xorg = axom::utilities::min(xorg, xb[nb]);
      yorg = axom::utilities::min(yorg, yb[nb]);
   }
   if (isym==1) {
      yorg = axom::utilities::max(yorg, 0.0);
   }

   area = 0.0;
   vol = 0.0;
   ycent = 0.0;
   if (xa_min>xb_max) {
      return;
   }
   if (xb_min>xa_max) {
      return;
   }
   if (ya_min>yb_max) {
      return;
   }
   if (yb_min>ya_max) {
      return;
   }

   // loop over faces of polygon a
   for (int na = 0 ; na < namax ; ++na) {
      int nap = (na+1)%namax;
      real xa1 = xa[na]  - xorg;
      real ya1 = ya[na]  - yorg;
      real xa2 = xa[nap] - xorg;
      real ya2 = ya[nap] - yorg;
      if (isym==1) {
         if (ya[na]<0.0 && ya[nap]<0.0) {
            continue;
         }
         if (ya[na]<0.0) {
            if (ya1!=ya2) {
               xa1 = xa1 - (ya1+yorg)*(xa2-xa1)/(ya2-ya1);
            }
            ya1 = -yorg;
         }
         else if (ya[nap]<0.0) {
            if (ya1!=ya2) {
               xa2 = xa2 - (ya2+yorg)*(xa1-xa2)/(ya1-ya2);
            }
            ya2 = -yorg;
         }
      }
      real dxa = xa2 - xa1;
      if (dxa==0.0) {
         continue;
      }
      real dya = ya2 - ya1;
      real slopea = dya/dxa;

      // loop over faces of polygon b
      for (int nb = 0 ; nb < nbmax ; ++nb) {
         int nbp = (nb+1)%nbmax;
         real xb1 = xb[nb]  - xorg;
         real yb1 = yb[nb]  - yorg;
         real xb2 = xb[nbp] - xorg;
         real yb2 = yb[nbp] - yorg;
         if (isym==1) {
            if (yb[nb]<0.0 && yb[nbp]<0.0) {
               continue;
            }
            if (yb[nb]<0.0) {
               if (yb1!=yb2) {
                  xb1 = xb1 - (yb1+yorg)*(xb2-xb1)/(yb2-yb1);
               }
               yb1 = -yorg;
            }
            else if (yb[nbp]<0.0) {
               if (yb1!=yb2) {
                  xb2 = xb2 - (yb2+yorg)*(xb1-xb2)/(yb1-yb2);
               }
               yb2 = -yorg;
            }
         }
         real dxb = xb2 - xb1;
         if (dxb==0.0) {
            continue;
         }
         real dyb = yb2 - yb1;
         real slopeb = dyb/dxb;

         // determine sign of volume of intersection
         real s = dxa * dxb;

         // calculate left and right coordinates of overlap
         real xl = axom::utilities::max(axom::utilities::min(xa1, xa2), axom::utilities::min(xb1, xb2) );
         real xr = axom::utilities::min(axom::utilities::max(xa1, xa2), axom::utilities::max(xb1, xb2) );
         if (xl>=xr) {
            continue;
         }
         real yla = ya1 + (xl-xa1)*slopea;
         real ylb = yb1 + (xl-xb1)*slopeb;
         real yra = ya1 + (xr-xa1)*slopea;
         real yrb = yb1 + (xr-xb1)*slopeb;
         real yl = axom::utilities::min(yla, ylb);
         real yr = axom::utilities::min(yra, yrb);

         real area1;
         real qy1;
         real ym;

         // check if lines intersect
         real dslope = slopea - slopeb;
         if (dslope!=0.0) {
            real xm = (yb1 - ya1 + slopea*xa1 - slopeb*xb1)/dslope;
            ym = ya1 + slopea*(xm-xa1);
            if (xm>xl && xm<xr) {
               // lines intersect, case ii
               area1 = 0.5 * copysign((yl+ym)*(xm-xl), s);
               real area2 = 0.5 * copysign((ym+yr)*(xr-xm), s);
               area  = area + area1 + area2;

               if (yl+ym>0) {
                  qy1 = 1.0/3.0*(ym+yl*yl/(yl+ym))*area1;
                  qy  = qy + qy1;
               }
               if (ym+yr>0) {
                  real qy2 = 1.0/3.0*(yr+ym*ym/(ym+yr))*area2;
                  qy  = qy + qy2;
               }

               if (isym==1) {
                  yl = yl + yorg;
                  ym = ym + yorg;
                  yr = yr + yorg;
                  vol = vol + copysign( (xm-xl)*(yl*yl+yl*ym+ym*ym)
                                        + (xr-xm)*(ym*ym+ym*yr+yr*yr), s) / 3.0;
               }
               continue;
            }
         }

         // lines do not intersect, case i
         area1 = 0.5 * copysign( (xr-xl)*(yr+yl), s);
         area  = area + area1;
         if (yl+yr>0) {
            qy1 = 1./3.0*(yr+yl*yl/(yl+yr))*area1;
            qy  = qy + qy1;
         }

         if (isym==1) {
            yl = yl + yorg;
            ym = ym + yorg;
            yr = yr + yorg;
            vol = vol + copysign( (xr-xl)*(yl*yl+yl*yr+yr*yr), s) / 3.0;
         }
      }
   }

   if (area != 0.0) {
      ycent = qy/area + yorg;
   }

   if (isym==0) {
      vol = area;
   }
   
   return;

} // end PolyInterYCentroid()

//------------------------------------------------------------------------------
void Local2DToGlobalCoords( real xloc, real yloc, 
                            real e1X, real e1Y, real e1Z,
                            real e2X, real e2Y, real e2Z,
                            real cX, real cY, real cZ,
                            real& xg, real& yg, real& zg )
{

   // This projection takes the two input local vector components and uses 
   // them as coefficients in a linear combination of local basis vectors. 
   // This gives a 3-vector with origin at the common plane centroid.
   real vx = xloc * e1X + yloc * e2X;
   real vy = xloc * e1Y + yloc * e2Y;
   real vz = xloc * e1Z + yloc * e2Z;
   
   // the vector in the global coordinate system requires the addition of the 
   // plane point vector (global Cartesian coordinate basis) to the previously 
   // computed vector
   xg = vx + cX;
   yg = vy + cY;
   zg = vz + cZ;

   return;

} // end Local2DToGlobalCoords()

//------------------------------------------------------------------------------
void GlobalTo2DLocalCoords( const real* const RESTRICT pX, 
                            const real* const RESTRICT pY, 
                            const real* const RESTRICT pZ,
                            real e1X, real e1Y, real e1Z,
                            real e2X, real e2Y, real e2Z,
                            real cX, real cY, real cZ,
                            real* const RESTRICT pLX, 
                            real* const RESTRICT pLY, int size )
{

   SLIC_ERROR_IF(size > 0 && (pLX == nullptr || pLY == nullptr),
                 "GlobalTo2DLocalCoords: local coordinate pointers are null");

   // loop over projected nodes
   for (int i=0; i<size; ++i) {

      // compute the vector between the point on the plane and the input plane point
      real vX = pX[i] - cX;
      real vY = pY[i] - cY;
      real vZ = pZ[i] - cZ;

      // project this vector onto the {e1,e2} local basis. This vector is 
      // in the plane so the out-of-plane component should be zero.
      pLX[i] = vX * e1X + vY * e1Y + vZ * e1Z; // projection onto e1
      pLY[i] = vX * e2X + vY * e2Y + vZ * e2Z; // projection onto e2
    
   }

   return;

} // end GlobalTo2DLocalCoords()

//------------------------------------------------------------------------------
void GlobalTo2DLocalCoords( real pX, real pY, real pZ,
                            real e1X, real e1Y, real e1Z,
                            real e2X, real e2Y, real e2Z,
                            real cX, real cY, real cZ,
                            real& pLX, real& pLY )
{
   // compute the vector between the point on the plane and the input plane point
   real vX = pX - cX;
   real vY = pY - cY;
   real vZ = pZ - cZ;

   // project this vector onto the {e1,e2} local basis. This vector is 
   // in the plane so the out-of-plane component should be zero.
   pLX = vX * e1X + vY * e1Y + vZ * e1Z; // projection onto e1
   pLY = vX * e2X + vY * e2Y + vZ * e2Z; // projection onto e2

   return;

} // end GlobalTo2DLocalCoords()

//------------------------------------------------------------------------------
void VertexAvgCentroid( const real* const RESTRICT x, 
                        const real* const RESTRICT y, 
                        const real* const RESTRICT z, 
                        const int numVert,
                        real& cX, real& cY, real& cZ )
{
   SLIC_ERROR_IF (numVert==0, "VertexAvgCentroid: numVert = 0.");

   // (re)initialize the input/output centroid components
   cX = 0.0;
   cY = 0.0;
   cZ = 0.0;

   // loop over nodes adding the position components
   real fac = 1.0 / numVert;
   for (int i=0; i<numVert; ++i) {
      cX += x[i];
      cY += y[i];
      if (z != nullptr)
      {
         cZ += z[i]; 
      }
   }

   // divide by the number of nodes to compute average
   cX *= fac;
   cY *= fac;
   cZ *= fac;

   return;

} // end VertexAvgCentroid()

//------------------------------------------------------------------------------
void VertexAvgCentroid( const real* const RESTRICT x, 
                        const int dim,
                        const int numVert,
                        real& cX, real& cY, real& cZ )
{
   SLIC_ERROR_IF(numVert==0, "VertexAvgCentroid: numVert = 0.");

   // (re)initialize the input/output centroid components
   cX = 0.0;
   cY = 0.0;
   cZ = 0.0;

   // loop over nodes adding the position components
   real fac = 1.0 / numVert;
   for (int i=0; i<numVert; ++i) {
      cX += x[dim*i];
      cY += x[dim*i+1];
      if (dim > 2)
      {
         cZ += x[dim*i+2];
      }
   }

   // divide by the number of nodes to compute average
   cX *= fac;
   cY *= fac;
   cZ *= fac;

   return;

} // end VertexAvgCentroid()

//------------------------------------------------------------------------------
void PolyAreaCentroid( const real* const RESTRICT x, 
                       const int dim,
                       const int numVert,
                       real& cX, real& cY, real& cZ )
{
   SLIC_ERROR_IF(numVert==0, "PolyAreaCentroid: numVert = 0.");

   // (re)initialize the input/output centroid components
   cX = 0.0;
   cY = 0.0;
   cZ = 0.0;

   // compute the vertex average centroid of the polygon in 
   // order to break it up into triangles
   real cX_poly, cY_poly, cZ_poly;
   VertexAvgCentroid( x, dim, numVert, cX_poly, cY_poly, cZ_poly);

   // loop over triangles formed from adjacent polygon vertices 
   // and the vertex averaged centroid
   real xTri[3] = { 0., 0., 0. };
   real yTri[3] = { 0., 0., 0. };
   real zTri[3] = { 0., 0., 0. };
   
   // assign all of the last triangle coordinates to the 
   // polygon's vertex average centroid
   xTri[2] = cX_poly;
   yTri[2] = cY_poly;
   zTri[2] = cZ_poly;
   real area_sum = 0.;
   for (int i=0; i<numVert; ++i) // loop over triangles
   {
      // group triangle coordinates
      int triId = i;
      int triIdPlusOne = (i == (numVert - 1)) ? 0 : triId + 1;
      xTri[0] = x[ dim * triId ];
      yTri[0] = x[ dim * triId + 1 ];
      zTri[0] = x[ dim * triId + 2 ];
      xTri[1] = x[ dim * triIdPlusOne ];
      yTri[1] = x[ dim * triIdPlusOne + 1 ];
      zTri[1] = x[ dim * triIdPlusOne + 2 ];

      // compute the area of the triangle
      real area_tri = Area3DTri( xTri, yTri, zTri );
      area_sum += area_tri;

      // compute the vertex average centroid of the triangle 
      real cX_tri, cY_tri, cZ_tri;
      VertexAvgCentroid( &xTri[0], &yTri[0], &zTri[0],
                         3, cX_tri, cY_tri, cZ_tri );

      cX += cX_tri * area_tri;
      cY += cY_tri * area_tri;
      cZ += cZ_tri * area_tri;
     
   }

   cX /= area_sum;
   cY /= area_sum;
   cZ /= area_sum;

   return;

} // end PolyAreaCentroid()

//------------------------------------------------------------------------------
void PolyCentroid( const real* const RESTRICT x, 
                   const real* const RESTRICT y,
                   const int numVert,
                   real& cX, real& cY )
{
   SLIC_ERROR_IF(numVert==0, "PolyAreaCentroid: numVert = 0.");

   // (re)initialize the input/output centroid components
   cX = 0.0;
   cY = 0.0;

   real area = 0.;

   for (int i=0; i<numVert; ++i)
   {
      int i_plus_one = (i == (numVert - 1)) ? 0 : i + 1;
      cX += (x[i] + x[i_plus_one])*(x[i]*y[i_plus_one] - x[i_plus_one]*y[i]);   
      cY += (y[i] + y[i_plus_one])*(x[i]*y[i_plus_one] - x[i_plus_one]*y[i]);
      area += (x[i]*y[i_plus_one] - x[i_plus_one]*y[i]);
   }

   area *= 1./2.;

   real fac = 1./(6.*area);
   cX *= fac;
   cY *= fac; 
 
} // end PolyCentroid()

//------------------------------------------------------------------------------
FaceGeomError Intersection2DPolygon( const real* const RESTRICT xA, 
                                     const real* const RESTRICT yA, 
                                     const int numVertexA, 
                                     const real* const RESTRICT xB, 
                                     const real* const RESTRICT yB, 
                                     const int numVertexB,
                                     real posTol, real lenTol, 
                                     real* RESTRICT * RESTRICT polyX, 
                                     real* RESTRICT * RESTRICT polyY, 
                                     int& numPolyVert, real& area, bool orientCheck )
{
   // for tribol, if you have called this routine it is because a positive area of 
   // overlap between two polygons (faces) exists. This routine does not perform a 
   // "proximity" check to determine if the faces are "close enough" to proceed with 
   // the full calculation. This can and probably should be added.

   // check to make sure the intersection polygon vertex pointers are null
   SLIC_ERROR_IF(*polyX != nullptr || *polyY != nullptr, 
                 "Intersection2DPolygon: expecting nullptr input arguments polyX, polyY.");

   // check numVertexA and numVertexB to make sure they are 3 (triangle) or more
   if (numVertexA < 3 || numVertexB < 3) 
   {
      SLIC_DEBUG( "Intersection2DPolygon(): one or more degenerate faces with < 3 vertices." );
      area = 0.0;
      return INVALID_FACE_INPUT; 
   }

   // check right hand rule ordering of polygon vertices. 
   // Note 1: This check is consistent with the ordering that comes from PolyReorder() 
   // of two faces with unordered vertices. 
   // Note 2: Intersection2DPolygon doesn't require consistent face vertex orientation
   // between faces, as long as each are 'ordered' (CW or CCW).
   if (orientCheck)
   {
      bool orientA = CheckPolyOrientation( xA, yA, numVertexA );
      bool orientB = CheckPolyOrientation( xB, yB, numVertexB );

      if (!orientA || !orientB)
      {
         SLIC_DEBUG( "Intersection2DPolygon(): check face orientations for face A." );
         return FACE_ORIENTATION;
      }
   }

   // determine minimum number of vertices (for use later)
   int numVertexMax = numVertexA >= numVertexB ? numVertexA : numVertexB;

   // allocate an array to hold ids of interior vertices 
   int interiorVAId[ numVertexA ];
   int interiorVBId[ numVertexB ];

   // initialize all entries in interior vertex array to -1
   initIntArray( &interiorVAId[0], numVertexA, -1 );
   initIntArray( &interiorVBId[0], numVertexB, -1 );

   // precompute the vertex averaged centroids for both polygons. 
   real xCA = 0.0;
   real yCA = 0.0;
   real xCB = 0.0;
   real yCB = 0.0;
   real zC = 0.0; // not required, only used as dummy argument in centroid routine

   VertexAvgCentroid( xA, yA, nullptr, numVertexA, xCA, yCA, zC );
   VertexAvgCentroid( xB, yB, nullptr, numVertexB, xCB, yCB, zC );

   // check to see if any of polygon A's vertices are in polygon B, and vice-versa. Track
   // which vertices are interior to the other polygon. Keep in mind that vertex 
   // coordinates are local 2D coordinates.
   int numVAI = 0;
   int numVBI = 0;

   // check A in B
   for (int i=0; i<numVertexA; ++i)
   {
      if (Point2DInFace( xA[i], yA[i], xB, yB, xCB, yCB, numVertexB ))
      {
         // interior A in B
         interiorVAId[i] = i;
         ++numVAI; 
      }
   }

   // check to see if ALL of A is in B; then A is the overlapping polygon.
   if (numVAI == numVertexA)
   {
      *polyX = new real[ numVertexA ];
      *polyY = new real[ numVertexA ];
      numPolyVert = numVertexA;
      for (int i=0; i<numVertexA; ++i)
      {
         (*polyX)[i] = xA[i];
         (*polyY)[i] = yA[i];
      }
      real* xVert = *polyX;
      real* yVert = *polyY;
      area = Area2DPolygon( xVert, yVert, numVertexA );
      return NO_FACE_GEOM_ERROR;
   }

   // check B in A
   for (int i=0; i<numVertexB; ++i) 
   {
      if (Point2DInFace( xB[i], yB[i], xA, yA, xCA, yCA, numVertexA) )
      {
         // interior B in A
         interiorVBId[i] = i;
         ++numVBI;
      }
   }

   // check to see if ALL of B is in A; then B is the overlapping polygon.
   if (numVBI == numVertexB)
   {
      *polyX = new real[ numVertexB ];
      *polyY = new real[ numVertexB ];
      numPolyVert = numVertexB;
      for (int i=0; i<numVertexB; ++i)
      {
         (*polyX)[i] = xB[i];
         (*polyY)[i] = yB[i];
      }
      real* xVert = *polyX;
      real* yVert = *polyY;
      area = Area2DPolygon( xVert, yVert, numVertexB );
      return NO_FACE_GEOM_ERROR;
   }

   // check for coincident interior vertices. That is, a vertex on A interior to 
   // B occupies the same point in space as a vertex on B interior to A. This is 
   // O(n^2), but the number of interior vertices is anticipated to be small
   // if we are at this location in the routine
   for (int i=0; i<numVertexA; ++i)
   {
      if (interiorVAId[i] != -1)
      {
         for (int j=0; j<numVertexB; ++j)
         {
            if (interiorVBId[j] != -1)
            {
              // compute the distance between interior vertices
              real distX = xA[i] - xB[j];
              real distY = yA[i] - yB[j];
              real distMag = magnitude( distX, distY );
              if (distMag < 1.E-15)
              {
                 // remove the interior designation for the vertex in polygon B
//                 SLIC_DEBUG( "Removing duplicate interior vertex id: " << j << ".\n" );
                 interiorVBId[j] = -1;
                 numVBI -= 1;
              }
            }
         }
      }
   }

   // determine the maximum number of intersection points

   // maximum number of vertices between the two polygons. 
   int maxSegInter = numVertexMax * numVertexMax;

   // allocate space to store the segment-segment intersection vertex coords. 
   // and a boolean array to indicate intersecting pairs
   real interX[ maxSegInter ];
   real interY[ maxSegInter ];
   bool intersect[ maxSegInter ];
   bool dupl; // boolean to indicate a segment-segment intersection that 
              // duplicates an existing interior vertex.
   bool interior [4];

   // initialize the interX and interY entries
   initRealArray( &interX[0], maxSegInter, 0. );
   initRealArray( &interY[0], maxSegInter, 0. );
   initBoolArray( &intersect[0], maxSegInter, false );
   dupl = false;

   // loop over segment-segment intersections to find the rest of the 
   // intersecting vertices. This is O(n^2), but segments defined by two 
   // nodes interior to the other polygon will be skipped. This will catch 
   // outlier cases.
   int interId = 0;

   // loop over A segments
   for (int ia=0; ia<numVertexA; ++ia)
   {
      int vAID1 = ia;
      int vAID2 = (ia == (numVertexA-1)) ? 0 : (ia + 1);
      
      // set boolean indicating which nodes on segment A are interior
      interior[0] = (interiorVAId[vAID1] != -1) ? true : false; 
      interior[1] = (interiorVAId[vAID2] != -1) ? true : false;
//      bool checkA = (interior[0] == -1 && interior[1] == -1) ? true : false;
      bool checkA = true;

      // loop over B segments
      for (int jb=0; jb<numVertexB; ++jb)
      {
         int vBID1 = jb;
         int vBID2 = (jb == (numVertexB-1)) ? 0 : (jb + 1);
         interior[2] = (interiorVBId[vBID1] != -1) ? true : false; 
         interior[3] = (interiorVBId[vBID2] != -1) ? true : false;
//         bool checkB = (interior[2] == -1 && interior[3] == -1) ? true : false;
         bool checkB = true;

         // if both segments are not defined by nodes interior to the other polygon
         if (checkA && checkB) 
         {
            if (interId >= maxSegInter) 
            {
               SLIC_DEBUG("Intersection2DPolygon: number of segment/segment intersections exceeds precomputed maximum; " << 
                          "check for degenerate overlap.");
               return DEGENERATE_OVERLAP;
            }

            intersect[interId] = SegmentIntersection2D( xA[vAID1], yA[vAID1], xA[vAID2], yA[vAID2],
                                                        xB[vBID1], yB[vBID1], xB[vBID2], yB[vBID2],
                                                        interior, interX[interId], interY[interId], 
                                                        dupl, posTol );
            ++interId;
         }
      } // end loop over A segments
   }  // end loop over B segments

   // count the number of segment-segment intersections
   int numSegInter = 0;
   for (int i=0; i<maxSegInter; ++i)
   {
      if (intersect[i]) ++numSegInter; 
   }

   // add check for case where there are no interior vertices or 
   // intersection vertices
   if (numSegInter == 0 && numVBI == 0 && numVAI == 0)
   {
      area = 0.0;
      return NO_FACE_GEOM_ERROR;
   }

   // allocate temp intersection polygon vertex coordinate arrays to consist 
   // of segment-segment intersections and number of interior points in A and B
   numPolyVert = numSegInter + numVAI + numVBI;
   real polyXTemp[ numPolyVert ];
   real polyYTemp[ numPolyVert ];

   // fill polyXTemp and polyYTemp with the intersection points
   int k = 0;
   for (int i=0; i<maxSegInter; ++i) 
   {
      if (intersect[i]) 
      {
         polyXTemp[k] = interX[i];
         polyYTemp[k] = interY[i];
         ++k;
      }
   }

   // fill polyX and polyY with the vertices on A that lie in B
   for (int i=0; i<numVertexA; ++i)
   {
      if (interiorVAId[i] != -1)
      {
         // debug
         if (k > numPolyVert)
         {
            SLIC_DEBUG("Intersection2DPolygon(): number of A vertices interior to B " << 
                       "polygon exceeds total number of overlap vertices. Check interior vertex id values.");
            return FACE_VERTEX_INDEX_EXCEEDS_OVERLAP_VERTICES;
         }

         polyXTemp[k] = xA[i];
         polyYTemp[k] = yA[i];
         ++k;
      }
   }

   for (int i=0; i<numVertexB; ++i)
   {
      if (interiorVBId[i] != -1)
      {
         // debug
         if (k > numPolyVert)
         {
            SLIC_DEBUG("Intersection2DPolygon(): number of B vertices interior to A " << 
                       "polygon exceeds total number of overlap vertices. Check interior vertex id values.");
            return FACE_VERTEX_INDEX_EXCEEDS_OVERLAP_VERTICES;
         }

         polyXTemp[k] = xB[i];
         polyYTemp[k] = yB[i];
         ++k;
      }
   }

   // order the unordered vertices (in counter clockwise fashion)
   PolyReorder( polyXTemp, polyYTemp, numPolyVert );

   // check length of segs against tolerance and collapse short segments if necessary
   // This is where polyX and polyY get allocated.
   int numFinalVert = 0; 

   FaceGeomError segErr = CheckPolySegs( polyXTemp, polyYTemp, numPolyVert, 
                                         lenTol, polyX, polyY, numFinalVert );

   if (segErr != 0)
   {
      return segErr;
   }

   numPolyVert = numFinalVert;

   // check to see if numPolyVert has degenerated to less than a triangle
   if (numPolyVert < 3)
   {
      area = 0.0;
      return NO_FACE_GEOM_ERROR; // don't return error here. We should tolerate 'collapsed' (zero area) overlaps
   }

   // compute the area of the polygon
   real* xVert = *polyX;
   real* yVert = *polyY;
   area = 0.0;
   area = Area2DPolygon( xVert, yVert, numPolyVert );

   return NO_FACE_GEOM_ERROR;

} // end Intersection2DPolygon()

//------------------------------------------------------------------------------
bool CheckPolyOrientation( const real* const RESTRICT x, 
                           const real* const RESTRICT y, 
                           const int numVertex )
{
   bool check = true;
   for (int i=0; i<numVertex; ++i)
   {
      // determine vertex indices of the segment
      int ia = i;
      int ib = (i == (numVertex-1)) ? 0 : (i+1);

      // compute segment vector
      real lambdaX = x[ib] - x[ia];
      real lambdaY = y[ib] - y[ia];
    
      // determine segment normal 
      real nrmlx = -lambdaY;
      real nrmly = lambdaX;

      // compute vertex-averaged centroid
      real* z = nullptr;
      real xc, yc, zc;
      VertexAvgCentroid( x, y, z, numVertex, xc, yc, zc );

      // compute vector between centroid and first vertex of current segment
      real vx = xc - x[ia];
      real vy = yc - y[ia];

      // compute dot product between segment normal and centroid-to-vertex vector.
      // the normal points inward toward the centroid
      real prod = vx * nrmlx + vy * nrmly;

      if (prod < 0.) // don't keep checking
      {
         check = false;
         return check;
      }
   }
   return check; // should equal true if here.

} // end CheckPolyOrientation()

//------------------------------------------------------------------------------
bool Point2DInFace( const real xPoint, const real yPoint, 
                    const real* const RESTRICT xPoly, 
                    const real* const RESTRICT yPoly,
                    const real xC, const real yC, 
                    const int numPolyVert )
{
   SLIC_ERROR_IF(numPolyVert<3, "Point2DInFace: number of face vertices is less than 3");

   SLIC_ERROR_IF(xPoly == nullptr || yPoly == nullptr, "Point2DInFace: input pointer not set");

   // if face is triangle (numPolyVert), call Point2DInTri once
   if (numPolyVert == 3)
   {
      return Point2DInTri( xPoint, yPoint, xPoly, yPoly );
   }

   // loop over triangles and determine if point is inside
   bool tri = false;
   for (int i=0; i<numPolyVert; ++i)
   {

      real xTri[3];
      real yTri[3];

      // construct polygon using i^th segment vertices and face centroid
      xTri[0] = xPoly[i];
      yTri[0] = yPoly[i];

      xTri[1] = (i == (numPolyVert-1)) ? xPoly[0] : xPoly[i+1];
      yTri[1] = (i == (numPolyVert-1)) ? yPoly[0] : yPoly[i+1];

      // last vertex of the triangle is the vertex averaged centroid of the polygonal face
      xTri[2] = xC;
      yTri[2] = yC;

      // call Point2DInTri for each triangle
      tri = Point2DInTri( xPoint, yPoint, xTri, yTri );

      if (tri) 
      {
         return true;
      }
   }
   return false;

} // end Point2DInFace()

//------------------------------------------------------------------------------
bool Point2DInTri( const real xp, const real yp, 
                   const real* const RESTRICT xTri, 
                   const real* const RESTRICT yTri )
{
   bool inside = false;

   // compute coordinate basis between the 1-2 and 1-3 vertices
   real e1x = xTri[1] - xTri[0];
   real e1y = yTri[1] - yTri[0];

   real e2x = xTri[2] - xTri[0];
   real e2y = yTri[2] - yTri[0];

   // compute vector components of vector between point and first vertex
   real p1x = xp - xTri[0];
   real p1y = yp - yTri[0];

   // compute dot products (e1,e1), (e1,e2), (e2,e2), (p1,e1), and (p1,e2)
   real e11 =  e1x * e1x + e1y * e1y;
   real e12 =  e1x * e2x + e1y * e2y;
   real e22 =  e2x * e2x + e2y * e2y;
   real p1e1 = p1x * e1x + p1y * e1y;
   real p1e2 = p1x * e2x + p1y * e2y;

   // compute the inverse determinant
   real invDet = 1.0 / (e11 * e22 - e12 * e12);

   // compute 2 local barycentric coordinates
   real u = invDet * (e22 * p1e1 - e12 * p1e2);
   real v = invDet * (e11 * p1e2 - e12 * p1e1); 
   
   // u or v may be negative, but numerically zero. Address this
   u = (std::abs(u) < 1.e-12) ? 0.0 : u;
   v = (std::abs(v) < 1.e-12) ? 0.0 : v;

   if ((u >= 0) && (v >= 0) && (u + v <= 1))
   {
      inside = true;
   }

   return inside;

} // end Point2DInTri()

//------------------------------------------------------------------------------
real Area2DPolygon( const real* const RESTRICT x, 
                    const real* const RESTRICT y, 
                    const int numPolyVert )
{

   real area = 0.;

   // compute vertex-averaged centroid to construct a triangle between segment 
   // vertices and centroid
   real* z = nullptr;
   real xc, yc, zc;
   VertexAvgCentroid(x, y, z, numPolyVert, xc, yc, zc);

   for (int i=0; i<numPolyVert; ++i)
   {
      // determine vertex indices of the segment
      int ia = i;
      int ib = (i == (numPolyVert-1)) ? 0 : (i+1);

      area += std::abs( 0.5 * (x[ia]*(y[ib]-yc) + 
                               x[ib]*(yc-y[ia]) + 
                               xc*(y[ia]-y[ib])));
   }
   return area;

} // end Area2DPolygon()

//------------------------------------------------------------------------------
real Area3DTri( const real* const RESTRICT x,
                const real* const RESTRICT y,
                const real* const RESTRICT z )
{
   real u[3] = { x[1] - x[0], y[1] - y[0], z[1] - z[0] };
   real v[3] = { x[2] - x[0], y[2] - y[0], z[2] - z[0] };

   return std::abs( 1./2. * magCrossProd( u, v ));

} // end Area3DTri()

//------------------------------------------------------------------------------
bool SegmentIntersection2D( const real xA1, const real yA1, const real xB1, const real yB1,
                            const real xA2, const real yA2, const real xB2, const real yB2,
                            const bool* const RESTRICT interior, real& x, real& y, 
                            bool& duplicate, const real tol )
{
   // note 1: this routine computes a unique segment-segment intersection, where two 
   // segments are assumed to intersect at a single point. A segment-segment overlap 
   // is a different computation and is not accounted for here. In the context of the 
   // use of this routine in the tribol polygon-polygon intersection calculation, 
   // two overlapping segments will have already registered the vertices that form 
   // the bounds of the overlapping length as vertices interior to the other polygon
   // and therefore will be in the list of overlapping polygon vertices prior to this 
   // routine. 
   //
   // note 2: any segment-segment intersection that occurs at a vertex of either segment 
   // will pass back the intersection coordinates, but will note a duplicate vertex. 
   // This is because that any vertex of polygon A that lies on a segment of polygon B 
   // will be caught and registered as a vertex interior to the other polygon and will 
   // be in the list of overlapping polygon vertices prior to calling this routine.

   // compute segment vectors
   real lambdaX1 = xB1 - xA1;
   real lambdaY1 = yB1 - yA1;

   real lambdaX2 = xB2 - xA2;
   real lambdaY2 = yB2 - yA2;

   real seg1Mag = magnitude( lambdaX1, lambdaY1 );
   real seg2Mag = magnitude( lambdaX2, lambdaY2 );

   // compute determinant of the lambda matrix, [ -lx1 -ly1, lx2 ly2 ]
   real det = -lambdaX1 * lambdaY2 + lambdaX2 * lambdaY1;

   // return false if det = 0. Check for numerically zero determinant
   real detTol = 1.E-12;
   if (det > -detTol && det < detTol)
   {
      x = 0.;
      y = 0.;
      duplicate = false;
      return false;
   }

   // compute intersection
   real invDet = 1.0 / det;
   real rX = xA1 - xA2;
   real rY = yA1 - yA2;
   real tA = invDet * (rX * lambdaY2 - rY * lambdaX2);
   real tB = invDet * (rX * lambdaY1 - rY * lambdaX1);

   // if tA and tB don't lie between [0,1] then return false.
   if ((tA < 0. || tA > 1.) || (tB < 0. || tB > 1.))
   {
      // no intersection
      x = 0.;
      y = 0.;
      duplicate = false;
      return false;
   }

   // TODO refine how these debug calculations are guarded
   {
     // debug check to make sure the intersection coordinates derived from 
     // each segment equation (scaled with tA and tB) are the same to some 
     // tolerance
     real xTest1 = xA1 + lambdaX1 * tA;
     real yTest1 = yA1 + lambdaY1 * tA;
     real xTest2 = xA2 + lambdaX2 * tB; 
     real yTest2 = yA2 + lambdaY2 * tB; 
     
     real xDiff = xTest1 - xTest2;
     real yDiff = yTest1 - yTest2;

     // make sure the differences are positive
     xDiff = (xDiff < 0.) ? -1.0 * xDiff : xDiff;
     yDiff = (yDiff < 0.) ? -1.0 * yDiff : yDiff;

     real diffTol = 1.0E-3;
     SLIC_DEBUG_IF( xDiff > diffTol || yDiff > diffTol, 
                   "SegmentIntersection2D(): Intersection coordinates are not equally derived." );
   }

   // if we get here then it means we have an intersection point.
   // Find the minimum distance of the intersection point to any of the segment 
   // vertices. 
   x = xA1 + lambdaX1 * tA;
   y = yA1 + lambdaY1 * tA;

   // for convenience, define an array of pointers that point to the 
   // input coordinates
   real xVert[4];
   real yVert[4];
   
   xVert[0] = xA1;
   xVert[1] = xB1;
   xVert[2] = xA2;
   xVert[3] = xB2;

   yVert[0] = yA1;
   yVert[1] = yB1;
   yVert[2] = yA2;
   yVert[3] = yB2;

   real distX[4];
   real distY[4];
   real distMag[4];

   for (int i=0; i<4; ++i)
   {
      distX[i] = x - xVert[i];
      distY[i] = y - yVert[i];
      distMag[i] = magnitude( distX[i], distY[i] );
   }

   real distMin = (seg1Mag > seg2Mag) ? seg1Mag: seg2Mag;
   int idMin;
   real xMinVert;
   real yMinVert;

   for (int i=0; i<4; ++i)
   {
      if (distMag[i] < distMin)
      {
         distMin = distMag[i];
         idMin = i;
         xMinVert = xVert[i];
         yMinVert = yVert[i];
      }
   }

   // check to see if the minimum distance is less than the position tolerance for 
   // the segments
   real distRatio = (idMin == 0 || idMin == 1) ? (distMin / seg1Mag) : (distMin / seg2Mag);

   // if the distRatio is less than the tolerance, or percentage cutoff of the original 
   // segment that we would like to keep, then check to see if the segment vertex closest 
   // to the computed intersection point is an interior point. If this is true, then collapse
   // the computed intersection point to the interior point and mark the duplicate boolean.
   // Also do this for the argument, interior, set to nullptr
   if (distRatio < tol)
   {
      if (interior == nullptr || interior[idMin])
      {
         x = xMinVert;
         y = yMinVert;
         duplicate = true;
         return false;
      }
   }

   // if we are here we are ready to return the true intersection point
   duplicate = false;
   return true;

} // end SegmentIntersection2D()

//------------------------------------------------------------------------------
FaceGeomError CheckPolySegs( const real* const RESTRICT x, const real* const RESTRICT y, 
                             const int numPoints, const real tol, 
                             real* RESTRICT * RESTRICT xnew, real* RESTRICT * RESTRICT ynew, 
                             int& numNewPoints )
{
   real newIDs[ numPoints ];

   // set newIDs[i] to original local ordering
   for (int i=0; i<numPoints; ++i)
   {
      newIDs[i] = i;
   }

   for (int i=0; i<numPoints; ++i)
   {
      // determine vertex indices of the segment
      int ia = i;
      int ib = (i == (numPoints-1)) ? 0 : (i+1);

      // compute segment vector magnitude
      real lambdaX = x[ib] - x[ia];
      real lambdaY = y[ib] - y[ia];
      real lambdaMag = magnitude( lambdaX, lambdaY );
     
      // check segment length against tolerance
      if (lambdaMag < tol)
      {
         // collapse second vertex to the first vertex of the current segment
         newIDs[ib] = i;
      }
   }

   // determine the number of new points
   numNewPoints = 0;
   for (int i=0; i<numPoints; ++i)
   {
      if (newIDs[i] == i)
      {
         ++numNewPoints;
      }
   }

   // allocate space for xnew and ynew. These are the input/output pointers 
   // to the overlapping polygon's vertex coordinates in the main polygon intersection 
   // routine
   *xnew = new real[ numNewPoints ];
   *ynew = new real[ numNewPoints ];
   
   // set the coordinates in xnew and ynew 
   int k = 0;
   for (int i=0; i<numPoints; ++i)
   {
      if (newIDs[i] == i)
      {
         if (k > numNewPoints)
         {
            SLIC_DEBUG("checkPolySegs(): index into polyX/polyY exceeds allocated space");
            return FACE_VERTEX_INDEX_EXCEEDS_OVERLAP_VERTICES;
         }

         (*xnew)[k] = x[i];
         (*ynew)[k] = y[i];
         ++k;
      }
   }

   return NO_FACE_GEOM_ERROR;

} // end CheckPolySegs()

//------------------------------------------------------------------------------
void PolyReorder( real* const RESTRICT x, real* const RESTRICT y, const int numPoints )
{

   SLIC_ERROR_IF(numPoints<3, "PolyReorder: numPoints < 3.");

   real xC, yC, zC;
   real * z = nullptr;
   real proj [numPoints - 2];

   int newIDs[ numPoints ];

   // initialize newIDs array to local ordering, 0,1,2,...,numPoints-1
   for (int i=0; i<numPoints; ++i)
   {
      newIDs[i] = i;
   }

   // compute vertex averaged centroid, in local coordinates
   VertexAvgCentroid( x, y, z, numPoints, xC, yC, zC );

   // using the first index into the x,y vertex coordinate arrays as 
   // the first vertex of the soon-to-be ordered list of vertices, determine 
   // the next vertex that will comprise the first segment in a counter 
   // clockwise ordering of vertices
   int id1 = -1;
   int id0 = 0;
   newIDs[0] = id0;
  
   for (int j=1; j<numPoints; ++j)
   {
      // determine segment vector and normal
      real lambdaX = x[j] - x[id0];
      real lambdaY = y[j] - y[id0];
      real nrmlx = -lambdaY;
      real nrmly = lambdaX;

      // project vectors that span from each point, except j,k, to first vertex (id0), onto the 
      // segment normal. There will always be numPoints-2 projections
      int pk = 0;
      for (int k=0; k<numPoints; ++k)
      {
         if (k != id0 && k != j)
         {
            proj[pk] = (x[k]-x[id0]) * nrmlx + (y[k]-y[id0]) * nrmly;
            ++pk;
         }
      }

      // check if all points are on one side of line defined by segment
      // (pk at this point should be equal to numPoints - 2)
      bool neg = false;
      bool pos = false;
      for (int ip=0; ip<pk; ++ip)
      {
         if (neg)
         {
            neg = true;
         }
         else if (!neg)
         {
            neg = (proj[ip] < 0.) ? true : false;
         }

         if (pos)
         {
            pos = true;
         }
         else if (!pos)
         {
            pos = (proj[ip] > 0.) ? true : false;
         }

         if (neg && pos)
         {
            break;
         }
      }

      // if one of the booleans is false then all points are on one side 
      // of line defined by i-j segment.
      if (!neg || !pos)
      {
         // check the orientation of the nodes to make sure we have the correct 
         // one of two segments that will pass the previous test.
         // Check the dot product between the normal and the vector
         // between the centroid and first (0th) vertex
         real vx = xC - x[id0];
         real vy = yC - y[id0];

         real prod = nrmlx * vx + nrmly * vy;

         // check if the two vertices are a segment on the convex hull and oriented CCW.
         // CCW orientation has prod > 0
         if (prod > 0) 
         {
            id1 = j;  
            break;
         }
      }

   } // end loop over j

   // swap ids
   if (id1 != -1)
   {
      newIDs[1] = id1;
      newIDs[id1] = 1;
   }

   // given the first (current) reference segment, compute the link vector between the jth vertex 
   // (j cannot be a vertex belonging to the reference segment) and the first vertex of 
   // the given reference segment. The next reference segment is between the second vertex of 
   // the current reference segment and the jth vertex whose link vector has the smallest 
   // dot product with the current reference segment.

   for (int i=0; i<(numPoints-3); ++i) // increment to (numPoints - 3) or (numPoints - 2)?
   {
      int jID;
      real cosThetaMax = -1.; // this handles angles up to 180 degrees. Not possible for convex polygons
      real cosTheta;
      real refMag, linkMag;

      // compute reference vector;
      real refx, refy;
      refx = x[newIDs[i+1]] - x[newIDs[i]];
      refy = y[newIDs[i+1]] - y[newIDs[i]];
      refMag = magnitude( refx, refy );

//      SLIC_ERROR_IF(refMag < 1.E-12, "PolyReorder: reference segment for link vector check is nearly zero length");

      // loop over link vectors of unassigned vertices
      int nextVertexID = 2+i;
      for (int j=nextVertexID; j<numPoints; ++j)
      {
         real lx, ly;

         lx = x[newIDs[j]] - x[newIDs[i]];
         ly = y[newIDs[j]] - y[newIDs[i]];
         linkMag = magnitude( lx, ly );

         cosTheta = ( lx * refx + ly * refy ) / (refMag * linkMag);
         if (cosTheta > cosThetaMax)
         {
            cosThetaMax = cosTheta;
            jID = j;
         } 
           
      } // end loop over j

      // we have found the minimum angle and the corresponding local vertex id.
      // swap ids
      int swapID = newIDs[ nextVertexID ];
      newIDs[ nextVertexID ] = newIDs[jID];
      newIDs[jID] = swapID;

   } // end loop over i

   // reorder x and y coordinate arrays based on newIDs id-array
   real xtemp[ numPoints ];
   real ytemp[ numPoints ];
   for (int i=0; i<numPoints; ++i)
   {
      xtemp[i] = x[i];
      ytemp[i] = y[i];
   }

   for (int i=0; i<numPoints; ++i)
   {
      x[i] = xtemp[ newIDs[i] ];
      y[i] = ytemp[ newIDs[i] ];
   }

   return;

} // end PolyReorder()

//------------------------------------------------------------------------------
void PolyReverse( real* const RESTRICT x, real* const RESTRICT y, const int numPoints )
{
   real xtemp[ numPoints ];
   real ytemp[ numPoints ];
   for (int i=0; i<numPoints; ++i)
   {
      xtemp[i] = x[i];
      ytemp[i] = y[i];
   }

   int k=1;
   for (int i=(numPoints-1); i>0; --i)
   {
      x[k] = xtemp[i];
      y[k] = ytemp[i];
      ++k;
   }
}

//------------------------------------------------------------------------------
void PolyReorderWithNormal( real* const RESTRICT x,
                            real* const RESTRICT y,
                            real* const RESTRICT z,
                            const int numPoints,
                            const real nX,
                            const real nY,
                            const real nZ )
{
   // form link vectors between second and first vertex and third and first
   // vertex
   real lv10X = x[1] - x[0];
   real lv10Y = y[1] - y[0];
   real lv10Z = z[1] - z[0];

   real lv20X = x[2] - x[0];
   real lv20Y = y[2] - y[0];
   real lv20Z = z[2] - z[0];

   // take the cross product of the vectors to get the normal
   real pNrmlX, pNrmlY, pNrmlZ;
   crossProd( lv10X, lv10Y, lv10Z,
              lv20X, lv20Y, lv20Z,
              pNrmlX, pNrmlY, pNrmlZ );

   // dot the computed plane normal based on vertex ordering with the 
   // input normal
   real v = dotProd( pNrmlX, pNrmlY, pNrmlZ, nX, nY, nZ ); 

   // check to see if v is negative. If so, reorient the vertices
   if (v < 0)
   {
      real xTemp[ numPoints ];
      real yTemp[ numPoints ];
      real zTemp[ numPoints ];

      xTemp[0] = x[0];
      yTemp[0] = y[0];
      zTemp[0] = z[0];

      for (int i=1; i<numPoints; ++i)
      {
         xTemp[i] = x[ numPoints - i ];
         yTemp[i] = y[ numPoints - i ];
         zTemp[i] = z[ numPoints - i ];
      }

      for (int i=0; i<numPoints; ++i)
      {
         x[i] = xTemp[i];
         y[i] = yTemp[i];
         z[i] = zTemp[i];
      }
   }

   return;

} // end PolyReorderWithNormal()

//------------------------------------------------------------------------------
bool LinePlaneIntersection( const real xA, const real yA, const real zA,
                            const real xB, const real yB, const real zB,
                            const real xP, const real yP, const real zP,
                            const real nX, const real nY, const real nZ,
                            real& x, real& y, real& z, bool& inPlane )
{

   // compute segment vector
   real lambdaX = xB - xA;
   real lambdaY = yB - yA;
   real lambdaZ = zB - zA;

   // check dot product with plane normal
   real prod = lambdaX*nX + lambdaY*nY + lambdaZ*nZ;

   if (prod == 0.) // line lies in plane
   {
      x = 0.;
      y = 0.;
      z = 0.;
      inPlane = true;
      return false;
   } 

   // compute vector difference between point on plane 
   // and first vertex on segment
   real vX = xP - xA;
   real vY = yP - yA;
   real vZ = zP - zA;

   // compute dot product between <vX, vY, vZ> and the plane normal
   real prodV = vX * nX + vY * nY + vZ * nZ;

   // compute the line segment parameter, t, and check to see if it is 
   // between 0 and 1, inclusive
   real t = prodV / prod;

   if (t >= 0 && t <= 1)
   {
      x = xA + lambdaX * t;
      y = yA + lambdaY * t;
      z = zA + lambdaZ * t;
      inPlane = false;
      return true;
   }
   else 
   {
      x = 0.;
      y = 0.;
      z = 0.;
      inPlane = false;
      return false;
   }

} // end LinePlaneIntersection()

//------------------------------------------------------------------------------
bool PlanePlaneIntersection( const real x1, const real y1, const real z1,
                             const real x2, const real y2, const real z2,
                             const real nX1, const real nY1, const real nZ1,
                             const real nX2, const real nY2, const real nZ2,
                             real& x, real& y, real& z )
{

   // note: this routine has not been tested

   // check dot product between two normals for coplanarity
   real coProd = nX1 * nX2 + nY1 * nY2 + nZ1 * nZ2;

   if (axom::utilities::isNearlyEqual(coProd, 1.0, 1.e-8)) 
   {
      x = 0.;
      y = 0.;
      z = 0.;
      return false;
   }

   // compute dot products between each plane's reference point and the normal
   real prod1 = nX1 * x1 + nY1 * y1 + nZ1 * z1;
   real prod2 = nX2 * x2 + nY2 * y2 + nZ2 * z2; 

   // form matrix of dot products between normals
   real A11 = nX1 * nX1 + nY1 * nY1 + nZ1 * nZ1;
   real A12 = nX1 * nX2 + nY1 * nY2 + nZ1 * nZ2;
   real A22 = nX2 * nX2 + nY2 * nY2 + nZ2 * nZ2;

   // form determinant and inverse determinant of 2x2 matrix
   real detA = A11 * A22 - A12 * A12;
   real invDetA = 1.0 / detA;

   // form inverse matrix components
   real invA11 =  A22;
   real invA12 = -A12;
   real invA22 =  A11;

   // compute two parameters for point on line of intersection
   real s1 = invDetA * (prod1 * invA11 + prod2 * invA12);
   real s2 = invDetA * (prod1 * invA12 + prod2 * invA22);

   // compute the point on the line of intersection
   x = s1 * nX1 + s2 * nX2;
   y = s1 * nY1 + s2 * nY2;
   z = s1 * nZ1 + s2 * nZ2;

   return true;

} // end PlanePlaneIntersection()

//------------------------------------------------------------------------------
void Vertex2DOrderToCCW( const real* const RESTRICT x, const real* const RESTRICT y,
                         real* RESTRICT xTemp, real* RESTRICT yTemp,
                         const int numVert )
{
   if (numVert <= 0)
   {
      SLIC_DEBUG("Vertex2DOrderToCCW: numVert <= 0; returning.");
      return;
   }

   SLIC_ERROR_IF(x == nullptr || y == nullptr || xTemp == nullptr || yTemp == nullptr,
                 "Vertex2DOrderToCCW: must set pointers prior to call to routine.");

   xTemp[0] = x[0];
   yTemp[0] = y[0];

  int k = 1;
  for (int i=numVert; i>0; --i)
  {
     xTemp[k] = x[i];
     yTemp[k] = y[i];
     ++k;
  }

  return;

} // end Vertex2DOrderToCCW()

//------------------------------------------------------------------------------

} // end namespace tribol
