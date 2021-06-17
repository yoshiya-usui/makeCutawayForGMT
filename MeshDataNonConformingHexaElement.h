//--------------------------------------------------------------------------
// This file is part of makeCutawayForGMT.
//
// makeCutawayForGMT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// makeCutawayForGMT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with makeCutawayForGMT. If not, see <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------

#ifndef DBLDEF_MESHDATA_NONCONFORMING_HEXA_ELEMENT
#define DBLDEF_MESHDATA_NONCONFORMING_HEXA_ELEMENT

#include <vector>
#include "MeshData.h"

// Class of FEM mesh for brick element
class MeshDataNonConformingHexaElement : public MeshData {

public:

	// Constructer
	MeshDataNonConformingHexaElement();

	// Destructer
	virtual ~MeshDataNonConformingHexaElement();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData();

	// Find element including a point
	int findElementIncludingPoint( const double locX, const double locY, const double locZ, double& xi, double& eta, double& zeta ) const;

	// Find elements including point on the surface of the earth
	int findElementIncludingPointOnSurface( const double locX, const double locY, int& faceID, double& xi, double& eta, double& zeta,
		const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod ) const;

	// Find element including a point on the Y-Z plane and return element ID of 2D mesh
	int findElementIncludingPointOnYZPlaneAndReturnElemID2D( const int iPlane, const double locY, const double locZ, double& xi, double& eta ) const;

	// Find element including a point on the Z-X plane and return element ID of 2D mesh
	int findElementIncludingPointOnZXPlaneAndReturnElemID2D( const int iPlane, const double locX, const double locZ, double& xi, double& eta ) const;

	// Get mesh type
	int getMeshType() const;

	// Get ID of a neighbor element
	int getIDOfNeighborElement( const int iElem, const int iFace, const int num ) const;

	// Get number of neighbor elements for an element-face
	int getNumNeighborElement( const int iElem, const int iFace ) const;

	// Get flag specifing whether an element face has slave faces
	bool faceSlaveElements( const int iElem, const int iFace ) const;

	// Get flag specifing whether an element face is outer boundary
	bool isOuterBoundary( const int iElem, const int iFace ) const;

	// Get local face ID of elements belonging to the boundary planes
	int getFaceIDLocalFromElementBoundaryPlanes( const int iPlane, const int iElem ) const;

	// Get global node ID of specified element and edge
	int getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const;

	// Get global node ID of specified element and face
	int getNodeIDGlobalFromElementAndFace( const int iElem, const int iFace, const int num ) const;

	// Get global node ID of specified element belonging to the boundary planes  
	int getNodeIDGlobalFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get global node ID from ID of element belonging to the boundary planes and its edge index
	int getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge, const int num ) const;

	// Get X coordinate of node of specified element belonging to the boundary planes  
	double getCoordXFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Y coordinate of node of specified element belonging to the boundary planes  
	double getCoordYFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Z coordinate of node of specified element belonging to the boundary planes  
	double getCoordZFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get local edge ID from local face ID
	int getEdgeIDLocalFromFaceIDLocal( const int iFace, const int num ) const;

	// Decide whether specified elements share same edges
	virtual bool shareSameEdges( const int elemID1, const int elemID2 ) const;

	// Calculate volume of a specified element
	virtual double calcVolume( const int elemID ) const;
	
	// Get ID of the nodes of elements belonging to the boundary planes
	virtual int getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const;

	// Calculate horizontal coordinate differences of edges of the elements on boundary planes
	double calcHorizontalCoordDifferenceBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Interpolate x coordinate on top or bottom face from local coordinate of horizontal plane
	double calcXCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Interpolate y coordinate on top or bottom face from local coordinate of horizontal plane
	double calcYCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Interpolate z coordinate on top or bottom face from local coordinate of horizontal plane
	double calcZCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Calculate length of edges of elements
	double calcEdgeLengthFromElementAndEdge( const int iElem, const int iEdge ) const;

	// Calculate length of edges of elements on boundary planes
	double calcEdgeLengthFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Get face index of neighbor element
	int getFaceIndexOfNeighborElement( const int iFace ) const;

	// Calculate area of face
	double calcAreaOfFace( const int iElem, const int iFace ) const;

	// Calculate area of face at bottom of mesh
	double calcAreaOfFaceAtBottomOfMesh( const int iElem ) const;

private:

	// Copy constructer
	MeshDataNonConformingHexaElement(const MeshDataNonConformingHexaElement& rhs);

	// Copy assignment operator
	MeshDataNonConformingHexaElement& operator=(const MeshDataNonConformingHexaElement& rhs);

	// Array of IDs of neighbor Elements
	std::vector<int>* m_neighborElementsForNonConformingHexa;

	// Array of faces of elements belonging to the boundary planes
	//   m_facesOfElementsBoundaryPlanes[0] : Y-Z Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[1] : Y-Z Plane ( Plus Side  )
	//   m_facesOfElementsBoundaryPlanes[2] : Z-X Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[3] : Z-X Plane ( Plus Side  )
	//   m_facesOfElementsBoundaryPlanes[4] : X-Y Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[5] : X-Y Plane ( Plus Side  )
	int* m_facesOfElementsBoundaryPlanes[6];

	// Number of elements belonging to the land surface
	int m_numElemOnLandSurface;

	// Array of elements belonging to the land surface
	int* m_elemOnLandSurface;

	// Array of faces belonging to the land surface
	int* m_faceLandSurface;

	// Array converting from face ID to node ID
	int m_faceID2NodeID[6][4];

	// Array converting from face ID to edge ID
	int m_faceID2EdgeID[6][4];

	// Array converting from edge ID to node ID
	int m_edgeID2NodeID[12][2];

	const static int m_numGauss = 2;

	const static int m_numIntegralPoints = m_numGauss * m_numGauss * m_numGauss;

	double m_integralPointXi[m_numIntegralPoints];

	double m_integralPointEta[m_numIntegralPoints];

	double m_integralPointZeta[m_numIntegralPoints];

	double m_weights[m_numIntegralPoints];

	// Array of reference coord xi values for each node
	double m_xiAtNode[8];

	// Array of reference coord eta values for each node
	double m_etaAtNode[8];

	// Array of reference coord zeta values for each node
	double m_zetaAtNode[8];

	// Check whether side element-faces are parallel to Z-X or Y-Z plane
	void checkWhetherSideFaceIsParallelToZXOrYZPlane() const;

	// Check whether the specified point is located in the specified element
	bool isLocatedInTheElement( const double x, const double y, const double z, const int iElem ) const; 

	// Calculate local coordinates
	void calcLocalCoordinates( const int iElem, const double x, const double y, const double z, double& xi, double& eta, double& zeta ) const;

	// Calculate horizontal local coordinates
	void calcHorizontalLocalCoordinates( const int iElem, const double x, const double y, double& xi, double& eta ) const;

	// Calculate determinant of jacobian matrix of the elements
	double calcDeterminantOfJacobianMatrix( const int iElem,  const double xi, const double eta, const double zeta ) const;

};

#endif
