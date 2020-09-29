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
#ifndef DBLDEF_MESHDATA_BRICKELEMENT
#define DBLDEF_MESHDATA_BRICKELEMENT

#include <vector>
#include "MeshData.h"

// Class of FEM mesh for brick element
class MeshDataBrickElement : public MeshData {

public:

	// Constructer
	MeshDataBrickElement();

	// Destructer
	virtual ~MeshDataBrickElement();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData();

	// Get global node ID of specified element and edge
	int getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const;

	// Get length of the edges parallel to X coordinate
	double getEdgeLengthX( const int iElem ) const;

	// Get length of the edges parallel to Y coordinate
	double getEdgeLengthY( const int iElem ) const;

	// Get length of the edges parallel to Z coordinate
	double getEdgeLengthZ( const int iElem ) const;

	// Get global X coordinate from local coordinate
	double calcGlobalCoordX( const int iElem, double localCoordX ) const;

	// Get global Y coordinate from local coordinate
	double calcGlobalCoordY( const int iElem, double localCoordY ) const;

	// Get global Z coordinate from local coordinate
	double calcGlobalCoordZ( const int iElem, double localCoordZ ) const;

	// Get number of Elements parallel to X direction
	int getNumElemX() const;

	// Get number of Elements parallel to Y direction
	int getNumElemY() const;

	// Get number of Elements parallel to Z direction
	int getNumElemZ() const;

	// Get number of the air layer
	int getNumAirLayer() const;	

	// Calculate number of edges of X-Y plane
	int calcNumEdgesOnXYPlane() const;

	// Calculate number of edges of Y-Z plane
	int calcNumEdgesOnYZPlane() const;

	// Calculate number of edges of Z-X plane
	int calcNumEdgesOnZXPlane() const;

	// Get ID of the nodes of elements belonging to the boundary planes
	virtual int getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const;

	// Decide whether specified elements share same edges
	virtual bool shareSameEdges( const int elemID1, const int elemID2 ) const;

	// Calculate volume of a specified element
	virtual double calcVolume( const int elemID ) const;

private:

	// Copy constructer
	MeshDataBrickElement(const MeshDataBrickElement& rhs);

	// Copy assignment operator
	MeshDataBrickElement& operator=(const MeshDataBrickElement& rhs);

	// Number of Elements parallel to X direction
	int m_numElemX;

	// Number of Elements parallel to Y direction
	int m_numElemY;

	// Number of Elements parallel to Z direction
	int m_numElemZ;

	// Number of the air layer
	int m_numAirLayer;

	// Array of edge lenghts
	double* m_edgeLength;

	// Array of nodes of elements belonging to the boundary planes
	//   m_nodesOfElementsBoundaryPlanes[0] : Y-Z Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[1] : Y-Z Plane ( Plus Side  )
	//   m_nodesOfElementsBoundaryPlanes[2] : Z-X Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[3] : Z-X Plane ( Plus Side  )
	//   m_nodesOfElementsBoundaryPlanes[4] : X-Y Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[5] : X-Y Plane ( Plus Side  )
	int* m_nodesOfElementsBoundaryPlanes[6];

	// Get local coordinate values from coordinate values
	virtual void getLocalCoordinateValues( const int iElem, const double coordX, const double coordY, const double coordZ,
		double& localCoordX, double& localCoordY, double& localCoordZ ) const;

};

#endif
