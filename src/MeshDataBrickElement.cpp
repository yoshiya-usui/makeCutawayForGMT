//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------
#include <stddef.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <map>
#include <cmath>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "MeshDataBrickElement.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
MeshDataBrickElement::MeshDataBrickElement():
	m_numElemX(NULL),
	m_numElemY(NULL),
	m_numElemZ(NULL),
	m_numAirLayer(NULL),
	m_edgeLength(NULL)
{
	m_numNodeOneElement = 8;
	m_numEdgeOneElement = 12;
	m_numNodeOnFaceOneElement = 4;
	m_numNeighborElement = 6;

	for ( int i = 0; i < 6; ++i ){
		m_nodesOfElementsBoundaryPlanes[i] = NULL;
	}

}

// Destructer
MeshDataBrickElement::~MeshDataBrickElement(){

	if( m_edgeLength != NULL){
		delete[] m_edgeLength;
		m_edgeLength = NULL;
	}

	for ( int i = 0; i < 6; ++i ){
		m_nodesOfElementsBoundaryPlanes[i] = NULL;
	}

	for ( int i = 0; i < 6; ++i ){
		if( m_nodesOfElementsBoundaryPlanes[i] != NULL){
			delete[] m_nodesOfElementsBoundaryPlanes[i];
			m_nodesOfElementsBoundaryPlanes[i] = NULL;
		}
	}

}

// Copy constructer
MeshDataBrickElement::MeshDataBrickElement(const MeshDataBrickElement& rhs){
	std::cerr << "Error : Copy constructer of the class MeshDataBrickElement is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
MeshDataBrickElement& MeshDataBrickElement::operator=(const MeshDataBrickElement& rhs){
	std::cerr << "Error : Assignment operator of the class MeshDataBrickElement is not implemented." << std::endl;
	exit(1);
}

// Input mesh data from "mesh.dat"
void MeshDataBrickElement::inputMeshData(){

	std::ifstream inFile( "mesh.dat", std::ios::in );
	if( inFile.fail() )
	{
		std::cerr << "File open error : mesh.dat !!" << std::endl;
		exit(1);
	}

	std::string sbuf;
	inFile >> sbuf;

	if( sbuf.substr(0,4).compare("HEXA") != 0 ){
		std::cerr << "Mesh data written in mesh.dat is different from the ones of hexahedral element !!" << std::endl;
		exit(1);
	}

	int nX(0);
	int nY(0);	
	int nZ(0);
	int nAirLayer(0);

	inFile >> nX;
	inFile >> nY;
	inFile >> nZ;
	inFile >> nAirLayer;
	if( nX <= 0 ){
		std::cerr << "Division number of X direction is less than or equal to zero ! : " << nX << std::endl;
		exit(1);
	}
	if( nY <= 0 ){
		std::cerr << "Division number of Y direction is less than or equal to zero ! : " << nY << std::endl;
		exit(1);
	}
	if( nZ <= 0 ){
		std::cerr << "Division number of Z direction is less than or equal to zero ! : " << nZ << std::endl;
		exit(1);
	}
	if( nAirLayer <= 0 ){
		std::cerr << "Number of the air layer is less than or equal to zero ! : " << nAirLayer << std::endl;
		exit(1);
	}
	
	//setMeshDivisons( nX, nY, nZ );// Set mesh division numbers
	m_numElemX = nX;
	m_numElemY = nY;
	m_numElemZ = nZ;
	m_numElemTotal = nX * nY * nZ;
	m_numNodeTotal = (nX+1) * (nY+1) * (nZ+1);
	m_numElemOnBoundaryPlanes[0] = nY * nZ;
	m_numElemOnBoundaryPlanes[1] = nY * nZ;
	m_numElemOnBoundaryPlanes[2] = nZ * nX;
	m_numElemOnBoundaryPlanes[3] = nZ * nX;
	m_numElemOnBoundaryPlanes[4] = nX * nY;
	m_numElemOnBoundaryPlanes[5] = nX * nY;

	//setNumAirLayer( nAirLayer );// Number of air layer
	m_numAirLayer = nAirLayer;

	if( m_xCoordinatesOfNodes != NULL ){
		delete[] m_xCoordinatesOfNodes;	
	}
	m_xCoordinatesOfNodes = new double[m_numNodeTotal];

	if( m_yCoordinatesOfNodes != NULL ){
		delete[] m_yCoordinatesOfNodes;	
	}
	m_yCoordinatesOfNodes = new double[m_numNodeTotal];

	if( m_zCoordinatesOfNodes != NULL ){
		delete[] m_zCoordinatesOfNodes;	
	}
	m_zCoordinatesOfNodes = new double[m_numNodeTotal];

	const int nNode = ( nX + 1 ) * ( nY + 1 ) * ( nZ + 1 );
	for( int iNode = 0; iNode < nNode; ++iNode ){
		int idum(0);
		inFile >> idum;

		// X coordinates of nodes
		double x(0);
		inFile >> x;
		//setXCoordinatesOfNodes( iNode, x );
		m_xCoordinatesOfNodes[iNode] = x;
		
		// Y coordinates of nodes
		double y(0);
		inFile >> y;
		//setYCoordinatesOfNodes( iNode, y );
		m_yCoordinatesOfNodes[iNode] = y;

		// Z coordinates of nodes
		double z(0);
		inFile >> z;
		//setZCoordinatesOfNodes( iNode, z );
		m_zCoordinatesOfNodes[iNode] = z;

#ifdef _DEBUG_WRITE
		std::cout << idum << " " << x << " " << y << " " << z << std::endl; // For debug
#endif

	}

	const int nElem = nX  * nY * nZ;
	if( m_neighborElements != NULL ){
		delete[] m_neighborElements;
	}
	m_neighborElements = new int[ m_numElemTotal * 6 ];

	if( m_nodesOfElements == NULL ){
		delete[] m_nodesOfElements;
	}
	m_nodesOfElements = new int[ m_numElemTotal * 8 ];

	for( int iElem = 0; iElem < nElem; ++iElem ){
		int idum(0);
		inFile >> idum;

		// IDs of neighbor Elements
		for( int i = 0; i < 6; ++i ){
			int neib(0);
			inFile >> neib;
			//setNeighborElements( iElem, i, neib );
			m_neighborElements[ iElem * 6 + i ] = neib;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << neib << std::endl; // For debug
#endif

		}

		// Nodes of the element
		for( int i = 0; i < 8; ++i ){
			int node(0);
			inFile >> node;
			//setNodesOfElements( iElem, i, node );
			m_nodesOfElements[ iElem * 8 + i ] = node;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << node << std::endl; // For debug
#endif

		}

	}
	

	for( int iPlane = 0; iPlane < 6; ++iPlane ){// Loop of boundary planes

		int nElemOnPlane;
		inFile >> nElemOnPlane;

#ifdef _DEBUG_WRITE
		std::cout << nElemOnPlane << std::endl; // For debug
#endif

		if( nElemOnPlane != getNumElemOnBoundaryPlanes( iPlane ) ){
			std::cerr << " Error : Total number of elements belonging to the plane "
				<< iPlane << " is inconsistent with division numbers." << std::endl;
			exit(1);
		}

		if( m_elemBoundaryPlanes[iPlane] != NULL ){
			delete [] m_elemBoundaryPlanes[iPlane];
		}
		m_elemBoundaryPlanes[iPlane] = new int[ nElemOnPlane ];

		if( m_nodesOfElementsBoundaryPlanes[iPlane] != NULL ){
			delete [] m_nodesOfElementsBoundaryPlanes[iPlane];	
		}
		m_nodesOfElementsBoundaryPlanes[iPlane] = new int[ nElemOnPlane * 4 ];

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){		
			// Set elements belonging to the boundary planes
			int elemID(0);
			inFile >> elemID;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << elemID << std::endl; // For debug
#endif

			//setElemBoundaryPlanes( iPlane, iElem, elemID );
			m_elemBoundaryPlanes[iPlane][iElem] = elemID;

			// Set nodes of elements belonging to the boundary planes
			for( int iNode = 0; iNode < 4; ++iNode ){		
				int nodeID(0);
				inFile >> nodeID;

#ifdef _DEBUG_WRITE
				std::cout << iNode << " " << nodeID << std::endl; // For debug
#endif

				//setNodesOfElementsBoundaryPlanes( iPlane, iElem, iNode, nodeID );
				m_nodesOfElementsBoundaryPlanes[iPlane][ iElem * 4 + iNode ] = nodeID;
				
			}
		}
		
	}

	inFile.close();

	//// Output mesh data to VTK file
	//if( OutputFiles::m_vtkFile.is_open() ){
	//	outputMeshDataToVTK();
	//}

}

// Get ID of the nodes of specified element and edge
int MeshDataBrickElement::getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iEdge >= 0 );
	assert( iEdge < 12 );
	assert( num == 0 || num == 1 );

	//if( iEdge < 4 ){
	//	// Edges parallel to X direction
	//	if( iEdge == 0 || iEdge == 2 ){
	//		if( num == 0 ){
	//			return getNodesOfElements( iElem, 2*iEdge );
	//		}else{
	//			return getNodesOfElements( iElem, 2*iEdge+1 );
	//		}
	//	}else{
	//		if( num == 0 ){
	//			return getNodesOfElements( iElem, 2*iEdge+1 );
	//		}else{
	//			return getNodesOfElements( iElem, 2*iEdge );
	//		}
	//	}
	//}else if( iEdge < 8 ){
	//	// Edges parallel to Y direction
	//	if( iEdge == 0 || iEdge == 2 ){
	//		if( num == 0 ){
	//			return getNodesOfElements( iElem, 2*iEdge     );
	//		}else{
	//			return getNodesOfElements( iElem, 2*iEdge + 3 );
	//		}
	//	}else{
	//		if( num == 0 ){
	//			return getNodesOfElements( iElem, 2*(iEdge-1) + 1 );
	//		}else{
	//			return getNodesOfElements( iElem, 2*(iEdge-1) + 2 );
	//		}
	//	}
	//}else{
	//	// Edges parallel to Z direction
	//	if( num == 0 ){
	//		return getNodesOfElements( iElem, iEdge );
	//	}else{
	//		return getNodesOfElements( iElem, iEdge+4 );
	//	}
	//}

	int nodeID[2] = { -1, -1 };
	switch (iEdge){
		case 0:
			nodeID[0] = 0;
			nodeID[1] = 1;
			break;
		case 1:
			nodeID[0] = 3;
			nodeID[1] = 2;
			break;
		case 2:
			nodeID[0] = 4;
			nodeID[1] = 5;
			break;
		case 3:
			nodeID[0] = 7;
			nodeID[1] = 6;
			break;
		case 4:
			nodeID[0] = 0;
			nodeID[1] = 3;
			break;
		case 5:
			nodeID[0] = 4;
			nodeID[1] = 7;
			break;
		case 6:
			nodeID[0] = 1;
			nodeID[1] = 2;
			break;
		case 7:
			nodeID[0] = 5;
			nodeID[1] = 6;
			break;
		case 8:
			nodeID[0] = 0;
			nodeID[1] = 4;
			break;
		case 9:
			nodeID[0] = 1;
			nodeID[1] = 5;
			break;
		case 10:
			nodeID[0] = 3;
			nodeID[1] = 7;
			break;
		case 11:
			nodeID[0] = 2;
			nodeID[1] = 6;
			break;
		default:
			std::cerr << "Edge number : " << iEdge << std::endl;
			exit(1);
			break;
	}	

	return getNodesOfElements( iElem, nodeID[num] );

}

// Get length of the edges parallel to X coordinate
double MeshDataBrickElement::getEdgeLengthX( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node1 = m_nodesOfElements[ iElem * 8 + 1 ];

	return std::fabs( m_xCoordinatesOfNodes[node1] - m_xCoordinatesOfNodes[node0] );

}

// Get length of the edges parallel to Y coordinate
double MeshDataBrickElement::getEdgeLengthY( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node3 = m_nodesOfElements[ iElem * 8 + 3 ];

	return std::fabs( m_yCoordinatesOfNodes[node3] - m_yCoordinatesOfNodes[node0] );

}

// Get length of the edges parallel to Z coordinate
double MeshDataBrickElement::getEdgeLengthZ( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node4 = m_nodesOfElements[ iElem * 8 + 4 ];

	return std::fabs( m_zCoordinatesOfNodes[node4] - m_zCoordinatesOfNodes[node0] );

}

// Get global X coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordX( const int iElem, double localCoordX ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_xCoordinatesOfNodes[node0];
	const double coordMax = m_xCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordX + ( coordMax + coordMin ) * 0.5 ;

}

// Get global Y coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordY( const int iElem, double localCoordY ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_yCoordinatesOfNodes[node0];
	const double coordMax = m_yCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordY + ( coordMax + coordMin ) * 0.5 ;

}

// Get global Z coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordZ( const int iElem, double localCoordZ ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_zCoordinatesOfNodes[node0];
	const double coordMax = m_zCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordZ + ( coordMax + coordMin ) * 0.5 ;

}

// Get number of Elements parallel to X direction
int MeshDataBrickElement::getNumElemX() const{
	return m_numElemX;
}

// Get number of Elements parallel to Y direction
int MeshDataBrickElement::getNumElemY() const{
	return m_numElemY;
}

// Get number of Elements parallel to Z direction
int MeshDataBrickElement::getNumElemZ() const{
	return m_numElemZ;
}

// Get number of the air layer
int MeshDataBrickElement::getNumAirLayer() const{
	return m_numAirLayer;
}

// Calculate number of edges of X-Y plane
int MeshDataBrickElement::calcNumEdgesOnXYPlane() const{
	return m_numElemX * ( m_numElemY + 1 ) + ( m_numElemX + 1 ) * m_numElemY;
}

// Calculate number of edges of Y-Z plane
int MeshDataBrickElement::calcNumEdgesOnYZPlane() const{
	return m_numElemY * ( m_numElemZ + 1 ) + ( m_numElemY + 1 ) * m_numElemZ;
}

// Calculate number of edges of Z-X plane
int MeshDataBrickElement::calcNumEdgesOnZXPlane() const{
	return m_numElemZ * ( m_numElemX + 1 ) + ( m_numElemZ + 1 ) * m_numElemX;
}

// Get array of nodes of elements belonging to the boundary planes
int MeshDataBrickElement::getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const{

	//if( iElem < 0 || iElem >= getNumElemOnBoundaryPlanes( iPlane ) ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getNodesOfElementsBoundaryPlanes. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iNode < 0 || iNode >= m_numNodeOnFaceOneElement ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getNodesOfElementsBoundaryPlanes. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < getNumElemOnBoundaryPlanes( iPlane ) );
	assert( iNode >= 0 );
	assert( iNode < m_numNodeOnFaceOneElement );

	return m_nodesOfElementsBoundaryPlanes[iPlane][ m_numNodeOnFaceOneElement*iElem + iNode ];
}

// Decide whether specified elements share same edges
bool MeshDataBrickElement::shareSameEdges( const int elemID1, const int elemID2 ) const{

	//if( elemID1 < 0 || elemID1 >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : elemID1 is out of range in shareSameNodes. elemID1 = " << elemID1 << std::endl;
	//	exit(1);
	//}

	//if( elemID2 < 0 || elemID2 >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : elemID2 is out of range in shareSameNodes. elemID2 = " << elemID2 << std::endl;
	//	exit(1);
	//}
	assert( elemID1 >= 0 );
	assert( elemID1 < m_numElemTotal );
	assert( elemID2 >= 0 );
	assert( elemID2 < m_numElemTotal );

	// Edges parallel to X direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		const int nodeID1[2] = { getNodesOfElements( elemID1, 2*iEdge1 ), getNodesOfElements( elemID1, 2*iEdge1+1 ) };

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			const int nodeID2[2] = { getNodesOfElements( elemID2, 2*iEdge2 ), getNodesOfElements( elemID2, 2*iEdge2+1 ) };

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	// Edges parallel to Y direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		int nodeID1[2];
		if( iEdge1 == 0 || iEdge1 == 2 ){
			nodeID1[0] = getNodesOfElements( elemID1, 2*iEdge1     );
			nodeID1[1] = getNodesOfElements( elemID1, 2*iEdge1 + 3 );
		}else{
			nodeID1[0] = getNodesOfElements( elemID1, 2*(iEdge1-1) + 1 );
			nodeID1[1] = getNodesOfElements( elemID1, 2*(iEdge1-1) + 2 );
		}

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			int nodeID2[2];
			if( iEdge2 == 0 || iEdge2 == 2 ){
				nodeID2[0] = getNodesOfElements( elemID2, 2*iEdge2     );
				nodeID2[1] = getNodesOfElements( elemID2, 2*iEdge2 + 3 );
			}else{
				nodeID2[0] = getNodesOfElements( elemID2, 2*(iEdge2-1) + 1 );
				nodeID2[1] = getNodesOfElements( elemID2, 2*(iEdge2-1) + 2 );
			}

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	// Edges parallel to Z direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		const int nodeID1[2] = { getNodesOfElements( elemID1, iEdge1 ), getNodesOfElements( elemID1, iEdge1+4 ) };

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			const int nodeID2[2] = { getNodesOfElements( elemID2, iEdge2 ), getNodesOfElements( elemID2, iEdge2+4 ) };

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	return false;

}

// Calculate volume of a specified element
double MeshDataBrickElement::calcVolume( const int elemID ) const{
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	return getEdgeLengthX(elemID) * getEdgeLengthY(elemID) * getEdgeLengthZ(elemID);
}

// Get local coordinate values from coordinate values
// <Input> : iElem, coordX, coordY, coordZ
// <Output>: localCoordX, localCoordY, localCoordZ
void MeshDataBrickElement::getLocalCoordinateValues( const int iElem, const double coordX, const double coordY, const double coordZ,
	double& localCoordX, double& localCoordY, double& localCoordZ ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	coordinateValue coordMin = { m_xCoordinatesOfNodes[node0], m_yCoordinatesOfNodes[node0], m_zCoordinatesOfNodes[node0] };
	coordinateValue coordMax = { m_xCoordinatesOfNodes[node6], m_yCoordinatesOfNodes[node6], m_zCoordinatesOfNodes[node6] };

	localCoordX = 2 * ( coordX - ( coordMin.X + coordMax.X ) * 0.5 ) / std::fabs( coordMax.X - coordMin.X );
	if( localCoordX > 1.0 ){
		localCoordX = 1.0;
	}
	if( localCoordX < -1.0 ){
		localCoordX = -1.0;
	}

	localCoordY = 2 * ( coordY - ( coordMin.Y + coordMax.Y ) * 0.5 ) / std::fabs( coordMax.Y - coordMin.Y );
	if( localCoordY > 1.0 ){
		localCoordY = 1.0;
	}
	if( localCoordY < -1.0 ){
		localCoordY = -1.0;
	}

	localCoordZ = 2 * ( coordZ - ( coordMin.Z + coordMax.Z ) * 0.5 ) / std::fabs( coordMax.Z - coordMin.Z );
	if( localCoordZ > 1.0 ){
		localCoordZ = 1.0;
	}
	if( localCoordZ < -1.0 ){
		localCoordZ = -1.0;
	}

}
