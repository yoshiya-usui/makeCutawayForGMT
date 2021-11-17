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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#include "MeshDataTetraElement.h"
#include "MeshDataBrickElement.h"
#include "MeshDataNonConformingHexaElement.h"
#include "ResistivityBlock.h"

enum PlaneType{
	UNKNOWN = -1,
	ZX_PLANE = 0,
	XY_PLANE,
};

enum ElementType{
	TETRA = 0,
	BRICK,
	NONCONFORMING_HEXA
};

struct Coord3D{
	double X;
	double Y;
	double Z;
};

struct Coord2D{
	double X;
	double Y;
};

ResistivityBlock m_ResistivityBlock;

int m_elementType = 0;
int m_numIteration = 0;
int m_planeType = UNKNOWN;
double m_centerCoord[3];
double m_rotationAngle = 0.0;
double m_normalVector[3];
std::vector<int> m_blockExcluded;

void run( const std::string& paramFile );
void readParameterFile( const std::string& paramFile );
void makeCutaway();
void makeCutawayBrick();
void makeCutawayNonconformingHexa();
bool calcIntersectPoint( const double* coord0, const double* coord1, std::vector<Coord3D>& coordIntersect );
double calcInnerProductWithNormalVec( const double* vec );
void calcNormalVector( const int planeType, const double rotationAngle, double* normalVector );
void deleteSamePoints( std::vector<Coord2D>& vec );
void reorderPoints( std::vector<Coord2D>& vec );
void meterToKilometer( std::vector<Coord2D>& vec );

int main( int argc, char* argv[] ){
	if( argc < 2 ){
		std::cerr << "You must specify parameter file  !!" << std::endl;
		exit(1);
	}
	run( argv[1] );
	return 0;
}

void run( const std::string& paramFile ){

	readParameterFile(paramFile);
	switch (m_elementType){
		case TETRA:
			makeCutaway();
			break;
		case BRICK:
			makeCutawayBrick();
			break;
		case NONCONFORMING_HEXA:
			makeCutawayNonconformingHexa();
			break;
		default:
			std::cerr << "Unknown plane type !!" << std::endl;
			exit(1);
			break;
	}

}

void readParameterFile( const std::string& paramFile ){

	std::ifstream ifs( paramFile.c_str(), std::ios::in );
	if( ifs.fail() ){
		std::cerr << "File open error : " << paramFile.c_str() << " !!" << std::endl;
		exit(1);
	}

	ifs >> m_elementType;
	switch (m_elementType){
		case TETRA:
			std::cout << "Element type : Tetra" << std::endl;
			break;
		case BRICK:
			std::cout << "Element type : Brick" << std::endl;
			break;
		case NONCONFORMING_HEXA:
			std::cout << "Element type : Nonconforming Hexa" << std::endl;
			break;
		default:
			std::cerr << "Unknown plane type !!" << std::endl;
			exit(1);
			break;
	}

	ifs >> m_numIteration;
	std::cout << "Iteration number : " << m_numIteration<< std::endl;

	ifs >> m_planeType;
	switch (m_planeType){
		case ZX_PLANE:
			std::cout << "Plane type : ZX plane" << std::endl;
			break;
		case XY_PLANE:
			std::cout << "Plane type : XY plane" << std::endl;
			break;
		default:
			std::cerr << "Unknown plane type !!" << std::endl;
			exit(1);
			break;
	}

	ifs >> m_centerCoord[0] >> m_centerCoord[1] >> m_centerCoord[2];
	std::cout << "Center coord (km) : (X, Y, Z) = (" << m_centerCoord[0] << ", " << m_centerCoord[1] << ", " << m_centerCoord[2] << ")" << std::endl;
	for( int i = 0; i < 3; ++i ){
		m_centerCoord[i] *= 1000.0;
	}

	ifs >> m_rotationAngle;
	std::cout << "Rotation angle [deg.] : " << m_rotationAngle << std::endl;
	m_rotationAngle *= CommonParameters::deg2rad;

	int numBlockExcepted(0);
	ifs >> numBlockExcepted;
	std::cout << "Number of the parameter cells excluded : " << numBlockExcepted << std::endl;
	m_blockExcluded.reserve(numBlockExcepted);
	for( int iBlk = 0; iBlk < numBlockExcepted; ++iBlk ){
		int blkID(0);
		ifs >> blkID;
		m_blockExcluded.push_back(blkID);
	}

	std::cout << "The following parameter cells are excluded : " << std::endl;
	for( std::vector<int>::iterator itr = m_blockExcluded.begin(); itr != m_blockExcluded.end(); ++itr ){
		std::cout << *itr << std::endl;
	}

	ifs.close();

	calcNormalVector( m_planeType, m_rotationAngle, m_normalVector );

}

void makeCutaway(){

	MeshDataTetraElement m_meshDataTetraElement;

	m_meshDataTetraElement.inputMeshData();
	m_ResistivityBlock.inputResisitivityBlock(m_numIteration);

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}

	ofs.precision(6);

	const int numElem = m_meshDataTetraElement.getNumElemTotal();
	const int numEdgeOnElem = 6;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;

		const int blockID = m_ResistivityBlock.getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = m_meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				m_meshDataTetraElement.getXCoordinatesOfNodes(nodeID0),
				m_meshDataTetraElement.getYCoordinatesOfNodes(nodeID0),
				m_meshDataTetraElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = m_meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				m_meshDataTetraElement.getXCoordinatesOfNodes(nodeID1),
				m_meshDataTetraElement.getYCoordinatesOfNodes(nodeID1),
				m_meshDataTetraElement.getZCoordinatesOfNodes(nodeID1)
			};

			std::vector<Coord3D> coordIntersect3D;
			double vec[2] = { 0.0, 0.0 };
			//if(calcIntersectPoint(coord0, coord1, coordIntersect)){
			calcIntersectPoint(coord0, coord1, coordIntersect3D);
			for( std::vector<Coord3D>::iterator itr = coordIntersect3D.begin(); itr != coordIntersect3D.end(); ++itr ){
				Coord2D coord = { 0.0, 0.0 };
				switch (m_planeType){
					case ZX_PLANE:
						coord.X = m_normalVector[1]*(itr->X - m_centerCoord[0]) - m_normalVector[0]*(itr->Y - m_centerCoord[1]);
						coord.Y = itr->Z;
						break;
					case XY_PLANE:
						vec[0] = itr->Y;
						vec[1] = itr->X;
						coord.X = vec[0] * cos( - m_rotationAngle ) - vec[1] * sin( - m_rotationAngle );
						coord.Y = vec[0] * sin( - m_rotationAngle ) + vec[1] * cos( - m_rotationAngle );
						break;
					default:
						std::cerr << "Unknown plane type !!" << std::endl;
						exit(1);
						break;
				}
				coordIntersectVec.push_back(coord);
			}
		}

		if( static_cast<int>(coordIntersectVec.size()) >= 3 ){
			deleteSamePoints(coordIntersectVec);
			reorderPoints(coordIntersectVec);
			meterToKilometer(coordIntersectVec);
			ofs << "> -Z " << std::setw(15) << std::scientific << log10(m_ResistivityBlock.getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

void makeCutawayBrick(){

	MeshDataBrickElement m_meshDataBrickElement;

	m_meshDataBrickElement.inputMeshData();
	m_ResistivityBlock.inputResisitivityBlock(m_numIteration);

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}

	ofs.precision(6);

	const int numElem = m_meshDataBrickElement.getNumElemTotal();
	const int numEdgeOnElem = 12;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;

		const int blockID = m_ResistivityBlock.getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = m_meshDataBrickElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				m_meshDataBrickElement.getXCoordinatesOfNodes(nodeID0),
				m_meshDataBrickElement.getYCoordinatesOfNodes(nodeID0),
				m_meshDataBrickElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = m_meshDataBrickElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				m_meshDataBrickElement.getXCoordinatesOfNodes(nodeID1),
				m_meshDataBrickElement.getYCoordinatesOfNodes(nodeID1),
				m_meshDataBrickElement.getZCoordinatesOfNodes(nodeID1)
			};

			std::vector<Coord3D> coordIntersect3D;
			double vec[2] = { 0.0, 0.0 };
			calcIntersectPoint(coord0, coord1, coordIntersect3D);
			for( std::vector<Coord3D>::iterator itr = coordIntersect3D.begin(); itr != coordIntersect3D.end(); ++itr ){
				Coord2D coord = { 0.0, 0.0 };
				switch (m_planeType){
					case ZX_PLANE:
						coord.X = m_normalVector[1]*(itr->X - m_centerCoord[0]) - m_normalVector[0]*(itr->Y - m_centerCoord[1]);
						coord.Y = itr->Z;
						break;
					case XY_PLANE:
						vec[0] = itr->Y;
						vec[1] = itr->X;
						coord.X = vec[0] * cos( - m_rotationAngle ) - vec[1] * sin( - m_rotationAngle );
						coord.Y = vec[0] * sin( - m_rotationAngle ) + vec[1] * cos( - m_rotationAngle );
						break;
					default:
						std::cerr << "Unknown plane type !!" << std::endl;
						exit(1);
						break;
				}
				coordIntersectVec.push_back(coord);
			}
		}

		if( static_cast<int>(coordIntersectVec.size()) >= 3 ){
			deleteSamePoints(coordIntersectVec);
			reorderPoints(coordIntersectVec);
			meterToKilometer(coordIntersectVec);
			ofs << "> -Z " << std::setw(15) << std::scientific << log10(m_ResistivityBlock.getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

void makeCutawayNonconformingHexa(){

	MeshDataNonConformingHexaElement m_meshDataNonConformingHexaElement;

	m_meshDataNonConformingHexaElement.inputMeshData();
	m_ResistivityBlock.inputResisitivityBlock(m_numIteration);

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}

	ofs.precision(6);

	const int numElem = m_meshDataNonConformingHexaElement.getNumElemTotal();
	const int numEdgeOnElem = 12;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;

		const int blockID = m_ResistivityBlock.getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = m_meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				m_meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID0),
				m_meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID0),
				m_meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = m_meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				m_meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID1),
				m_meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID1),
				m_meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID1)
			};
			std::vector<Coord3D> coordIntersect3D;
			double vec[2] = { 0.0, 0.0 };
			calcIntersectPoint(coord0, coord1, coordIntersect3D);
			for( std::vector<Coord3D>::iterator itr = coordIntersect3D.begin(); itr != coordIntersect3D.end(); ++itr ){
				Coord2D coord = { 0.0, 0.0 };
				switch (m_planeType){
					case ZX_PLANE:
						coord.X = m_normalVector[1]*(itr->X - m_centerCoord[0]) - m_normalVector[0]*(itr->Y - m_centerCoord[1]);
						coord.Y = itr->Z;
						break;
					case XY_PLANE:
						vec[0] = itr->Y;
						vec[1] = itr->X;
						coord.X = vec[0] * cos( - m_rotationAngle ) - vec[1] * sin( - m_rotationAngle );
						coord.Y = vec[0] * sin( - m_rotationAngle ) + vec[1] * cos( - m_rotationAngle );
						break;
					default:
						std::cerr << "Unknown plane type !!" << std::endl;
						exit(1);
						break;
				}
				coordIntersectVec.push_back(coord);
			}
		}

		if( static_cast<int>(coordIntersectVec.size()) >= 3 ){
			deleteSamePoints(coordIntersectVec);
			reorderPoints(coordIntersectVec);
			meterToKilometer(coordIntersectVec);
			ofs << "> -Z " << std::setw(15) << std::scientific << log10(m_ResistivityBlock.getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

bool calcIntersectPoint( const double* coord0, const double* coord1, std::vector<Coord3D>& coordIntersect ){

	double vec0[3] = { 0.0, 0.0 ,0.0 };
	double vec1[3] = { 0.0, 0.0 ,0.0 };

	for( int i = 0; i < 3; ++i ){
		vec0[i] = coord0[i] - m_centerCoord[i];
		vec1[i] = coord1[i] - m_centerCoord[i];
	}

	const double innerProduct0 = calcInnerProductWithNormalVec(vec0);
	const double innerProduct1 = calcInnerProductWithNormalVec(vec1);

	if( innerProduct0 * innerProduct1 > 0 ){//This segment does not intersect with the plane
		return false;
	}

	const double EPS = 1.0e-12;
	const double sum = fabs(innerProduct0) + fabs(innerProduct1);
	if( sum > EPS ){
		const double ratio = fabs(innerProduct0) / sum;
		const Coord3D coordAdded = {
			coord0[0] + ratio * ( coord1[0] - coord0[0] ),
			coord0[1] + ratio * ( coord1[1] - coord0[1] ),
			coord0[2] + ratio * ( coord1[2] - coord0[2] )
		};
		coordIntersect.push_back(coordAdded);
	}
	else{
		const Coord3D coordAdded0 = { coord0[0], coord0[1], coord0[2] };
		coordIntersect.push_back(coordAdded0);
		const Coord3D coordAdded1 = { coord0[1], coord1[1], coord1[2] };
		coordIntersect.push_back(coordAdded1);
	}

	return true;

}

double calcInnerProductWithNormalVec( const double* vec ){
	double val(0.0);
	for( int i = 0; i < 3; ++i ){
		val += vec[i]*m_normalVector[i];
	}
	return val;
}

void calcNormalVector( const int planeType, const double rotationAngle, double* vec ){

	double vecOrg[3] = {0.0, 0.0, 0.0};	
	switch (planeType){
		case ZX_PLANE:
			vecOrg[0] = 0.0;
			vecOrg[1] = 1.0;
			vecOrg[2] = 0.0;
			vec[0] = vecOrg[0] * cos( - rotationAngle ) - vecOrg[1] * sin( - rotationAngle );
			vec[1] = vecOrg[0] * sin( - rotationAngle ) + vecOrg[1] * cos( - rotationAngle );
			vec[2] = vecOrg[2];
			return;
			break;
		case XY_PLANE:
			vec[0] =  0.0;
			vec[1] =  0.0;
			vec[2] = -1.0;
			return;
			break;
		default:
			std::cerr << "Unknown plane type !!" << std::endl;
			exit(1);
			break;
	}	
	
	return;
}

void deleteSamePoints( std::vector<Coord2D>& vec ){

	std::vector< std::pair<double, double> > vecPair;
	for( std::vector<Coord2D>::iterator itr = vec.begin(); itr != vec.end(); ++itr ){
		vecPair.push_back( std::make_pair(itr->X, itr->Y) );
	}

	sort(vecPair.begin(), vecPair.end());
	vecPair.erase( unique(vecPair.begin(), vecPair.end()), vecPair.end() );

	std::vector<Coord2D> vecNew;
	for( std::vector< std::pair<double, double> >::iterator itr = vecPair.begin(); itr != vecPair.end(); ++itr ){
		Coord2D coord = { itr->first, itr->second };
		vecNew.push_back(coord);
	}

	vec.swap(vecNew);

}

void reorderPoints( std::vector<Coord2D>& vec ){

	assert(!vec.empty());

	Coord2D center = { 0.0, 0.0 };
	for( std::vector<Coord2D>::const_iterator itr = vec.begin(); itr != vec.end(); ++itr ){
		center.X += itr->X;
		center.Y += itr->Y;
	}
	center.X /= static_cast<double>(vec.size());
	center.Y /= static_cast<double>(vec.size());
	
	std::vector< std::pair<double,int> > ang2serial;
	ang2serial.reserve(vec.size());
	int icount(0);
	for( std::vector<Coord2D>::const_iterator itr = vec.begin(); itr != vec.end(); ++itr, ++icount ){
		const double angle = atan2(itr->Y - center.Y, itr->X - center.X);
		ang2serial.push_back( std::make_pair(angle, icount) );
	}
	std::sort(ang2serial.begin(), ang2serial.end());

	std::vector<Coord2D> vecNew;
	vecNew.reserve(vec.size());
	for(std::vector< std::pair<double,int> >::const_iterator itr = ang2serial.begin(); itr != ang2serial.end(); ++itr ){		
		vecNew.push_back( vec[itr->second] );
	}

	vec.swap(vecNew);

}

void meterToKilometer( std::vector<Coord2D>& vec ){

	for( std::vector<Coord2D>::iterator itr = vec.begin(); itr != vec.end(); ++itr ){
		itr->X *= 0.001;
		itr->Y *= 0.001;
	}

}