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
#include "ResistivityBlockIsotropic.h"
#include "ResistivityBlockAnisotropic.h"

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

ResistivityBlock* m_ptrResistivityBlock(NULL);

int m_elementType = 0;
int m_numIteration = 0;
int m_planeType = UNKNOWN;
double m_centerCoord[3];
double m_rotationAngle = 0.0;
double m_normalVector[3];
std::vector<int> m_blockExcluded;
double m_anisotropyThreshold(0.0);

void run(const bool isAnisotropic, const std::string& paramFile);
void readParameterFile(const bool isAnisotropic, const std::string& paramFile);
void makeCutawayTetraIsotropic();
void makeCutawayTetraAnisotropic();
void makeCutawayBrick();
void makeCutawayNonconformingHexaIsotropic();
void makeCutawayNonconformingHexaAnisotropic();
bool calcIntersectPoint( const double* coord0, const double* coord1, std::vector<Coord3D>& coordIntersect );
double calcInnerProductWithNormalVec( const double* vec );
void calcNormalVector( const int planeType, const double rotationAngle, double* normalVector );
void deleteSamePoints( std::vector<Coord2D>& vec );
void reorderPoints( std::vector<Coord2D>& vec );
void meterToKilometer( std::vector<Coord2D>& vec );

int main( int argc, char* argv[] ){

	if (argc < 2) {
		std::cerr << "You must specify parameter file  !!" << std::endl;
		exit(1);
	}
	bool isAnisotropic(false);
	std::string inputFileName;
	for (int iarg = 1; iarg < argc; ++iarg) {
		const std::string sbuf = argv[iarg];
		if (sbuf.substr(0, 6).compare("-aniso") == 0) {
			isAnisotropic = true;
		} else {
			inputFileName = sbuf;
		}
	}
	run(isAnisotropic, inputFileName);
	return 0;

}

void run(const bool isAnisotropic, const std::string& paramFile) {

	readParameterFile(isAnisotropic, paramFile);
	if (isAnisotropic) {
		m_ptrResistivityBlock = new ResistivityBlockAnisotropic();
	}
	else {
		m_ptrResistivityBlock = new ResistivityBlockIsotropic();
	}
	m_ptrResistivityBlock->inputResistivityBlock(m_numIteration);
	switch (m_elementType) {
	case TETRA:
		isAnisotropic ? makeCutawayTetraAnisotropic() :	makeCutawayTetraIsotropic();
		break;
	case BRICK:
		makeCutawayBrick();
		break;
	case NONCONFORMING_HEXA:
		isAnisotropic ? makeCutawayNonconformingHexaAnisotropic() : makeCutawayNonconformingHexaIsotropic();
		break;
	default:
		std::cerr << "Unknown element type !!" << std::endl;
		exit(1);
		break;
	}
	delete m_ptrResistivityBlock;

}

void readParameterFile(const bool isAnisotropic, const std::string& paramFile) {

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

	if (isAnisotropic) {
		ifs >> m_anisotropyThreshold;
		std::cout << "The threshold of anisotropy above which strike and dip are outputted : " << m_anisotropyThreshold << std::endl;
	}

	ifs.close();

	calcNormalVector( m_planeType, m_rotationAngle, m_normalVector );

}

void makeCutawayTetraIsotropic() {

	MeshDataTetraElement meshDataTetraElement;

	meshDataTetraElement.inputMeshData();

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}
	ofs.precision(6);

	const std::string splitter = CommonParameters::isGMT5OrHigher ? "> -Z" : "> -Z ";

	const int numElem = meshDataTetraElement.getNumElemTotal();
	const int numEdgeOnElem = 6;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;
		const int blockID = m_ptrResistivityBlock->getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				meshDataTetraElement.getXCoordinatesOfNodes(nodeID0),
				meshDataTetraElement.getYCoordinatesOfNodes(nodeID0),
				meshDataTetraElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				meshDataTetraElement.getXCoordinatesOfNodes(nodeID1),
				meshDataTetraElement.getYCoordinatesOfNodes(nodeID1),
				meshDataTetraElement.getZCoordinatesOfNodes(nodeID1)
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
			ofs << splitter << log10(dynamic_cast<ResistivityBlockIsotropic*>(m_ptrResistivityBlock)->getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

void makeCutawayTetraAnisotropic() {

	MeshDataTetraElement meshDataTetraElement;
	meshDataTetraElement.inputMeshData();

	std::ostringstream ofileRhoXX;
	ofileRhoXX << "rhoXX_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoXX(ofileRhoXX.str().c_str(), std::ios::out);
	if (ofsRhoXX.fail()) {
		std::cerr << "File open error : " << ofileRhoXX.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoXX.precision(6);

	std::ostringstream ofileRhoYY;
	ofileRhoYY << "rhoYY_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoYY(ofileRhoYY.str().c_str(), std::ios::out);
	if (ofsRhoYY.fail()) {
		std::cerr << "File open error : " << ofileRhoYY.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoYY.precision(6);

	std::ostringstream ofileRhoZZ;
	ofileRhoZZ << "rhoZZ_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoZZ(ofileRhoZZ.str().c_str(), std::ios::out);
	if (ofsRhoZZ.fail()) {
		std::cerr << "File open error : " << ofileRhoZZ.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoZZ.precision(6);

	std::ostringstream ofileAnisotropy;
	ofileAnisotropy << "anisotropy_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsAnisotropy(ofileAnisotropy.str().c_str(), std::ios::out);
	if (ofsAnisotropy.fail()) {
		std::cerr << "File open error : " << ofileAnisotropy.str() << " !!" << std::endl;
		exit(1);
	}
	ofsAnisotropy.precision(6);

	std::ostringstream ofileStrike;
	ofileStrike << "strike_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsStrike(ofileStrike.str().c_str(), std::ios::out);
	if (ofsStrike.fail()) {
		std::cerr << "File open error : " << ofileStrike.str() << " !!" << std::endl;
		exit(1);
	}
	ofsStrike.precision(6);

	std::ostringstream ofileDip;
	ofileDip << "dip_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsDip(ofileDip.str().c_str(), std::ios::out);
	if (ofsDip.fail()) {
		std::cerr << "File open error : " << ofileDip.str() << " !!" << std::endl;
		exit(1);
	}
	ofsDip.precision(6);

	std::ostringstream ofileSlant;
	ofileSlant << "slant_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsSlant(ofileSlant.str().c_str(), std::ios::out);
	if (ofsSlant.fail()) {
		std::cerr << "File open error : " << ofileSlant.str() << " !!" << std::endl;
		exit(1);
	}
	ofsSlant.precision(6);

	const std::string splitter = CommonParameters::isGMT5OrHigher ? "> -Z" : "> -Z ";

	const int numElem = meshDataTetraElement.getNumElemTotal();
	const int numEdgeOnElem = 6;
	for (int iElem = 0; iElem < numElem; ++iElem) {
		std::vector<Coord2D> coordIntersectVec;
		const int blockID = m_ptrResistivityBlock->getBlockIDFromElemID(iElem);
		if (find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end()) {
			continue;// Excluded
		};

		for (int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge) {
			const int nodeID0 = meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				meshDataTetraElement.getXCoordinatesOfNodes(nodeID0),
				meshDataTetraElement.getYCoordinatesOfNodes(nodeID0),
				meshDataTetraElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = meshDataTetraElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				meshDataTetraElement.getXCoordinatesOfNodes(nodeID1),
				meshDataTetraElement.getYCoordinatesOfNodes(nodeID1),
				meshDataTetraElement.getZCoordinatesOfNodes(nodeID1)
			};

			std::vector<Coord3D> coordIntersect3D;
			double vec[2] = { 0.0, 0.0 };
			//if(calcIntersectPoint(coord0, coord1, coordIntersect)){
			calcIntersectPoint(coord0, coord1, coordIntersect3D);
			for (std::vector<Coord3D>::iterator itr = coordIntersect3D.begin(); itr != coordIntersect3D.end(); ++itr) {
				Coord2D coord = { 0.0, 0.0 };
				switch (m_planeType) {
				case ZX_PLANE:
					coord.X = m_normalVector[1] * (itr->X - m_centerCoord[0]) - m_normalVector[0] * (itr->Y - m_centerCoord[1]);
					coord.Y = itr->Z;
					break;
				case XY_PLANE:
					vec[0] = itr->Y;
					vec[1] = itr->X;
					coord.X = vec[0] * cos(-m_rotationAngle) - vec[1] * sin(-m_rotationAngle);
					coord.Y = vec[0] * sin(-m_rotationAngle) + vec[1] * cos(-m_rotationAngle);
					break;
				default:
					std::cerr << "Unknown plane type !!" << std::endl;
					exit(1);
					break;
				}
				coordIntersectVec.push_back(coord);
			}
		}

		if (static_cast<int>(coordIntersectVec.size()) >= 3) {
			deleteSamePoints(coordIntersectVec);
			reorderPoints(coordIntersectVec);
			meterToKilometer(coordIntersectVec);
			const ResistivityBlockAnisotropic::AnistropicResistivityParameters params = dynamic_cast<ResistivityBlockAnisotropic*>(m_ptrResistivityBlock)->getAnisotropicResistivityParametersFromBlockID(blockID);

			ofsRhoXX << splitter << log10(params.rhoXX) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoXX << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoXX << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			ofsRhoYY << splitter << log10(params.rhoYY) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoYY << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoYY << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			ofsRhoZZ << splitter << log10(params.rhoZZ) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoZZ << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoZZ << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			const double anisotropy = fabs(log10(params.rhoXX) - log10(params.rhoYY));
			ofsAnisotropy << splitter << anisotropy << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsAnisotropy << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsAnisotropy << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			if (anisotropy > m_anisotropyThreshold) {
				ofsStrike << splitter << params.strike * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsStrike << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsStrike << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

				ofsDip << splitter << params.dip * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsDip << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsDip << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

				ofsSlant << splitter << params.slant * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsSlant << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsSlant << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
			}
		}
	}

	ofsRhoXX.close();
	ofsRhoYY.close();
	ofsRhoZZ.close();
	ofsAnisotropy.close();
	ofsStrike.close();
	ofsDip.close();
	ofsSlant.close();

}

void makeCutawayBrick(){

	MeshDataBrickElement meshDataBrickElement;
	meshDataBrickElement.inputMeshData();

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}

	ofs.precision(6);

	const std::string splitter = CommonParameters::isGMT5OrHigher ? "> -Z" : "> -Z ";

	const int numElem = meshDataBrickElement.getNumElemTotal();
	const int numEdgeOnElem = 12;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;
		const int blockID = m_ptrResistivityBlock->getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = meshDataBrickElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				meshDataBrickElement.getXCoordinatesOfNodes(nodeID0),
				meshDataBrickElement.getYCoordinatesOfNodes(nodeID0),
				meshDataBrickElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = meshDataBrickElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				meshDataBrickElement.getXCoordinatesOfNodes(nodeID1),
				meshDataBrickElement.getYCoordinatesOfNodes(nodeID1),
				meshDataBrickElement.getZCoordinatesOfNodes(nodeID1)
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
			ofs << splitter << log10(dynamic_cast<ResistivityBlockIsotropic*>(m_ptrResistivityBlock)->getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

void makeCutawayNonconformingHexaIsotropic() {

	MeshDataNonConformingHexaElement meshDataNonConformingHexaElement;
	meshDataNonConformingHexaElement.inputMeshData();

	std::ostringstream ofile;
	ofile << "resistivity_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofs( ofile.str().c_str(), std::ios::out );
	if( ofs.fail() ){
		std::cerr << "File open error : " << ofile.str() << " !!" << std::endl;
		exit(1);
	}

	ofs.precision(6);

	const std::string splitter = CommonParameters::isGMT5OrHigher ? "> -Z" : "> -Z ";

	const int numElem = meshDataNonConformingHexaElement.getNumElemTotal();
	const int numEdgeOnElem = 12;
	for( int iElem = 0; iElem < numElem; ++iElem ){
		std::vector<Coord2D> coordIntersectVec;
		const int blockID = m_ptrResistivityBlock->getBlockIDFromElemID(iElem);
		if( find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end() ){
			continue;// Excluded
		};

		for( int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge ){
			const int nodeID0 = meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID0),
				meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID0),
				meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID1),
				meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID1),
				meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID1)
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
			ofs << splitter << log10(dynamic_cast<ResistivityBlockIsotropic*>(m_ptrResistivityBlock)->getResistivityValuesFromBlockID(blockID)) << std::endl;
			for( std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr ){
				ofs << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofs << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
		}
	}

	ofs.close();

}

void makeCutawayNonconformingHexaAnisotropic() {

	MeshDataNonConformingHexaElement meshDataNonConformingHexaElement;
	meshDataNonConformingHexaElement.inputMeshData();

	std::ostringstream ofileRhoXX;
	ofileRhoXX << "rhoXX_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoXX(ofileRhoXX.str().c_str(), std::ios::out);
	if (ofsRhoXX.fail()) {
		std::cerr << "File open error : " << ofileRhoXX.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoXX.precision(6);

	std::ostringstream ofileRhoYY;
	ofileRhoYY << "rhoYY_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoYY(ofileRhoYY.str().c_str(), std::ios::out);
	if (ofsRhoYY.fail()) {
		std::cerr << "File open error : " << ofileRhoYY.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoYY.precision(6);

	std::ostringstream ofileRhoZZ;
	ofileRhoZZ << "rhoZZ_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsRhoZZ(ofileRhoZZ.str().c_str(), std::ios::out);
	if (ofsRhoZZ.fail()) {
		std::cerr << "File open error : " << ofileRhoZZ.str() << " !!" << std::endl;
		exit(1);
	}
	ofsRhoZZ.precision(6);

	std::ostringstream ofileAnisotropy;
	ofileAnisotropy << "anisotropy_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsAnisotropy(ofileAnisotropy.str().c_str(), std::ios::out);
	if (ofsAnisotropy.fail()) {
		std::cerr << "File open error : " << ofileAnisotropy.str() << " !!" << std::endl;
		exit(1);
	}
	ofsAnisotropy.precision(6);

	std::ostringstream ofileStrike;
	ofileStrike << "strike_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsStrike(ofileStrike.str().c_str(), std::ios::out);
	if (ofsStrike.fail()) {
		std::cerr << "File open error : " << ofileStrike.str() << " !!" << std::endl;
		exit(1);
	}
	ofsStrike.precision(6);

	std::ostringstream ofileDip;
	ofileDip << "dip_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsDip(ofileDip.str().c_str(), std::ios::out);
	if (ofsDip.fail()) {
		std::cerr << "File open error : " << ofileDip.str() << " !!" << std::endl;
		exit(1);
	}
	ofsDip.precision(6);

	std::ostringstream ofileSlant;
	ofileSlant << "slant_GMT_iter" << m_numIteration << ".dat";
	std::ofstream ofsSlant(ofileSlant.str().c_str(), std::ios::out);
	if (ofsSlant.fail()) {
		std::cerr << "File open error : " << ofileSlant.str() << " !!" << std::endl;
		exit(1);
	}
	ofsSlant.precision(6);

	const std::string splitter = CommonParameters::isGMT5OrHigher ? "> -Z" : "> -Z ";

	const int numElem = meshDataNonConformingHexaElement.getNumElemTotal();
	const int numEdgeOnElem = 12;
	for (int iElem = 0; iElem < numElem; ++iElem) {
		std::vector<Coord2D> coordIntersectVec;
		const int blockID = m_ptrResistivityBlock->getBlockIDFromElemID(iElem);
		if (find(m_blockExcluded.begin(), m_blockExcluded.end(), blockID) != m_blockExcluded.end()) {
			continue;// Excluded
		};

		for (int iEdge = 0; iEdge < numEdgeOnElem; ++iEdge) {
			const int nodeID0 = meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 0);
			const double coord0[3] = {
				meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID0),
				meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID0),
				meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID0)
			};
			const int nodeID1 = meshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge(iElem, iEdge, 1);
			const double coord1[3] = {
				meshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID1),
				meshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID1),
				meshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID1)
			};
			std::vector<Coord3D> coordIntersect3D;
			double vec[2] = { 0.0, 0.0 };
			calcIntersectPoint(coord0, coord1, coordIntersect3D);
			for (std::vector<Coord3D>::iterator itr = coordIntersect3D.begin(); itr != coordIntersect3D.end(); ++itr) {
				Coord2D coord = { 0.0, 0.0 };
				switch (m_planeType) {
				case ZX_PLANE:
					coord.X = m_normalVector[1] * (itr->X - m_centerCoord[0]) - m_normalVector[0] * (itr->Y - m_centerCoord[1]);
					coord.Y = itr->Z;
					break;
				case XY_PLANE:
					vec[0] = itr->Y;
					vec[1] = itr->X;
					coord.X = vec[0] * cos(-m_rotationAngle) - vec[1] * sin(-m_rotationAngle);
					coord.Y = vec[0] * sin(-m_rotationAngle) + vec[1] * cos(-m_rotationAngle);
					break;
				default:
					std::cerr << "Unknown plane type !!" << std::endl;
					exit(1);
					break;
				}
				coordIntersectVec.push_back(coord);
			}
		}

		if (static_cast<int>(coordIntersectVec.size()) >= 3) {
			deleteSamePoints(coordIntersectVec);
			reorderPoints(coordIntersectVec);
			meterToKilometer(coordIntersectVec);
			const ResistivityBlockAnisotropic::AnistropicResistivityParameters params = dynamic_cast<ResistivityBlockAnisotropic*>(m_ptrResistivityBlock)->getAnisotropicResistivityParametersFromBlockID(blockID);

			ofsRhoXX << splitter << log10(params.rhoXX) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoXX << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoXX << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			ofsRhoYY << splitter << log10(params.rhoYY) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoYY << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoYY << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			ofsRhoZZ << splitter << log10(params.rhoZZ) << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsRhoZZ << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsRhoZZ << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			const double anisotropy = fabs(log10(params.rhoXX) - log10(params.rhoYY));
			ofsAnisotropy << splitter << anisotropy << std::endl;
			for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
				ofsAnisotropy << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
			}
			ofsAnisotropy << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

			if (anisotropy > m_anisotropyThreshold) {
				ofsStrike << splitter << params.strike * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsStrike << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsStrike << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

				ofsDip << splitter << params.dip * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsDip << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsDip << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;

				ofsSlant << splitter << params.slant * CommonParameters::rad2deg << std::endl;
				for (std::vector<Coord2D>::const_iterator itr = coordIntersectVec.begin(); itr != coordIntersectVec.end(); ++itr) {
					ofsSlant << std::setw(15) << std::scientific << itr->X << std::setw(15) << std::scientific << itr->Y << std::endl;
				}
				ofsSlant << std::setw(15) << std::scientific << coordIntersectVec.front().X << std::setw(15) << std::scientific << coordIntersectVec.front().Y << std::endl;
			}
		}
	}

	ofsRhoXX.close();
	ofsRhoYY.close();
	ofsAnisotropy.close();
	ofsStrike.close();
	ofsDip.close();

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

	const double EPS = 1.0e-9;
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
		const Coord3D coordAdded1 = { coord1[1], coord1[1], coord1[2] };
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