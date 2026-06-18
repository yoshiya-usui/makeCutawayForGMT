//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
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
//-------------------------------------------------------------------------------------------------------
#include "ResistivityBlockIsotropic.h"

#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

// Constructer
ResistivityBlockIsotropic::ResistivityBlockIsotropic():
	m_numResistivityBlockNotFixed(0),
	m_resistivityValues(NULL),
	m_resistivityValuesPre(NULL),
	m_resistivityValuesUpdatedFull(NULL),
	m_weightingConstants(NULL),
	m_fixResistivityValues(NULL),
	m_isolated(NULL)
{}

// Destructer
ResistivityBlockIsotropic::~ResistivityBlockIsotropic(){

	if( m_resistivityValues != NULL){
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}

	if( m_resistivityValuesPre != NULL){
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}

	if( m_resistivityValuesUpdatedFull != NULL){
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}

	if( m_weightingConstants != NULL){
		delete[] m_weightingConstants;
		m_weightingConstants = NULL;
	}

	if( m_fixResistivityValues != NULL){
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}

	if( m_isolated != NULL){
		delete[] m_isolated;
		m_isolated = NULL;
	}

}

// Read isotropic reslstivity values from input file
void ResistivityBlockIsotropic::inputResistivityValues(const int nElem, std::ifstream& inFile){

	if (m_blockID2modelID != NULL) {
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}
	m_blockID2modelID = new int[m_numResistivityBlockTotal];

	if (m_resistivityValues != NULL) {
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}
	m_resistivityValues = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesPre != NULL) {
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}
	m_resistivityValuesPre = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesUpdatedFull != NULL) {
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}
	m_resistivityValuesUpdatedFull = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesMin != NULL) {
		delete[] m_resistivityValuesMin;
		m_resistivityValuesMin = NULL;
	}
	m_resistivityValuesMin = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesMax != NULL) {
		delete[] m_resistivityValuesMax;
		m_resistivityValuesMax = NULL;
	}
	m_resistivityValuesMax = new double[m_numResistivityBlockTotal];

	if (m_weightingConstants != NULL) {
		delete[] m_weightingConstants;
		m_weightingConstants = NULL;
	}
	m_weightingConstants = new double[m_numResistivityBlockTotal];

	if (m_fixResistivityValues != NULL) {
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}
	m_fixResistivityValues = new bool[m_numResistivityBlockTotal];

	if (m_isolated != NULL) {
		delete[] m_isolated;
		m_isolated = NULL;
	}
	m_isolated = new bool[m_numResistivityBlockTotal];

	for (int i = 0; i < m_numResistivityBlockTotal; ++i) {
		m_resistivityValues[i] = 0.0;
		m_resistivityValuesPre[i] = 0.0;
		m_resistivityValuesUpdatedFull[i] = 0.0;
		m_resistivityValuesMin[i] = 0.0;
		m_resistivityValuesMax[i] = 0.0;
		m_weightingConstants[i] = 1.0;
		m_fixResistivityValues[i] = false;
		m_isolated[i] = false;
	}

	m_numResistivityBlockNotFixed = 0;
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		int idum(0);
		int itype(0);
		inFile >> idum >> m_resistivityValues[iBlk] >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk] >> m_weightingConstants[iBlk] >> itype;
		if (idum != iBlk) {
			std::cerr << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		if (m_resistivityValues[iBlk] <= 0.0) {
			std::cerr << "Error : Resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValues[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMin[iBlk] <= 0.0) {
			std::cerr << "Error : Minimum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMin[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMax[iBlk] <= 0.0) {
			std::cerr << "Error : Maximum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMax[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMax[iBlk] < m_resistivityValues[iBlk]) {
			std::cerr << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMin[iBlk] > m_resistivityValues[iBlk]) {
			std::cerr << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
		if (m_weightingConstants[iBlk] <= 0.0) {
			std::cerr << "Error : Weighting constant of block " << iBlk << " is less than or equal to zero !! : " << m_weightingConstants[iBlk] << std::endl;
			exit(1);
		}
		switch (itype) {
		case ResistivityBlockIsotropic::FREE_AND_CONSTRAINED:// Go through
		case ResistivityBlockIsotropic::FREE_AND_ISOLATED:
			m_blockID2modelID[iBlk] = m_numResistivityBlockNotFixed++;
			break;
		case ResistivityBlockIsotropic::FIXED_AND_ISOLATED:// Go through
		case ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED:
			m_fixResistivityValues[iBlk] = true;
			m_blockID2modelID[iBlk] = -1;
			break;
		default:
			std::cerr << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
			exit(1);
			break;
		}
		switch (itype) {
		case ResistivityBlockIsotropic::FREE_AND_CONSTRAINED:// Go through
		case ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED:
			m_isolated[iBlk] = false;
			break;
		case ResistivityBlockIsotropic::FIXED_AND_ISOLATED:// Go through
		case ResistivityBlockIsotropic::FREE_AND_ISOLATED:
			m_isolated[iBlk] = true;
			break;
		default:
			std::cerr << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
			exit(1);
			break;
		}
	}

	if (!m_fixResistivityValues[0]) {
		std::cerr << "Error : Resistivity block 0 must be the air. And, its resistivity must be fixed." << std::endl;
		exit(1);
	}

	inFile.close();

	memcpy(m_resistivityValuesPre, m_resistivityValues, sizeof(double) * (m_numResistivityBlockTotal));

	if (m_modelID2blockID != NULL) {
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}
	m_modelID2blockID = new int[m_numResistivityBlockNotFixed];

	int icount(0);
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {

		if (!m_fixResistivityValues[iBlk]) {
			m_modelID2blockID[icount] = iBlk;
			++icount;
		}

	}

	if (icount != m_numResistivityBlockNotFixed) {
		std::cerr << "Error : icount is not equal to m_numResistivityBlockNotFixed. icount = " << icount << " m_numResistivityBlockNotFixed = " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}
	if (m_numResistivityBlockNotFixed <= 0) {
		std::cerr << "Error : Total number of modifiable resisitivity value is zero or negative !! : " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

}

// Get number of unfixed resistivity parameters
int ResistivityBlockIsotropic::getNumberOfUnfixedResistivityParameters() const {
	return m_numResistivityBlockNotFixed;
}

// Get resistivity values from resisitivity block ID
double ResistivityBlockIsotropic::getResistivityValuesFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValues[ iblk ];
}

// Get previous resistivity values from resisitivity block ID
double ResistivityBlockIsotropic::getResistivityValuesPreFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValuesPre[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValuesPre[ iblk ];

}

// Get conductivity values from resisitivity block ID
double ResistivityBlockIsotropic::getConductivityValuesFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 0.0;
	}

	return 1.0/m_resistivityValues[ iblk ];
}

// Get resistivity values from element ID
double ResistivityBlockIsotropic::getResistivityValuesFromElemID( const int ielem ) const{

	const int iblk = getBlockIDFromElemID(ielem);
	return getResistivityValuesFromBlockID( iblk );
}

// Get conductivity values from element ID
double ResistivityBlockIsotropic::getConductivityValuesFromElemID( const int ielem ) const{

	const int iblk = getBlockIDFromElemID(ielem);
	return getConductivityValuesFromBlockID( iblk );
}

// Get model ID from block ID
int ResistivityBlockIsotropic::getModelIDFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_blockID2modelID[ iblk ];
}

// Get block ID from model ID
int ResistivityBlockIsotropic::getBlockIDFromModelID( const int imdl ) const{

	assert( imdl >= 0 );
	assert( imdl < m_numResistivityBlockNotFixed );

	return m_modelID2blockID[ imdl ];

}

// Get flag specifing whether resistivity value of each block is fixed or not
bool ResistivityBlockIsotropic::isFixedResistivityValue( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_fixResistivityValues[iblk];
}

// Get flag specifing whether resistivity block is excluded from roughing matrix
bool ResistivityBlockIsotropic::isolated( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_isolated[iblk];

}

// Get type of resistivity block
int ResistivityBlockIsotropic::getTypeOfResistivityBlock( const bool fixed, const bool isolated ) const{
	if( fixed ){
		if( isolated ){
			return ResistivityBlockIsotropic::FIXED_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED;
		}
	}
	else{// Free
		if( isolated ){
			return ResistivityBlockIsotropic::FREE_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlockIsotropic::FREE_AND_CONSTRAINED;
		}
	}
}
