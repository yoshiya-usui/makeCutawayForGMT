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
#include "ResistivityBlock.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <assert.h>

// Constructer
ResistivityBlock::ResistivityBlock():
	m_elementID2blockID(NULL),
	m_blockID2modelID(NULL),
	m_modelID2blockID(NULL),
	m_numResistivityBlockTotal(0),
	m_numResistivityBlockNotFixed(0),
	m_resistivityValues(NULL),
	m_resistivityValuesPre(NULL),
	m_resistivityValuesUpdatedFull(NULL),
	m_fixResistivityValues(NULL),
	m_blockID2Elements(NULL)
{}

// Destructer
ResistivityBlock::~ResistivityBlock(){

	if( m_elementID2blockID != NULL){
		delete[] m_elementID2blockID;
		m_elementID2blockID = NULL;
	}

	if( m_blockID2modelID != NULL){
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}

	if( m_modelID2blockID != NULL){
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}

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

	if( m_fixResistivityValues != NULL){
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}

	if( m_blockID2Elements != NULL){
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}

}

// Read data of resisitivity block model from input file
void ResistivityBlock::inputResisitivityBlock(const int iterNum){

	std::ostringstream inputFile;
	inputFile << "resistivity_block_iter" << iterNum << ".dat";
	std::ifstream inFile( inputFile.str().c_str(), std::ios::in );

	if( inFile.fail() )
	{
		std::cerr << "File open error : " << inputFile.str().c_str() << " !!" << std::endl;
		exit(1);
	}

	int nElem(0);
	inFile >> nElem;
	if( m_elementID2blockID != NULL ){
		delete[] m_elementID2blockID;
		m_elementID2blockID = NULL;
	}
	m_elementID2blockID = new int[ nElem ];

	int nBlk(0);
	inFile >> nBlk;	
	m_numResistivityBlockTotal = nBlk;

	if( m_blockID2modelID != NULL){
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}
	m_blockID2modelID = new int[ m_numResistivityBlockTotal ];

	if( m_resistivityValues != NULL ){
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}
	m_resistivityValues = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesPre != NULL){
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}
	m_resistivityValuesPre = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesUpdatedFull != NULL){
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}
	m_resistivityValuesUpdatedFull = new double[ m_numResistivityBlockTotal ];

	if( m_fixResistivityValues != NULL ){
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}
	m_fixResistivityValues = new bool[ m_numResistivityBlockTotal ];

	for( int i = 0; i < m_numResistivityBlockTotal; ++i ){
		m_resistivityValues[i] = 0.0;
		m_resistivityValuesPre[i] = 0.0;
		m_resistivityValuesUpdatedFull[i] = 0.0;
		m_fixResistivityValues[i] = false;
	}

#ifdef _DEBUG_WRITE
	std::cout << nElem << " " << m_numResistivityBlockTotal << std::endl; // For debug
#endif
	
	for( int iElem = 0; iElem < nElem; ++iElem ){
		int idum(0);
		int iblk(0);// Resistivity block ID
		inFile >> idum >> iblk;
		m_elementID2blockID[ iElem ] = iblk;

		if( iblk >= m_numResistivityBlockTotal || iblk < 0 )	{
			std::cerr << "Error : Resistivity block ID " << iblk << " of element " << iElem << " is improper !!" << std::endl;
			exit(1);
		}		

#ifdef _DEBUG_WRITE
		std::cout << iElem << " " << iblk << std::endl; // For debug
#endif
	}
	
	m_numResistivityBlockNotFixed = 0;
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		int idum(0);
		int ifix(0);
#ifdef _OLD
		inFile >> idum >> m_resistivityValues[iBlk] >> ifix;
#else
		double ddum(0.0);
		inFile >> idum >> m_resistivityValues[iBlk] >> ddum >> ddum >> ddum >> ifix;
#endif

		if( idum != iBlk ){
			std::cerr << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}

		if( ifix == 1 ){
			m_fixResistivityValues[iBlk] = true;
			m_blockID2modelID[iBlk] = -1;
		}else{
			m_blockID2modelID[iBlk] = m_numResistivityBlockNotFixed++;
		}

#ifdef _DEBUG_WRITE
		std::cout << iBlk << " " << m_resistivityValues[iBlk] << " " << m_fixResistivityValues[iBlk] << " " << m_blockID2modelID[iBlk] << std::endl; // For debug
#endif

	}

	if( !m_fixResistivityValues[0] ){
		std::cerr << "Error : Resistivity block 0 must be the air. And, its resistivity must be fixed." << std::endl;
		exit(1);
	}

	inFile.close();

	memcpy( m_resistivityValuesPre, m_resistivityValues, sizeof(double)*(m_numResistivityBlockTotal) );

	if( m_modelID2blockID != NULL){
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}
	m_modelID2blockID = new int[ m_numResistivityBlockNotFixed ];

	int icount(0);
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){

		if( !m_fixResistivityValues[iBlk] ){
			m_modelID2blockID[icount] = iBlk;
			++icount;
		}

	}

	if( icount != m_numResistivityBlockNotFixed ){
		std::cerr << "Error : icount is not equal to m_numResistivityBlockNotFixed. icount = " << icount << " m_numResistivityBlockNotFixed = " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	for( int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl ){
		std::cout << " iMdl m_modelID2blockID[iMdl] : " << iMdl << " " << m_modelID2blockID[iMdl] << std::endl;
	}
#endif

	if( m_blockID2Elements != NULL){
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}
	m_blockID2Elements= new std::vector< std::pair<int,double> >[m_numResistivityBlockTotal];

	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_blockID2Elements[ m_elementID2blockID[ iElem ] ].push_back( std::make_pair(iElem,1.0) );
	}
#ifndef _LINUX
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		m_blockID2Elements[ iBlk ].shrink_to_fit();
	}
#endif

	if( m_numResistivityBlockNotFixed <= 0 ){
		std::cerr << "Error : Total number of modifiable resisitivity value is zero or negative !! : " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		const int num = static_cast<int>( m_blockID2Elements[ iBlk ].size() ); 
		for( int i = 0; i < num; ++i ){
			std::cout << " m_blockID2Elements[ " << iBlk << " ][ " << i << "] : " << m_blockID2Elements[iBlk][i].first << std::endl;
		}
	}
#endif

}

// Get resistivity values from resisitivity block ID
double ResistivityBlock::getResistivityValuesFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValues[ iblk ];
}

//// Get resisitivity block ID from element ID
//inline int ResistivityBlock::getBlockIDFromElemID( const int ielem ) const{
//	return m_elementID2blockID[ ielem ];
//}

// Get conductivity values from resisitivity block ID
double ResistivityBlock::getConductivityValuesFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 0.0;
	}

	return 1.0/m_resistivityValues[ iblk ];
}

// Get resistivity values from element ID
double ResistivityBlock::getResistivityValuesFromElemID( const int ielem ) const{

	//const int iblk = m_elementID2blockID[ ielem ];
	const int iblk = getBlockIDFromElemID(ielem);
	return getResistivityValuesFromBlockID( iblk );
}

// Get conductivity values from element ID
double ResistivityBlock::getConductivityValuesFromElemID( const int ielem ) const{

	//const int iblk = m_elementID2blockID[ ielem ];
	const int iblk = getBlockIDFromElemID(ielem);
	return getConductivityValuesFromBlockID( iblk );
}

// Get model ID from block ID
int ResistivityBlock::getModelIDFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_blockID2modelID[ iblk ];
}

// Get block ID from model ID
int ResistivityBlock::getBlockIDFromModelID( const int imdl ) const{

	//if( imdl < 0 || imdl >= m_numResistivityBlockNotFixed ){
	//	OutputFiles::m_logFile << "Error : Specified model ID is out of range. imdl = " << imdl << std::endl;
	//	exit(1);
	//}
	assert( imdl >= 0 );
	assert( imdl < m_numResistivityBlockNotFixed );

	return m_modelID2blockID[ imdl ];

}

// Get total number of resistivity blocks
int ResistivityBlock::getNumResistivityBlockTotal() const{
	return m_numResistivityBlockTotal;
}

// Get number of resistivity blocks whose resistivity values are fixed
int ResistivityBlock::getNumResistivityBlockNotFixed() const{
	return m_numResistivityBlockNotFixed;
}

// Get flag specifing whether resistivity value of each block is fixed or not
bool ResistivityBlock::isFixedResistivityValue( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_fixResistivityValues[iblk];
}

// Get arrays of elements belonging to each resistivity block
const std::vector< std::pair<int,double> >&  ResistivityBlock::getBlockID2Elements( const int iBlk ) const{

	return m_blockID2Elements[iBlk];

}
