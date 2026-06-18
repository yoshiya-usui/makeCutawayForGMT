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

#include "ResistivityBlock.h"

#include <string.h>
#include <assert.h>
#include <iomanip>
#include <sstream>
#include <fstream>

// Constructer
ResistivityBlock::ResistivityBlock():
	m_numResistivityBlockTotal(0),
	m_elementID2blockID(NULL),
	m_blockID2Elements(NULL),
	m_modelID2blockID(NULL),
	m_blockID2modelID(NULL),
	m_resistivityValuesMin(NULL),
	m_resistivityValuesMax(NULL)
{}

// Destructer
ResistivityBlock::~ResistivityBlock(){

	if( m_elementID2blockID != NULL){
		delete[] m_elementID2blockID;
		m_elementID2blockID = NULL;
	}

	if (m_blockID2Elements != NULL) {
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}

	if (m_modelID2blockID != NULL) {
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}

	if( m_blockID2modelID != NULL){
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}

	if( m_resistivityValuesMin != NULL){
		delete[] m_resistivityValuesMin;
		m_resistivityValuesMin = NULL;
	}

	if( m_resistivityValuesMax != NULL){
		delete[] m_resistivityValuesMax;
		m_resistivityValuesMax = NULL;
	}

}

// Read data of resisitivity block model from input file
void ResistivityBlock::inputResistivityBlock( const int iterationNumInit ){

	std::ostringstream inputFile;
	inputFile << "resistivity_block_iter" << iterationNumInit << ".dat";
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

	for (int iElem = 0; iElem < nElem; ++iElem) {
		int idum(0);
		int iblk(0);// Resistivity block ID
		inFile >> idum >> iblk;
		m_elementID2blockID[iElem] = iblk;
		if (iblk >= m_numResistivityBlockTotal || iblk < 0) {
			std::cerr << "Error : Resistivity block ID " << iblk << " of element " << iElem << " is improper !!" << std::endl;
			exit(1);
		}
	}

	if (m_blockID2Elements != NULL) {
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}
	m_blockID2Elements = new std::vector<int>[m_numResistivityBlockTotal];
	for (int iElem = 0; iElem < nElem; ++iElem) {
		m_blockID2Elements[m_elementID2blockID[iElem]].push_back(iElem);
	}
#ifndef _LINUX
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		m_blockID2Elements[iBlk].shrink_to_fit();
	}
#endif

	inputResistivityValues(nElem, inFile);

	inFile.close();

}

// Get resisitivity block ID from element ID
int ResistivityBlock::getBlockIDFromElemID(const int ielem) const {

	return m_elementID2blockID[ielem];

}


// Get model ID from block ID
int ResistivityBlock::getModelIDFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_blockID2modelID[ iblk ];
}

// Get block ID from model ID
int ResistivityBlock::getBlockIDFromModelID( const int imdl ) const{

	assert( imdl >= 0 );
	assert( imdl < getNumberOfUnfixedResistivityParameters());

	return m_modelID2blockID[ imdl ];

}

// Get total number of resistivity blocks
int ResistivityBlock::getNumResistivityBlockTotal() const{
	return m_numResistivityBlockTotal;
}

// Get arrays of elements belonging to each resistivity block
const std::vector<int>&  ResistivityBlock::getBlockID2Elements( const int iBlk ) const{

	return m_blockID2Elements[iBlk];

}
