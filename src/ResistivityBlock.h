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
#ifndef DBLDEF_RESISTIVITY_BLOCK
#define DBLDEF_RESISTIVITY_BLOCK

#include <vector>
#include <iostream>
#include <cstdlib>

// Class of resistivity blocks
class ResistivityBlock{

public:

	// Constructer
	ResistivityBlock();

	// Destructer
	virtual ~ResistivityBlock();

	// Get number of unfixed resistivity parameters
	virtual int getNumberOfUnfixedResistivityParameters() const = 0;

	// Read data of resisitivity block model from input file
	void inputResistivityBlock(const int iterationNumInit);

	// Get resisitivity block ID from element ID
	int getBlockIDFromElemID(const int ielem) const;

	// Get arrays of elements belonging to a specified resistivity block
	const std::vector<int>& getBlockID2Elements(const int iBlk) const;

	// Get model ID from block ID
	int getModelIDFromBlockID( const int iblk ) const;

	// Get block ID from model ID
	int getBlockIDFromModelID( const int imdl ) const;

	// Get total number of resistivity blocks
	int getNumResistivityBlockTotal() const;
	
protected:

	// Total number of resistivity blocks
	int m_numResistivityBlockTotal;

	// Array of the resistivity block IDs of each element
	int* m_elementID2blockID;

	// Arrays of elements belonging to each resistivity model
	std::vector<int>* m_blockID2Elements;

	// Array convert model IDs to IDs of resistivity block  
	int* m_modelID2blockID;

	// Array convert IDs of resistivity block to model IDs 
	int* m_blockID2modelID;

	// Array of minimum resistivity values of each block
	double* m_resistivityValuesMin;

	// Array of maximum resistivity values of each block
	double* m_resistivityValuesMax;

private:

	// Copy constructer
	ResistivityBlock(const ResistivityBlock& rhs){
		std::cerr << "Error : Copy constructer of the class ResistivityBlock is not implemented." << std::endl;
		exit(1);
	};

	// Assignment operator
	ResistivityBlock& operator=(const ResistivityBlock& rhs){
		std::cerr << "Error : Assignment operator of the class ResistivityBlock is not implemented." << std::endl;
		exit(1);
	};

	// Read reslstivity values from input file
	virtual void inputResistivityValues( const int nElem, std::ifstream& inFile ) = 0;

};

#endif
