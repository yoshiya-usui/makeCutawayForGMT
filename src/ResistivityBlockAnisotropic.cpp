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
#include "ResistivityBlockAnisotropic.h"
#include "Util.h"

#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>
#include <assert.h>


// Constructer
ResistivityBlockAnisotropic::ResistivityBlockAnisotropic():
	m_numAnisotropicResistivityParametersNotFixed(0)
{}

// Destructer
ResistivityBlockAnisotropic::~ResistivityBlockAnisotropic(){
}

// Read reslstivity values from input file
void ResistivityBlockAnisotropic::inputResistivityValues(const int nElem, std::ifstream& inFile) {

	m_anisotropicResistivityParameters.reserve(m_numResistivityBlockTotal);
	m_anisotropicResistivityParametersPre.reserve(m_numResistivityBlockTotal);
	m_anisotropicResistivityParametersUpdatedFull.reserve(m_numResistivityBlockTotal);
	for (int iParam = 0; iParam < 6; ++iParam) {
		m_anisotropicResistivityParameter2modelID[iParam].reserve(m_numResistivityBlockTotal);
		m_fixedAnisotropicResistivityParameters[iParam].reserve(m_numResistivityBlockTotal);
	}

	if (m_resistivityValuesMin != NULL) {
		delete[] m_resistivityValuesMin;
	}
	m_resistivityValuesMin = new double[m_numResistivityBlockTotal];
	if (m_resistivityValuesMax != NULL) {
		delete[] m_resistivityValuesMax;
	}
	m_resistivityValuesMax = new double[m_numResistivityBlockTotal];
	for (int i = 0; i < m_numResistivityBlockTotal; ++i) {
		m_resistivityValuesMin[i] = 0.0;
		m_resistivityValuesMax[i] = 0.0;
	}

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		int idum(0);
		inFile >> idum;
		if (idum != iBlk) {
			std::cerr << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		inFile >> idum;
		const int typeOfAnistropy = idum;
		m_typeOfAnisotropy.push_back(typeOfAnistropy);
		AnistropicResistivityParameters anistropicResistivityParam = { -1.0, -1.0, -1.0, 0.0, 0.0, 0.0 };
		int isFixedArray[6] = { -1, -1, -1, -1, -1, -1 };
		if (typeOfAnistropy == ResistivityBlockAnisotropic::ISOTROPY) {
			inFile >> anistropicResistivityParam.rhoXX;
			if (anistropicResistivityParam.rhoXX <= 0.0) {
				std::cerr << "Error : The resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoXX << std::endl;
				exit(1);
			}
			anistropicResistivityParam.rhoYY = anistropicResistivityParam.rhoXX;
			anistropicResistivityParam.rhoZZ = anistropicResistivityParam.rhoXX;
			inFile >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk];
			if (m_resistivityValuesMin[iBlk] <= 0.0) {
				std::cerr << "Error : Minimum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMin[iBlk] << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] <= 0.0) {
				std::cerr << "Error : Maximum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMax[iBlk] << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoXX) {
				std::cerr << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoXX) {
				std::cerr << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			int isFixed(0);
			inFile >> isFixed;
			isFixedArray[RHO_XX] = isFixed;
			isFixedArray[RHO_YY] = 1;
			isFixedArray[RHO_ZZ] = 1;
			isFixedArray[STRIKE] = 1;
			isFixedArray[DIP] = 1;
			isFixedArray[SLANT] = 1;
		}
		else if (typeOfAnistropy == TRANSVERSE_ISOTROPY || typeOfAnistropy == GENERAL_ANISOTROPY) {
			inFile >> anistropicResistivityParam.rhoXX >> anistropicResistivityParam.rhoYY;
			if (anistropicResistivityParam.rhoXX <= 0.0) {
				std::cerr << "Error : The xx-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoXX << std::endl;
				exit(1);
			}
			if (anistropicResistivityParam.rhoYY <= 0.0) {
				std::cerr << "Error : The yy-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoYY << std::endl;
				exit(1);
			}
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				anistropicResistivityParam.rhoZZ = anistropicResistivityParam.rhoXX;
			}
			else{
				inFile >> anistropicResistivityParam.rhoZZ;
				if (anistropicResistivityParam.rhoZZ <= 0.0) {
					std::cerr << "Error : The zz-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoZZ << std::endl;
					exit(1);
				}
			}
			inFile >> anistropicResistivityParam.strike >> anistropicResistivityParam.dip;
			anistropicResistivityParam.strike *= CommonParameters::deg2rad;
			anistropicResistivityParam.dip *= CommonParameters::deg2rad;
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				anistropicResistivityParam.slant = 0.0;
			}
			else{
				inFile >> anistropicResistivityParam.slant;
				anistropicResistivityParam.slant *= CommonParameters::deg2rad;
			}
			inFile >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk];
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoXX) {
				std::cerr << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the xx-component ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoXX) {
				std::cerr << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the xx-component ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoYY) {
				std::cerr << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the yy-component ( " << anistropicResistivityParam.rhoYY << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoYY) {
				std::cerr << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the yy-component ( " << anistropicResistivityParam.rhoYY << " )." << std::endl;
				exit(1);
			}
			if (typeOfAnistropy == GENERAL_ANISOTROPY) {
				if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoZZ) {
					std::cerr << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the zz-component ( " << anistropicResistivityParam.rhoZZ << " )." << std::endl;
					exit(1);
				}
				if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoZZ) {
					std::cerr << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the zz-component ( " << anistropicResistivityParam.rhoZZ << " )." << std::endl;
					exit(1);
				}
			}
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				inFile >> isFixedArray[RHO_XX] >> isFixedArray[RHO_YY];
				isFixedArray[RHO_ZZ] = 1;
				inFile >> isFixedArray[STRIKE] >> isFixedArray[DIP];
				isFixedArray[SLANT] = 1;
			}
			else {
				inFile >> isFixedArray[RHO_XX] >> isFixedArray[RHO_YY] >> isFixedArray[RHO_ZZ];
				inFile >> isFixedArray[STRIKE] >> isFixedArray[DIP] >> isFixedArray[SLANT];
			}
		}
		else {
			std::cerr << "Error : Unsupported type of anisotropy (" << typeOfAnistropy << ")" << std::endl;
			exit(1);
		}
		m_anisotropicResistivityParameters.push_back(anistropicResistivityParam);
		m_anisotropicResistivityParametersPre.push_back(anistropicResistivityParam);
		m_anisotropicResistivityParametersUpdatedFull.push_back(anistropicResistivityParam);
		for (int iParam = 0; iParam < 6; ++iParam) {
			m_fixedAnisotropicResistivityParameters[iParam].push_back(isFixedArray[iParam]);
			if (isFixedArray[iParam]) {
				// Fixed
				m_anisotropicResistivityParameter2modelID[iParam].push_back(-1);
			}
			else {
				// Changable
				m_anisotropicResistivityParameter2modelID[iParam].push_back(m_numAnisotropicResistivityParametersNotFixed);
				++m_numAnisotropicResistivityParametersNotFixed;
			}
		}
	}

	inFile.close();

	if (m_fixedAnisotropicResistivityParameters[0][0] < 1) {
		std::cerr << "Error : Resistivity block 0 must be the air, and its resistivity must be fixed." << std::endl;
		exit(1);
	}
	if (m_numAnisotropicResistivityParametersNotFixed <= 0) {
		std::cerr << "Error : Total number of modifiable anistropic conductivity parameters is zero or negative !! : " << m_numAnisotropicResistivityParametersNotFixed << std::endl;
		exit(1);
	}
	m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.reserve(m_numAnisotropicResistivityParametersNotFixed);
	int icount(0);
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		for (int iParam = 0; iParam < 6; ++iParam) {
			if (!m_fixedAnisotropicResistivityParameters[iParam][iBlk]) {
				m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.push_back(std::make_pair(iBlk, iParam));
				++icount;
			}
		}
	}
	assert(icount == m_numAnisotropicResistivityParametersNotFixed);

}


// Get number of anisotropic resistivity parameters whose values are NOT fixed
int ResistivityBlockAnisotropic::getNumberOfUnfixedResistivityParameters() const {

	return m_numAnisotropicResistivityParametersNotFixed;

}

// Get anisotropy type
int ResistivityBlockAnisotropic::getTypeOfAnisotropy(const int iblk) const {

	return m_typeOfAnisotropy[iblk];

}

// Get anisotropic resistivity parameter from resisitivity block ID
double ResistivityBlockAnisotropic::getAnisotropicResistivityParameterFromBlockID(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	const AnistropicResistivityParameters paramValues = getAnisotropicResistivityParametersFromBlockID(iblk);
	switch (iparam)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		return paramValues.rhoXX;
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		return paramValues.rhoYY;
		break;
	case ResistivityBlockAnisotropic::RHO_ZZ:
		return paramValues.rhoZZ;
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		return paramValues.strike;
		break;
	case ResistivityBlockAnisotropic::DIP:
		return paramValues.dip;
		break;
	case ResistivityBlockAnisotropic::SLANT:
		return paramValues.slant;
		break;
	default:
		std::cerr << "Error : Unsupported type of anisotropy (" << iparam << ")" << std::endl;
		exit(1);
		break;
	}

}

// Get anisotropic previous resistivity parameter from resisitivity block ID
double ResistivityBlockAnisotropic::getAnisotropicResistivityParameterPreFromBlockID(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	const AnistropicResistivityParameters paramValues = getAnisotropicResistivityParametersPreFromBlockID(iblk);
	switch (iparam)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		return paramValues.rhoXX;
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		return paramValues.rhoYY;
		break;
	case ResistivityBlockAnisotropic::RHO_ZZ:
		return paramValues.rhoZZ;
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		return paramValues.strike;
		break;
	case ResistivityBlockAnisotropic::DIP:
		return paramValues.dip;
		break;
	case ResistivityBlockAnisotropic::SLANT:
		return paramValues.slant;
		break;
	default:
		std::cerr << "Error : Unsupported type of anisotropy (" << iparam << ")" << std::endl;
		exit(1);
		break;
	}

}

// Get anisotropic resistivity from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParameters[iblk];


}

// Get anisotropic previous resistivity parameters from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersPreFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParametersPre[iblk];

}

// Get anisotropic fully-updated resistivity parameters from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersUpdatedFullFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParametersUpdatedFull[iblk];

}

// Get model ID from block ID and anisotropic parameter
int ResistivityBlockAnisotropic::getModelIDFromBlockIDAndAnisotropicParameter(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	return m_anisotropicResistivityParameter2modelID[iparam][iblk];

}

// Get flag specifing whether anisotrpic conductivity parameters is fixed or not
bool ResistivityBlockAnisotropic::isFixedAnisotropicResistivityParameters(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_fixedAnisotropicResistivityParameters[iparam][iblk];

}

// Get indexes of block and anistropic resistivity parameter from model index
std::pair<int, int> ResistivityBlockAnisotropic::getIndexesOfBlockAndAnisotropicResistivityParameterFromModelIndex(const int iMdl) const {

	return m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes[iMdl];

}

// Confirm that there is no element having general anisotropy
void ResistivityBlockAnisotropic::confirmNoElementHavingGeneralAnisotropy() const {

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::GENERAL_ANISOTROPY) {
			std::cerr << "Error : Resistivity block (" << iBlk << ") have general anisotropy !!" << std::endl;
			exit(1);
		}
	}

}

// Calculate array of fully updated anisotropic resistivity parameters obtained by inversion
void ResistivityBlockAnisotropic::calcUnfixedAnisotropicResistivityParametersUpdatedFull(const double* const updatedModel) {

	int index(0);
	for (std::vector< std::pair<int, int> >::const_iterator itr = m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.begin();
		itr != m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.end(); ++itr, ++index) {
		const int iBlk = itr->first;
		const AnistropicResistivityParameters anisoParamsPre = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
		const int iParam = itr->second;
		switch (iParam)
		{
		case ResistivityBlockAnisotropic::RHO_XX:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoXX = pow(10.0, log10(anisoParamsPre.rhoXX) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::RHO_YY:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoYY = pow(10.0, log10(anisoParamsPre.rhoYY) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::RHO_ZZ:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoZZ = pow(10.0, log10(anisoParamsPre.rhoZZ) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::STRIKE:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].strike = anisoParamsPre.strike + updatedModel[index];
			break;
		case ResistivityBlockAnisotropic::DIP:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].dip = anisoParamsPre.dip + updatedModel[index];
			break;
		case ResistivityBlockAnisotropic::SLANT:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].slant = anisoParamsPre.slant + updatedModel[index];
			break;
		default:
			std::cerr << "Error : Unsupported type of anisotropy (" << iParam << ")" << std::endl;
			exit(1);
			break;
		}
	}

}
