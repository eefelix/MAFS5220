#include "DataImport.h"

//EquityForwardData InstrumentData::LoadEquityForward(string tradeRef) {
//
//	EquityForwardData oEquityForwardData;
//
//	ifstream filein(tradeRef + ".csv");
//
//	for (string line; getline(filein, line);){
//		stringstream lineStream(line);
//		string cell;
//		int cell_idx = 0;
//
//		string item_;
//
//		while (getline(lineStream, cell, ',')){
//			switch (cell_idx) {
//			case 0:		// Field Name
//				item_ = cell; break;
//			case 1:		// Field Value
//				if (item_ == "Code") {
//					oEquityForwardData.Code_ = cell;
//				}
//				else if (item_ == "Underlying") {
//					oEquityForwardData.Underlying_ = cell;
//				}
//				else if (item_ == "Strike") {
//					oEquityForwardData.Strike_ = atof(cell.c_str());
//				}
//				else if (item_ == "SettlementDay") {
//					oEquityForwardData.SettlementDay_ = atof(cell.c_str());
//				}
//				else if (item_ == "PositionType") {
//					oEquityForwardData.PositionType_ = cell;
//				}
//				else if (item_ == "ValueDay") {
//					oEquityForwardData.ValueDay_ = atof(cell.c_str());
//				}
//				else if (item_ == "MaturityYear") {
//					oEquityForwardData.MaturityYear_ = atof(cell.c_str());
//				}
//				break;
//			}
//			cell_idx++;
//		}
//	}
//	return oEquityForwardData;
//}


std::string MarketData::GetCurrentWorkingDirectory()
{
	char temp_working_directory[MAX_PATH + 1];
	GetCurrentDirectoryA(sizeof(temp_working_directory), temp_working_directory); // **** win32 specific ****
	std::string workingDirectory = temp_working_directory;
	if (workingDirectory.substr(workingDirectory.size() - 1) != "\\") {
		workingDirectory.append("\\");
	}
	return workingDirectory;
}

bool MarketData::IsFileExist(string filePath) {
	if (std::ifstream(filePath)) {
		return true;
	}
	else {
		return false;
	}
}

EquityMarketData MarketData::LoadEquity(string bbgCode) {

	EquityMarketData oEquityMarketData;
	
	if (!IsFileExist(GetCurrentWorkingDirectory() + bbgCode + ".csv")) {
		cout << GetCurrentWorkingDirectory() + bbgCode << ".csv not exist. Program terminated." << endl;
		system("pause");
		std::exit;
	}

	ifstream filein(GetCurrentWorkingDirectory() + bbgCode + ".csv");

	for (string line; getline(filein, line);){
		stringstream lineStream(line);
		string cell;
		int cell_idx = 0;

		string item_;

		while (getline(lineStream, cell, ',')){
			switch (cell_idx) {
				case 0:		// Field Name
					item_ = cell; break;
				case 1:		// Field Value
					if (item_ == "Code") {
						oEquityMarketData.Code_ = cell;
					}
					else if (item_ == "Spot") {
						oEquityMarketData.Spot_ = atof(cell.c_str());
					}
					else if (item_ == "Dividend") {
						oEquityMarketData.Dividend_ = atof(cell.c_str()); 
					}
					else if (item_ == "Drift") {
						oEquityMarketData.Drift_ = atof(cell.c_str());
					}
					else if (item_ == "Volatility") {
						oEquityMarketData.Volatility_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds6Mths") {
						oEquityMarketData.Cds6Mths_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds1Yrs") {
						oEquityMarketData.Cds1Yrs_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds2Yrs") {
						oEquityMarketData.Cds2Yrs_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds3Yrs") {
						oEquityMarketData.Cds3Yrs_ = atof(cell.c_str());
					}
					else if (item_ == "Cds4Yrs") {
						oEquityMarketData.Cds4Yrs_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds5Yrs") {
						oEquityMarketData.Cds5Yrs_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds7Yrs") {
						oEquityMarketData.Cds7Yrs_ = atof(cell.c_str()); 
					}
					else if (item_ == "Cds10Yrs") {
						oEquityMarketData.Cds10Yrs_ = atof(cell.c_str());
					}
					break;
				}
				cell_idx++;
			}
		}
	return oEquityMarketData;
}