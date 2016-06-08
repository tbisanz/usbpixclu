#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <memory>

static const size_t FE_ROWS = 336;
static const size_t FE_COLS = 80;

std::vector<std::string> split(const std::string& str, const std::string& delim)
{
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

/** The totDecoder class holds the calibration parameters for the 
 *  ToT to charge calibration. It is used as a (bad) static variable
 *  for the hit class. */
class totDecoder {
  public:
	/** Get the charge from ToT */
	double getQ(int x, int y, int tot) {
		return parA.at(x-1).at(y-1) + parB.at(x-1).at(y-1)*tot + parC.at(x-1).at(y-1)*tot*tot;
	}

	/** Default constructor, needed because the hit class has a static totDecoder.
 	  * WARNING: THIS IS A MAJOR DESIGN FLAW!!! DO NOT DO THIS */
	totDecoder() = default;

	/** Actual constructor which should be used */
	totDecoder(std::string filename) {

		//The calibration parameters are stored in a ROOT file
		TFile f(filename.c_str());

		TH1F* ParA = (TH1F*)f.Get("ParA_000");
		TH1F* ParB = (TH1F*)f.Get("ParB_00");
		TH1F* ParC = (TH1F*)f.Get("ParC_00");
	
		//The 2D-histogram values are cached for fast access	
		for(size_t i = 0; i < 80; ++i){
			for(size_t j = 0; j < 336; ++j){
				parA.at(i).at(j) = ParA->GetBinContent(i+1, j+1);
				parB.at(i).at(j) = ParB->GetBinContent(i+1, j+1);
				parC.at(i).at(j) = ParC->GetBinContent(i+1, j+1);
			}
		}
	}
  private:
	std::array<std::array<double,FE_ROWS>, FE_COLS> parA;
	std::array<std::array<double,FE_ROWS>, FE_COLS> parB;
	std::array<std::array<double,FE_ROWS>, FE_COLS> parC;
};

/** A rawHit is the (POD) object representing a data record (DR) in the
 *  FE-I4B datastream. It does not necessarily correspond to only a hit,
 *  as (depending on HitDiscConfig) a DR can also indicate a delayed
 *  hit or contain two hits for neighbouring pixels as well. The ToT 
 *  values correspond to ToT code, not real ToT. */
struct rawHit{
	int x;
	int y;
	int tot1;
	int tot2;
	int lvl1;
	rawHit(int x, int y, int tot1, int tot2, int lvl1):
	x(x), y(y), tot1(tot1), tot2(tot2), lvl1(lvl1){}
};

/** A hit is the entity describing an actually hit pixel. The
 *  ToT it contains, corresponds to the real ToT (not ToT code). 
 *  Using a static instance of a totDecoder, the real ToT is also
 *  transformed into charge (in electrons) */
class hit{
  private:
	static totDecoder decoder;
  public:	
	static void setDecoder(std::string path);// 
	int x;
	int y;
	int tot;
	bool smallTot;
	double charge;
	int lvl1;
	hit(int x, int y, int tot, int lvl1, bool isSmall = false): 
	x(x), y(y), tot(tot), smallTot(isSmall), lvl1(lvl1), charge(0.) {
		charge = decoder.getQ(x, y, tot);
	}
};

//Static member(s/ functions)  of the hit class
totDecoder hit::decoder = totDecoder();
void hit::setDecoder(std::string path) { std::cout << "Setting decoder..." << std::endl; decoder = totDecoder(path);} 

//Helper function for cantor encoding pixel indices
inline int cantorPair(int x, int y){
	return y+(x+y)*(x+y+1)/2.0;
}


std::vector<std::vector<hit>> clusterHits( std::vector<hit> const & hits ){
	if(hits.size() == 1){
		std::vector<std::vector<hit>> result;
		result.push_back(hits);
		return result;
	} else {
		std::vector<hit> hitPixelVec = hits;
		std::vector<hit> newlyAdded;		
		std::vector<std::vector<hit>> clusters;

		while( !hitPixelVec.empty() ) {
			clusters.push_back( std::vector<hit>() );			
			newlyAdded.push_back( hitPixelVec.front() );
			clusters.back().push_back( hitPixelVec.front() );
			hitPixelVec.erase( hitPixelVec.begin() );
			
			while( !newlyAdded.empty() ) {
				bool newlyDone = true;
				int  x1, x2, y1, y2, lv1, lv2, dX, dY, dLv;

				for( auto candidate = hitPixelVec.begin(); candidate != hitPixelVec.end(); ++candidate ){

					//get the relevant infos from the newly added pixel
					x1 = newlyAdded.front().x;
					y1 = newlyAdded.front().y;
					lv1 = newlyAdded.front().lvl1;

					//and the pixel we test against
					x2 = candidate->x;
					y2 = candidate->y;
					lv2 = candidate->lvl1;

					dX = x1-x2;
					dY = y1-y2;
					dLv = lv1-lv2;

					int spatDist = dX*dX+dY*dY;
					int tempDist = dLv*dLv;
	
					if( spatDist <= 8 && tempDist <= 9) {
						newlyAdded.push_back( *candidate );	
						clusters.back().push_back( *candidate );
						hitPixelVec.erase( candidate );
						newlyDone = false;
						break;
					}
				}
				if(newlyDone) newlyAdded.erase( newlyAdded.begin() );
			}	
		}	
	return std::move(clusters);
	}			
}

/** Decoder for decoding rawHits into hits with HitDiscConfig=0,
 *  The ToT Codes correspond to following values:
 *  15(E): No Hit
 *  14(D): Delayed Hit
 *  0-13: ToT-1 (i.e. ToT Code 7 corresponds to a hit with real ToT 8) */
std::vector<hit> decodeHitsHitDisc0( std::vector<rawHit>& rawHits ) {
	std::vector<hit> result;
	std::map<int, int> delayedHit;

	//There are several cases we have to consider:
	//ToT2 = 15 - no hit in adjacent pixel
	//ToT2 = 14 - Corresponds to a delayed hit, the delayed hit has the LvL1 of
	//the DR which reported the delayed hit
	for(auto& raw: rawHits){
		if(raw.tot2 == 15){
			auto delayedHitIt = delayedHit.find(cantorPair(raw.x, raw.y));
			int lvl1 = 0;

			if(delayedHitIt != delayedHit.end()){
					lvl1 = delayedHitIt->second;
					delayedHit.erase(delayedHitIt);
			} else {
					lvl1 = raw.lvl1;
			}
			result.emplace_back(raw.x, raw.y, raw.tot1+1, lvl1);

		} else if(raw.tot2 == 14) {
			delayedHit[cantorPair(raw.x, raw.y+1)] = raw.lvl1;
		} else {
			result.emplace_back(raw.x, raw.y, raw.tot1+1, raw.lvl1);
			result.emplace_back(raw.x, raw.y+1, raw.tot2+1, raw.lvl1);
		}
	}
	return std::move(result);
}

std::vector<hit> decodeHitsHitDisc1( std::vector<rawHit>& rawHits ) {
	std::vector<hit> result;

	for(auto& raw: rawHits){
		if(raw.tot2 == 15){
			int trueTot = 0;
			if(raw.tot1 == 14) trueTot = 1;
			else trueTot = raw.tot1+2;
			result.emplace_back(raw.x, raw.y, trueTot, raw.lvl1);
		} else {
			int trueTot1 = 0;
			int trueTot2 = 0;

			if(raw.tot1 == 14) trueTot1 = 1;
			else trueTot1 = raw.tot1+2;
			if(raw.tot2 == 14) trueTot2 = 1;
			else trueTot2 = raw.tot2+2;

			result.emplace_back(raw.x, raw.y, trueTot1, raw.lvl1);
			result.emplace_back(raw.x, raw.y+1, trueTot2, raw.lvl1);
		}
	}
	return std::move(result);
}


std::vector<hit> decodeHitsHitDisc2( std::vector<rawHit>& rawHits ) {
	std::vector<hit> result;
	
	bool dumpHit = false; //

	for(auto& raw: rawHits){
		if(raw.tot2 == 15){
			int trueTot = 0;
			bool isSmall = false;
			if(raw.tot1 == 14) {
				trueTot = 1;
				dumpHit = true; //
				//std::cout << "ToT 1 hit" << std::endl; // 
				isSmall= true;
			} else {
				trueTot = raw.tot1+3;
			}
			result.emplace_back(raw.x, raw.y, trueTot, raw.lvl1, isSmall);
		} else {
			int trueTot1 = 0;
			int trueTot2 = 0;
			bool isSmall1 = false;
			bool isSmall2 = false;

			if(raw.tot1 == 14) {
				trueTot1 = 1;
				isSmall1 = true;
			} else {
				trueTot1 = raw.tot1+3;
			}		
			if(raw.tot2 == 14) {
				trueTot2 = 1;
				isSmall2 = true;
			} else {
				trueTot2 = raw.tot2+3;
			}
			result.emplace_back(raw.x, raw.y, trueTot1, raw.lvl1, isSmall1);
			result.emplace_back(raw.x, raw.y+1, trueTot2, raw.lvl1, isSmall2);
		}
		//if(dumpHit) std::cout << "Hit!" << std::endl; //
		//if(dumpHit && rawHits.size() == 1) std::cout << "It happened!" << std::endl; //
	}
	return std::move(result);
}

int main(int argc, char* argv[]) {

	const std::string DH = "DH";
	const std::string DR = "DR";
	const std::string TD = "TD";
	const std::string CHANNEL = "CHANNEL";

	if( argc != 3)  {
		std::cout << "Usage: ./main data_file par_file" << std::endl; 
		return 0;
	}
	std::string inFilename = std::string(argv[1]);
	std::string parFilename = std::string(argv[2]);
	
	std::string suffix = split(inFilename, "/").back();
	suffix.resize(suffix.size()-4);

	std::cout << "Output suffix: " << suffix << std::endl;

	std::string outFileName = "out_"+suffix+".root";
	TFile outFile(outFileName.c_str(), "RECREATE");
	outFile.cd();
	//outFile.SetDirectory(gDirectory);

	TH1D* clusterToTHist = new TH1D("totalToT", "Total ToT for all cluster sizes;true ToT;entries", 30, -0.5, 29.5);
	TH1D* clusterQHist = new TH1D("totalQ", "Total charge for all cluster sizes;charge/e;entries", 40, 0, 50000);
	TH1D* clusterToTHist1 = new TH1D("totalToTsize1", "Total ToT for cluster with size 1;true ToT;entries", 30, -0.5, 29.5);
	TH1D* clusterToTHist2 = new TH1D("totalToTsize2", "Total ToT for clusters with size 2;true ToT;entries", 35, -0.5, 34.5);
	TH1D* clusterToTHist3 = new TH1D("totalToTsize3", "Total ToT for clusters with size 3;true ToT;entries", 60, -0.5, 59.5);
	TH1D* clusterToTHist4 = new TH1D("totalToTsize4", "Total ToT for clusters with size 4 or larger;true ToT;entries", 60, -0.5, 59.5);
	
	TH1D* totHit = new TH1D("totHit", "ToT distribution for all hits (w/o clustering);true ToT;entries", 18, -0.5, 17.5);
	TH1D* lvl1Hit = new TH1D("lvl1", "LvL1 distribution for all hits (w/o clustering);lvl1;entries", 16, -0.5, 15.5);
	
	TH1D* cluSize = new TH1D("cluSize", "Cluster size;Size/pixels;entries", 20, -0.5, 19.5);
	TH1D* noClu = new TH1D("noClu", "Number of clusters per event (=per trigger or read-out block);No of clusters;entries", 10, -0.5, 9.5);

	int _pHitDiscConf = 0;
	size_t _pLv1ReadOut = 16;
/*
	std::cout << "FE-I4B Clustering Script!\nPlease configure your setup, \ndefault values are indicated in [].\nSimply pressing return confirms those\nif you don't want to change them!" << std::endl;
	std::cout << "HitDiscConf [" << _pHitDiscConf << "]: ";
		
	char s[10];
	std::cin.get( s, 10);

	//We have to do some ugly ASCII parsing, numbers start at 48 (for 0) and
	//10 corresponds to a NEWLINE character
	switch(atoi(s)){
		case 48:
			_pHitDiscConf = 0;
			break;
		case 49:
			_pHitDiscConf = 1;
			break;
		case 50:
			_pHitDiscConf = 2;
			break;
		case 10:
			break;
		default:
			return 1;
	}

	std::cout << "HitDiscConf set to: " << _pHitDiscConf << std::endl;
*/	
	hit::setDecoder(parFilename);
 	
	std::string cmd; 
	std::string cmdDec; 
	

	std::ifstream infile(inFilename);

	size_t DHCount = 0;
	size_t TRCount = 0;
	size_t trigger = 0;
	
	std::vector<rawHit> hitVec;

	//We always read two lines, one correspods to the decoded word (cmdDec) and one to the
	//raw undecoded hit, which we usually just discard
	while (std::getline(infile, cmdDec)) { //&& std::getline(infile, cmdDec)){
			
		//Sometimes there is a "CHANNEL X" in the data file, this confuses our parser and
		//we need to correct for this, in this case the cmd is actually the second line we
		//read and the actual cmdDec is the third. The first line - i.e. "CHANNEL X" can
		//safely be discarded
		
		if( std::equal( cmdDec.cbegin(), cmdDec.cbegin()+7, CHANNEL.begin() ) ) {
			//cmd = cmdDec;
			//std::getline(infile, cmdDec);	
			if( DHCount =! 0) {
				DHCount = _pLv1ReadOut;
			}
		}
		
		//The LvL1 is determined by counting the data header (DH), indicated by an entry
		//starting with "DH". Records containing data are so called data records, they start
		//with a "DR". In case an external trigger was provided, we also have trigger words
		//("TD") in the datastream. They are not used but can be counted for statistics and 
		//verification
		if( std::equal( cmdDec.cbegin(), cmdDec.cbegin()+2, DH.begin() ) ) {
				DHCount++;
		} else if( std::equal( cmdDec.cbegin(), cmdDec.cbegin()+2, DR.begin() ) ) {
			int x, y, tot1, tot2;
			std::string word;
			std::stringstream line(cmdDec);
			if(line >> word >> x >> y >> tot1 >> tot2) {
				hitVec.emplace_back(x, y, tot1, tot2, DHCount-1);
			} else {
				//PANIC!
			}
		} else if( std::equal( cmdDec.cbegin(), cmdDec.cbegin()+2, TD.begin() ) ) {
				TRCount++;
				if(DHCount != 0) {
					DHCount = _pLv1ReadOut;
				}
		}

		//We always read out a fixed number of data headers - if we reach this number, we will
		//start clustering and increment the (internal) trigger counter 
		if(DHCount == _pLv1ReadOut) {
			DHCount = 0;
			trigger++;
			if(trigger%250000 == 0) std::cout << "Processed " << trigger << " triggers" << std::endl;

			std::vector<hit> decodedHits;
			switch(_pHitDiscConf) {
				case 0:
					decodedHits = decodeHitsHitDisc0(hitVec);
					break;
				case 1:
					decodedHits = decodeHitsHitDisc1(hitVec);
					break;
				case 2:
					decodedHits = decodeHitsHitDisc2(hitVec);
					break;			
			}
				
			//	if(!decodedHits.empty()) std::cout << "New event block: " << std::endl;
			auto clusters = clusterHits(decodedHits);
			/*		
				if(!decodedHits.empty()) std::cout << "All hits: " << std::endl;
				for(auto& pixel: decodedHits){
					std::cout << pixel.x << "|" << pixel.y << "|" << pixel.tot << std::endl;
				}
			*/
			noClu->Fill(clusters.size());	
			for(auto& cluster: clusters){
				double qtot = 0;
				double tot = 0;
				for(auto& pixel: cluster) {
					qtot += pixel.charge;
					tot += pixel.tot;
					totHit->Fill(pixel.tot);
					lvl1Hit->Fill(pixel.lvl1);
				}
				clusterToTHist->Fill(tot);
				clusterQHist->Fill(qtot);
				cluSize->Fill(cluster.size());
				if(cluster.size() == 1) clusterToTHist1->Fill(tot);
				else if(cluster.size() == 2) clusterToTHist2->Fill(tot);
				else if(cluster.size() == 3) clusterToTHist3->Fill(tot);
				else clusterToTHist4->Fill(tot);
/*
				if(cluster.size() == 2) {
					std::cout << "Two hit cluster: " << std::endl;
					for(auto& pixel: cluster) { std::cout << "Pixel: " << pixel.x << "|" << pixel.y << "|" << pixel.tot << "|" << pixel.lvl1 << std::endl;}
				}*/
			}
			hitVec.clear();
		}
	}
	//outFile->cd();
	auto c1 = std::unique_ptr<TCanvas>(new TCanvas());

	clusterQHist->Draw();	
	std::string clusterQName = "clusterQ_"+suffix+".pdf";
	c1->SaveAs(clusterQName.c_str(),"pdf");

	clusterToTHist->Draw();	
	std::string clusterToTName = "clusterToT_"+suffix+".pdf";
	c1->SaveAs(clusterToTName.c_str(),"pdf");

	clusterToTHist1->Draw();	
	std::string clusterToTSize1 = "clusterToT_size1_"+suffix+".pdf";
	c1->SaveAs(clusterToTSize1.c_str(),"pdf");

	clusterToTHist2->Draw();	
	std::string clusterToTSize2 = "clusterToT_size2_"+suffix+".pdf";
	c1->SaveAs(clusterToTSize2.c_str(),"pdf");

	clusterToTHist3->Draw();	
	std::string clusterToTSize3 = "clusterToT_size3_"+suffix+".pdf";
	c1->SaveAs(clusterToTSize3.c_str(),"pdf");

	clusterToTHist4->Draw();	
	std::string clusterToTSize4 = "clusterToT_size4orlarger_"+suffix+".pdf";
	c1->SaveAs(clusterToTSize4.c_str(),"pdf");

	totHit->Draw();	
	std::string totAllHit = "totAllHit_"+suffix+".pdf";
	c1->SaveAs(totAllHit.c_str(),"pdf");

	lvl1Hit->Draw();	
	std::string lvl1AllHit = "lvl1AllHit_"+suffix+".pdf";
	c1->SaveAs(lvl1AllHit.c_str(),"pdf");

	c1->SetLogy();
	cluSize->Draw();
	std::string clusterSize = "clusterSize_"+suffix+".pdf";
	c1->SaveAs(clusterSize.c_str(),"pdf");

	noClu->Draw();	
	std::string noClusters = "numberClusters_"+suffix+".pdf";
	c1->SaveAs(noClusters.c_str(),"pdf");

	outFile.Write();
	outFile.Close();

	std::cout << "Processed " << trigger << " triggers" << std::endl;
	std::cout << "Processed " << TRCount << " TDs" << std::endl;
	return 1;
}
