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

		TH1F* ParA = (TH1F*)f.Get("ParA_00");
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
	
					if( spatDist <= 2 && tempDist <= 9) {
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

int main() {

	const std::string DH = "DH";
	const std::string DR = "DR";
	const std::string TD = "TD";
	const std::string CHANNEL = "CHANNEL";

	std::string outFileName = "out.root";
	TFile outFile(outFileName.c_str(), "RECREATE");
	outFile.cd();
	//outFile.SetDirectory(gDirectory);

	TH1D* totClu = new TH1D("test", "test", 20, -0.5, 19.5);
	TH1D* totClu1 = new TH1D("test1", "test1", 20, -0.5, 19.5);
	TH1D* totClu2 = new TH1D("test2", "test2", 25, -0.5, 24.5);
	TH1D* totClu3 = new TH1D("test3", "test3", 60, -0.5, 59.5);
	TH1D* totClu4 = new TH1D("test4", "test4", 60, -0.5, 59.5);
	
	TH1D* totHit = new TH1D("totHit", "totHit", 18, -0.5, 17.5);
	TH1D* lvl1Hit = new TH1D("lvl1", "lvl1", 16, -0.5, 15.5);
	
	TH1D* cluSize = new TH1D("cluSize", "cluSize", 20, -0.5, 19.5);
	TH1D* noClu = new TH1D("noClu", "noClu", 10, -0.5, 9.5);

	int _pHitDiscConf = 2;
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
	hit::setDecoder("par.root");
 	
	std::string cmd; 
	std::string cmdDec; 

	//std::ifstream infile("quellenscan.21.04_SOURCE_SCAN_5_0_0_0.raw");
	std::ifstream infile("messung1_SOURCE_SCAN_17_0_0_0.raw");

	size_t DHCount = 0;
	size_t TRCount = 0;
	size_t trigger = 0;
	
	std::vector<rawHit> hitVec;

	//We always read two lines, one correspods to the decoded word (cmdDec) and one to the
	//raw undecoded hit, which we usually just discard
	while (std::getline(infile, cmd) && std::getline(infile, cmdDec)){
			
		//Sometimes there is a "CHANNEL X" in the data file, this confuses our parser and
		//we need to correct for this, in this case the cmd is actually the second line we
		//read and the actual cmdDec is the third. The first line - i.e. "CHANNEL X" can
		//safely be discarded
		if( std::equal( cmd.cbegin(), cmd.cbegin()+7, CHANNEL.begin() ) ) {
			cmd = cmdDec;
			std::getline(infile, cmdDec);	
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
				hitVec.emplace_back(x, y, tot1, tot2, DHCount);
			} else {
				//PANIC!
			}
		} else if( std::equal( cmdDec.cbegin(), cmdDec.cbegin()+2, TD.begin() ) ) {
				TRCount++;
		}

		//We always read out a fixed number of data headers - if we reach this number, we will
		//start clustering and increment the (internal) trigger counter 
		if(DHCount == _pLv1ReadOut) {
			DHCount = 0;
			trigger++;

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
				totClu->Fill(tot);
				cluSize->Fill(cluster.size());
				if(cluster.size() == 1) totClu1->Fill(tot);
				else if(cluster.size() == 2) totClu2->Fill(tot);
				else if(cluster.size() == 3) totClu3->Fill(tot);
				else totClu4->Fill(tot);
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
	totClu->Draw();	
	c1->SaveAs("test.pdf","pdf");
	totClu1->Draw();	
	c1->SaveAs("test1.pdf","pdf");
	totClu2->Draw();	
	c1->SaveAs("test2.pdf","pdf");
	totClu3->Draw();	
	c1->SaveAs("test3.pdf","pdf");
	totClu4->Draw();	
	c1->SaveAs("test4.pdf","pdf");
	totHit->Draw();	
	c1->SaveAs("totHit.pdf","pdf");
	lvl1Hit->Draw();	
	c1->SaveAs("lvl1Hit.pdf","pdf");
	c1->SetLogy();
	cluSize->Draw();	
	c1->SaveAs("cluSize.pdf","pdf");
	noClu->Draw();	
	c1->SaveAs("noClu.pdf","pdf");

	outFile.Write();
	outFile.Close();

	std::cout << "Processed " << trigger << " triggers" << std::endl;
	std::cout << "Processed " << TRCount << " triggers" << std::endl;
	return 1;
}