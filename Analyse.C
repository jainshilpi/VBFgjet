#define Analyse_cxx
#include "Analyse.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TRandom3.h"
#include "Utilities.h"
//#ifdef __CLING__
/*R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
//#else
class ExRootTreeReader;
class ExRootResult;
//#endif
*/

constexpr unsigned int str2int(const char* str, int h = 0)
{
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

void Analyse::Loop(string fname, string fnameout)
{
  
  cout<<"Now inside loop, fname : fout : "<<fname<<" "<<fnameout<<endl;
  

  //bool debug = true;
  bool debug = false;
  
  TRandom3 gran;
  double fracTrain = 0.8;

  //string weightOptionFWM = "s";
  string weightOptionFWM = "p";
  //string weightOptionFWM = "T";
  //string weightOptionFWM = "z";
  //string weightOptionFWM = "y";
  //string weightOptionFWM = "1";

  //TTree *tin = list_files(fname.c_str());

  TChain *chain = new TChain("Delphes");
  chain->Add(fname.c_str());


  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  Init(treeReader);


  TFile *fout = new TFile(fnameout.c_str(), "RECREATE");
  map<string,TH1F*> hmap;
  hmap["phoPt"] = new TH1F("phoPt","",200,0,500);
  hmap["jet1Pt"] = new TH1F("jet1Pt","",200,0,500);
  hmap["jet2Pt"] = new TH1F("jet2Pt","",200,0,500);

  hmap["phoEta"] = new TH1F("phoEta","",100,-2.5,2.5);
  hmap["jet1Eta"] = new TH1F("jet1Eta","",200,-5,5);
  hmap["jet2Eta"] = new TH1F("jet2Eta","",200,-5,5);

  hmap["eta1TimesEta2"] = new TH1F("eta1TimesEta2", "", 200, -25, 25);

  hmap["dEta"] = new TH1F("dEta","",200, -10, 10);

  hmap["jetMass"] = new TH1F("jetMass","",500, 0, 3000);
  hmap["jet1Mass"] = new TH1F("jet1Mass","",500, 0, 100);
  hmap["jet2Mass"] = new TH1F("jet2Mass","",500, 0, 100);

  hmap["eta_pho_jetSys"] = new TH1F("eta_pho_jetSys", "", 200, -10,10);
  hmap["eta_pho_jet1"] = new TH1F("eta_pho_jet1", "", 200, -10,10);
  hmap["eta_pho_jet2"] = new TH1F("eta_pho_jet2", "", 200, -10,10);
  
  hmap["eta_pho_jetSys_norm"] = new TH1F("eta_pho_jetSys_norm", "", 200, -10,10);

  hmap["centrality"] = new TH1F("centrality", "", 200, -1,1);

  hmap["phi_pho_jet1"] = new TH1F("phi_pho_jet1", "", 200, -3.2,3.2);
  hmap["phi_pho_jet2"] = new TH1F("phi_pho_jet2", "", 200, -3.2,3.2);

  hmap["phi_pho_jetmin"] = new TH1F("phi_pho_jetmin", "", 200, -3.2,3.2);
  hmap["phi_pho_jetmax"] = new TH1F("phi_pho_jetmax", "", 200, -3.2,3.2);

  hmap["phi_jet1_jet2"] = new TH1F("phi_jet1_jet2", "", 200, -3.2,3.2);
  hmap["phi_pho_jetSys"] = new TH1F("phi_pho_jetSys", "", 200, -3.2,3.2);

  hmap["dr_jet1jet2"] = new TH1F("dr_jet1jet2", "", 200, 0,20);
  hmap["dr_phojet1"] = new TH1F("dr_phojet1", "", 200, 0,20);
  hmap["dr_phojet2"] = new TH1F("dr_phojet2", "", 200, 0,20);

  hmap["dr_phojetmin"] = new TH1F("dr_phojetmin", "", 200, 0,20);
  hmap["dr_phojetmax"] = new TH1F("dr_phojetmax", "", 200, 0,20);

  hmap["pt_phoTojet1"] = new TH1F("pt_phoTojet1", "", 200, 0,20);
  hmap["pt_phoTojet2"] = new TH1F("pt_phoTojet2", "", 200, 0,20);
  hmap["pt_phoTojetSys"] = new TH1F("pt_phoTojetSys", "", 200, 0,20);

  hmap["mass3body"] = new TH1F("mass3body", "", 1000, 0,5000);

  hmap["pt_phobrem_jet"] = new TH1F("pt_phobrem_jet", "", 100, 0,5);

  hmap["sumTrkPt"] = new TH1F("sumTrkPt","",500,0,200);
  hmap["ntracks"] = new TH1F("ntracks","",200,0,200);
  
  hmap["trkPt"] = new TH1F("trkPt","",500,0,50);
  hmap["trkEta"] = new TH1F("trkEta","",200,-3,3);

  hmap["diJetScalarPt"] = new TH1F("diJetScalarPt","",500,0,1500);
  hmap["diJetVectorPt"] = new TH1F("diJetVectorPt","",200,0,800);

  hmap["ratio_trackp_To_Jetp"] = new TH1F("ratio_trackp_To_Jetp","",50,0,1);

  hmap["centrality_trk"] = new TH1F("centrality_trk","",500,-0.1,100);

  ///29th Nov, 2021
  hmap["momEquaVar"] = new TH1F("momEquaVar","",500,-20,20);

  ///9th Jan, 2022
  hmap["fwm0"] = new TH1F("fwm0","",200,-0.1,1.1);
  hmap["fwm1"] = new TH1F("fwm1","",200,-0.1,1.1);
  hmap["fwm2"] = new TH1F("fwm2","",200,-0.1,1.1);
  hmap["fwm3"] = new TH1F("fwm3","",200,-0.1,1.1);
  hmap["fwm4"] = new TH1F("fwm4","",200,-0.1,1.1);
  hmap["fwm5"] = new TH1F("fwm5","",200,-0.1,1.1);
  
  
  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    
    hmap[it->first]->Sumw2();
  }
  
  ////Tree variables
  double phoPt_, phoEta_, jet1Eta_, jet2Eta_, jet1Pt_, jet2Pt_, jetEta1TimesEta2_, jetdEta_, jetMass_, zeppenfeld_, zeppenfeldNorm_, dEta_phoJet1_, dEta_phoJet2_, dPhi_phoJet1_, dPhi_phoJet2_, dPhi_jet1Jet2_, dPhi_phoJetSys_, dR_jet1Jet2_, dR_phoJet1_, dR_phoJet2_, pT_phoTojet1_, pT_phoTojet2_, pT_phoTojetSys_, invMass_3body_, pT_phoBremJet_, dPhi_phoJetMin_, dPhi_phoJetMax_, dR_phoJetMin_, dR_phoJetMax_, sumTrkPt_, centrality_, diJetScalarPt_, diJetVectorPt_, centrality_pt_trk_;
  int nTracks_, isTrain_;
  
  double fwm0_, fwm1_, fwm2_, fwm3_, fwm4_, fwm5_;
  
  TTree *tout = new TTree("vars","Tree of variables");
  tout->Branch("phoPt",&phoPt_);
  tout->Branch("phoEta",&phoEta_);
  tout->Branch("jet1Eta",&jet1Eta_);
  tout->Branch("jet2Eta",&jet2Eta_);
  tout->Branch("jet1Pt",&jet1Pt_);
  tout->Branch("jet2Pt",&jet2Pt_);
  tout->Branch("jetEta1TimesEta2",&jetEta1TimesEta2_);
  tout->Branch("jetdEta",&jetdEta_);
  tout->Branch("jetMass",&jetMass_);
  tout->Branch("zeppenfeld",&zeppenfeld_);
  tout->Branch("zeppenfeldNorm",&zeppenfeldNorm_);
  tout->Branch("centrality",&centrality_);
  tout->Branch("diJetScalarPt",&diJetScalarPt_);
  tout->Branch("diJetVectorPt",&diJetVectorPt_);
  tout->Branch("dEta_phoJet1",&dEta_phoJet1_);
  tout->Branch("dEta_phoJet2",&dEta_phoJet2_);

  tout->Branch("dPhi_phoJet1",&dPhi_phoJet1_);
  tout->Branch("dPhi_phoJet2",&dPhi_phoJet2_);
  tout->Branch("dPhi_jet1Jet2",&dPhi_jet1Jet2_);
  tout->Branch("dPhi_phoJetSys",&dPhi_phoJetSys_);
  tout->Branch("dPhi_phoJetMin",&dPhi_phoJetMin_);
  tout->Branch("dPhi_phoJetMax",&dPhi_phoJetMax_);

  tout->Branch("pT_phoTojet1",&pT_phoTojet1_);
  tout->Branch("pT_phoTojet2",&pT_phoTojet2_);
  tout->Branch("pT_phoTojetSys",&pT_phoTojetSys_);

  tout->Branch("dR_phoJet1",&dR_jet1Jet2_);
  tout->Branch("dR_phoJet2",&dR_phoJet1_);
  tout->Branch("dR_jet1Jet2",&dR_phoJet2_);
  tout->Branch("dR_phoJetMin",&dR_phoJetMin_);
  tout->Branch("dR_phoJetMax",&dR_phoJetMax_);

  tout->Branch("invMass_3body",&invMass_3body_);
  tout->Branch("pT_phoBremJet",&pT_phoBremJet_);

  tout->Branch("sumTrkPt",&sumTrkPt_);
  tout->Branch("nTracks",&nTracks_);
  tout->Branch("isTrain", &isTrain_);  

  tout->Branch("centrality_pt_trk",&centrality_pt_trk_);

  tout->Branch("fwm0",&fwm0_);
  tout->Branch("fwm1",&fwm1_);
  tout->Branch("fwm2",&fwm2_);
  tout->Branch("fwm3",&fwm3_);
  tout->Branch("fwm4",&fwm4_);
  tout->Branch("fwm5",&fwm5_);


  //if (fChain == 0) return;

   Long64_t nentries = treeReader->GetEntries();
   
   std::cout<<"Number of entries are "<<nentries<<std::endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     nb = treeReader->ReadEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int phoInd = selPho();
      int jet1Ind = -99, jet2Ind = -99;
      
      selVBFJets(jet1Ind, jet2Ind);

      

      ////this is the event selection so do all here now
      if(phoInd>=0 && jet1Ind>=0 && jet2Ind>=0){

	//cout<<"=========================="<<endl;
	int nMoreJets = 0;
	bool foundMorejets = foundThirdJet(jet1Ind, jet2Ind, nMoreJets);


	if( foundMorejets ) 
	  continue;

	Photon *pho = (Photon*) branchPhoton->At(phoInd);

	Jet *jet1b = (Jet*) branchJet->At(jet1Ind);
	Jet *jet2b = (Jet*) branchJet->At(jet2Ind);
	
	if(debug){
	  cout<<"Found pho and jets, ind "<<phoInd<<" "<<jet1Ind<<" "<<jet2Ind<<endl;
	}
	
	
	TLorentzVector phov;
	phov.SetPtEtaPhiE(pho->PT, pho->Eta, pho->Phi, pho->E);

	///jet mass
	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jet1b->PT, jet1b->Eta, jet1b->Phi, jet1b->Mass);
	jet2.SetPtEtaPhiM(jet2b->PT, jet2b->Eta, jet2b->Phi, jet2b->Mass);


	
	TLorentzVector jetSys = (jet1 + jet2);

	hmap["phoPt"]->Fill(photon->PT);
	hmap["phoEta"]->Fill(photon->Eta);

	hmap["jet1Pt"]->Fill(jet1b->PT);
	hmap["jet2Pt"]->Fill(jet2b->PT);

	double jet1eta = jet1b->Eta;
	double jet2eta = jet2b->Eta;

	hmap["jet1Eta"]->Fill(jet1eta);
	hmap["jet2Eta"]->Fill(jet2eta);
	
	hmap["eta1TimesEta2"]->Fill(jet1eta*jet2eta);
	hmap["dEta"]->Fill(jet1eta-jet2eta);



	double jetmass = (jet1 + jet2).M();
	hmap["jetMass"]->Fill(jetmass); 

	hmap["jet1Mass"]->Fill(jet1.M()); 
	hmap["jet2Mass"]->Fill(jet2.M()); 

	double eta_pho_jetSys = pho->Eta - (jet1eta+jet2eta)/2.;
	hmap["eta_pho_jetSys"]->Fill(eta_pho_jetSys); 

	double eta_pho_jetSys_norm = -99;
	double deta = jet1eta-jet2eta; 
	eta_pho_jetSys_norm = ( pho->Eta - (jet1eta+jet2eta)/2. )/deta;

	hmap["eta_pho_jetSys_norm"]->Fill(eta_pho_jetSys_norm); 

	double exp_eta_pho_jetSys_norm = exp( -pow(eta_pho_jetSys_norm,2) );
	hmap["centrality"]->Fill(exp_eta_pho_jetSys_norm);

	double eta_pho_jet1 = pho->Eta - jet1eta;
	double eta_pho_jet2 = pho->Eta - jet2eta;
	
	hmap["eta_pho_jet1"]->Fill(eta_pho_jet1); 
	hmap["eta_pho_jet2"]->Fill(eta_pho_jet2); 

	
	
	double phi_pho_jet1 = phov.DeltaPhi(jet1);
	double phi_pho_jet2 = phov.DeltaPhi(jet2);

	hmap["phi_pho_jet1"]->Fill(phi_pho_jet1);
	hmap["phi_pho_jet2"]->Fill(phi_pho_jet2);


	double dphi_pho_jetmin = min(phi_pho_jet1,phi_pho_jet2);
	hmap["phi_pho_jetmin"]->Fill(dphi_pho_jetmin);

	double dphi_pho_jetmax = max(phi_pho_jet1,phi_pho_jet2);
	hmap["phi_pho_jetmax"]->Fill(dphi_pho_jetmax);

	double phi_jet1_jet2 = jet1.DeltaPhi(jet2);
	hmap["phi_jet1_jet2"]->Fill(phi_jet1_jet2);

	double phi_pho_jetSys = phov.DeltaPhi(jetSys);
	hmap["phi_pho_jetSys"]->Fill(phi_pho_jetSys);

	double dr_jet1jet2 = jet1.DeltaR(jet2);
	double dr_phojet1 = phov.DeltaR(jet1);
	double dr_phojet2 = phov.DeltaR(jet2);
	
	hmap["dr_jet1jet2"]->Fill(dr_jet1jet2);
	hmap["dr_phojet1"]->Fill(dr_phojet1);
	hmap["dr_phojet2"]->Fill(dr_phojet2);
	
	double mindR_phojet = min(dr_phojet1, dr_phojet2);
	hmap["dr_phojetmin"]->Fill(mindR_phojet);

	double maxdR_phojet = max(dr_phojet1, dr_phojet2);
	hmap["dr_phojetmax"]->Fill(maxdR_phojet);

	double pt_phoTojet1 = phov.Pt()/jet1.Pt();
	double pt_phoTojet2 = phov.Pt()/jet2.Pt();
	double pt_phoTojetSys = phov.Pt()/jetSys.Pt();

	hmap["pt_phoTojet1"]->Fill(pt_phoTojet1);
	hmap["pt_phoTojet2"]->Fill(pt_phoTojet2);
	hmap["pt_phoTojetSys"]->Fill(pt_phoTojetSys);

	TLorentzVector phoJetSys = (phov + jetSys);
	double mass3body = phoJetSys.M();
	hmap["mass3body"]->Fill(mass3body);

	////check which photon is closer to which jet (a bkg has photon as brem from the outgoing quarks) so form the 'brem' + jet sys and then take the pT ratio with the other jet

	TLorentzVector phoBremJet = (phov + jet1);
	double pt_phobrem_jet = phoBremJet.Pt()/jet2.Pt();
	
	if(dr_phojet2 < dr_phojet1)
	  {
	    phoBremJet = (phov + jet2);
	    pt_phobrem_jet = phoBremJet.Pt()/jet1.Pt();
	  }

	hmap["pt_phobrem_jet"]->Fill(pt_phobrem_jet);

	///29th Nov, 2021 - Variable got from momentum equation 
	
	double quark1_momeqn = jet1.Energy() * (1-cos(jet1.Theta()));

	double quark2_momeqn = jet2.Energy() * (1-cos(jet2.Theta()));

	double diff = (quark1_momeqn - quark2_momeqn)/(jet1.Energy()+jet2.Energy());
	hmap["momEquaVar"]->Fill(diff);
	//std::cout<<"momEquaVar "<<diff<<std::endl;
	////constituents inside the jet: https://github.com/delphes/delphes/blob/master/examples/Example3.C
	///https://cp3.irmp.ucl.ac.be/projects/delphes/ticket/999

	
	////make sure to call this for both the jets - used by selTracks
	vector<Track*> consti_jet1;
	vector<Tower*> towerconsti_jet1;
	getJetConstituents(jet1b,consti_jet1, towerconsti_jet1);

	vector<Track*> consti_jet2;
	vector<Tower*> towerconsti_jet2;
	getJetConstituents(jet2b,consti_jet2, towerconsti_jet2);

	
	// Loop over all jet's constituents ---> highest pT jet for now
	TLorentzVector momentum;
	
	//cout<<"Constitutes of jets "<<jet1b->Constituents.GetEntries()<<endl;
	bool foundMin1Const = false;
	

	for(int itk=0; itk<consti_jet1.size(); itk++){

	  Track* trk = (Track*)consti_jet1[itk];
	  momentum += trk->P4(); 
	  foundMin1Const = true;
	}

	double ratio_trackp_To_Jetp = momentum.P()/jet1.P();
	
	//cout<<"ratio_trackp_To_Jetp "<<ratio_trackp_To_Jetp<<endl;
	if(foundMin1Const) hmap["ratio_trackp_To_Jetp"]->Fill(ratio_trackp_To_Jetp);
	
	
	///nTracks between the two jets 

	
	double sumTrkPt = 0;
	double centrality_trk = -999;
	int ntracks = selTracks(jet1Ind, jet2Ind, consti_jet1, consti_jet2, towerconsti_jet1, towerconsti_jet2, sumTrkPt, hmap, centrality_trk);
	hmap["ntracks"]->Fill(ntracks);
	hmap["sumTrkPt"]->Fill(sumTrkPt);
	hmap["centrality_trk"]->Fill(centrality_trk);
	
	double dijet_scalarPt = jet1.Pt() + jet2.Pt();
	double dijet_vectorPt = (jet1+jet2).Pt();
	
	hmap["diJetScalarPt"]->Fill(dijet_scalarPt);
	hmap["diJetVectorPt"]->Fill(dijet_vectorPt);

	///Fox Wolfram Moments : 
	///1. https://arxiv.org/pdf/1212.4436.pdf
	///2. https://www.thphys.uni-heidelberg.de/~plehn/includes/theses/butter_b.pdf

	double fwm0 = calculateFWM(0, &jet1, &jet2, weightOptionFWM);
	double fwm1 = calculateFWM(1, &jet1, &jet2, weightOptionFWM);
	double fwm2 = calculateFWM(2, &jet1, &jet2, weightOptionFWM);

	double fwm3 = calculateFWM(3, &jet1, &jet2, weightOptionFWM);

	double fwm4 = calculateFWM(4, &jet1, &jet2, weightOptionFWM);

	double fwm5 = calculateFWM(5, &jet1, &jet2, weightOptionFWM);

	hmap["fwm0"]->Fill(fwm0);
	hmap["fwm1"]->Fill(fwm1);
	hmap["fwm2"]->Fill(fwm2);
	hmap["fwm3"]->Fill(fwm3);
	hmap["fwm4"]->Fill(fwm4);
	hmap["fwm5"]->Fill(fwm5);

	////Fill tree variables
	phoPt_ = photon->PT;
	phoEta_ = photon->Eta;
	jet1Eta_ = jet1b->Eta;
	jet2Eta_ = jet2b->Eta;
	jet1Pt_ = jet1b->PT;
	jet2Pt_ = jet2b->PT;
	jetEta1TimesEta2_ = jet1eta*jet2eta;
	jetdEta_ = jet1eta-jet2eta;
	jetMass_ = jetmass;

	dEta_phoJet1_ = eta_pho_jet1;
	dEta_phoJet2_ = eta_pho_jet2;
	
	dPhi_phoJet1_ = phi_pho_jet1;
	dPhi_phoJet2_ = phi_pho_jet2;
	dPhi_jet1Jet2_ = phi_jet1_jet2;
	dPhi_phoJetSys_ = phi_pho_jetSys;

	dPhi_phoJetMin_ = dphi_pho_jetmin;
	dPhi_phoJetMax_ = dphi_pho_jetmax;
	
	dR_jet1Jet2_ = dr_jet1jet2;
	dR_phoJet1_ = dr_phojet1;
	dR_phoJet2_ = dr_phojet2;
	
	dR_phoJetMin_ = mindR_phojet;
	dR_phoJetMax_ = maxdR_phojet;
	

	pT_phoTojet1_  = pt_phoTojet1;
	pT_phoTojet2_  = pt_phoTojet2;
	pT_phoTojetSys_  = pt_phoTojetSys;

	invMass_3body_ = mass3body;
	pT_phoBremJet_ = pt_phobrem_jet;

	zeppenfeld_ = eta_pho_jetSys;
	zeppenfeldNorm_ = eta_pho_jetSys_norm;
	nTracks_ = ntracks;
	sumTrkPt_ = sumTrkPt;
	
	centrality_ = exp_eta_pho_jetSys_norm;
	diJetScalarPt_ = dijet_scalarPt;
	diJetVectorPt_ = dijet_vectorPt;
	
	centrality_pt_trk_ = centrality_trk;

	///FWMs
	fwm0_ = fwm0;
	fwm1_ = fwm1;
	fwm2_ = fwm2;
	fwm3_ = fwm3;
	fwm4_ = fwm4;
	fwm5_ = fwm5;

	//isTrain
	double ran = gran.Rndm();
	if(ran<fracTrain)
	  isTrain_ = 1;
	else
	  isTrain_ = 0;
	///                  

	tout->Fill();
      }//if(phoInd>=0 && jet1Ind>=0 && jet2Ind>=0)



      
   }

  for(map<string,TH1F*>::iterator it = hmap.begin(); it != hmap.end(); ++it) {
    hmap[it->first]->Write();
  }

  tout->Write();

}


int Analyse::selPho(){

  int phoind = -99;
  double maxPt = -99;
  
  for(int ipho=0; ipho<branchPhoton->GetEntriesFast(); ipho++){

    photon = (Photon*) branchPhoton->At(ipho);
    
    //// skip photons with references to multiple particles
    //if(photon->Particles.GetEntriesFast() != 1) continue;

    
    if( fabs(photon->Eta) > 2.5 ) continue;
    //if( Photon_EhadOverEem[ipho] > 0.05 ) continue;
    if( photon->PT < 10. ) continue;
    //if( photon->PT < 20. ) continue;
    
    if(maxPt < photon->PT)
      {
	maxPt = photon->PT;
	phoind = ipho;
      }
  }

  return phoind;
}


void Analyse::selVBFJets(int &jet1ind, int &jet2ind){

  jet1ind = -99;
  jet2ind = -99;
  
  double max1Pt = -99.;
  double max2Pt = -99.;

  ///1st jet selection
  for(int ijet=0; ijet<branchJet->GetEntriesFast(); ijet++){

    jet = (Jet*) branchJet->At(ijet);
    if( fabs(jet->Eta) > 4.7 ) continue;
    if( jet->PT < 30 ) continue;
      
    if(max1Pt < jet->PT)
      {
	max1Pt = jet->PT;
	jet1ind = ijet;
      }
  }///for(int ijet=0; ijet<Jet_; ijet++)

  ///2nd jet selection
  for(int ijet=0; ijet<branchJet->GetEntriesFast(); ijet++){

    if(ijet==jet1ind) continue;
    
    if( fabs(jet->Eta) > 4.7 ) continue;
    if( jet->PT < 30 ) continue;
    
    if(max2Pt < jet->PT)
      {
	max2Pt = jet->PT;
	jet2ind = ijet;
      }
  }///for(int ijet=0; ijet<Jet_; ijet++)

  
}

int Analyse::selTracks(int jet1Ind, int jet2Ind, vector<Track*> consti_jet1, vector<Track*> consti_jet2, vector<Tower*> towerconsti_jet1, vector<Tower*> towerconsti_jet2, double &sumTrkPt, map<string, TH1F*> hmap, double &centrality_trk){

  centrality_trk = 0;
  
  int nTrks = 0;

  sumTrkPt = 0;
  
  Jet *jet1 = (Jet*) branchJet->At(jet1Ind); 
  Jet *jet2 = (Jet*) branchJet->At(jet2Ind); 

  ///jet mass
  TLorentzVector jet1v, jet2v;
  jet1v.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
  jet2v.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);

  double avJetEta = (jet1v.Eta() + jet2v.Eta())/2.;
  
  TLorentzVector momentum1, momentum2;
  
  /////////////////////Track consti of jet 1
  //cout<<"#track consti in jet 1 "<<consti_jet1.size()<<endl;
  for(int itk=0; itk<consti_jet1.size(); itk++){
    
    Track* trk = (Track*)consti_jet1[itk];

    momentum1 += trk->P4();

    //cout<<"For jet 1, this constituent track eta and Pt "<<trk->Eta<<" "<<trk->PT<<endl;
    
  }


  ///////////////////////Tower consti of jet 1
  //cout<<"#tower consti in jet 1 "<<towerconsti_jet1.size()<<endl;
  for(int itk=0; itk<towerconsti_jet1.size(); itk++){
    
    Tower *tower = (Tower*)towerconsti_jet1[itk];

    momentum1 += tower->P4();

    //cout<<"For jet 1, this constituent tower eta and Pt "<<tower->Eta<<" "<<tower->ET<<endl;
    
  }

  /////////////////////Track consti of jet 2
  //cout<<"#track consti in jet 2 "<<consti_jet2.size()<<endl;
  for(int itk=0; itk<consti_jet2.size(); itk++){
    
    Track* trk = (Track*)consti_jet2[itk];

    momentum2 += trk->P4();
    //cout<<"For jet 2, this constituent track eta and Pt "<<trk->Eta<<" "<<trk->PT<<endl;

  }

  ///////////////////////Tower consti of jet 2
  //cout<<"#tower consti in jet 2 "<<towerconsti_jet2.size()<<endl;
  for(int itk=0; itk<towerconsti_jet2.size(); itk++){
    
    Tower *tower = (Tower*)towerconsti_jet2[itk];

    momentum2 += tower->P4();

    //cout<<"For jet 2, this constituent tower eta and Pt "<<tower->Eta<<" "<<tower->ET<<endl;
    
  }

  //cout<<"Total summed up momentum of jet1 : jet1 pt : Total summed up momentum of jet2 : jet2 pt "<<momentum1.Pt()<<" "<<jet1->PT<<" "<<momentum2.Pt()<<" "<<jet2->PT<<endl;

  //cout<<"Now tracks collection variables"<<endl;

  
  for(int itrk=0; itrk<branchTrack->GetEntries(); itrk++){
    
    //cout<<"inside trk loop i trk "<<itrk<<endl;
    
    track = (Track*) branchTrack->At(itrk);
    
    //std::cout<<" track eta : jet 1 eta : jet 2 eta : "<<track->Eta<<" "<<jet1->Eta<<" "<<jet2->Eta<<std::endl;

    bool sel_midTracks = (track->Eta > jet1->Eta && track->Eta<jet2->Eta) || (track->Eta>jet2->Eta && track->Eta<jet1->Eta);
    
    if( !sel_midTracks ) continue;

    TLorentzVector *tr = new TLorentzVector();
    tr->SetPtEtaPhiM(track->PT, track->Eta, track->Phi, track->Mass);
    
    double dR_jet1 = tr->DeltaR(jet1v);
    double dR_jet2 = tr->DeltaR(jet2v);
    
    ///also make sure that they are 1 unit away in eta from each jet
    //sel_midTracks = dR_jet1 > 4 && dR_jet2 > 4;  ///this makes #tracks in noISRFSR in sig is <= 3

    //sel_midTracks = dR_jet1 > 0.5 && dR_jet2 > 0.5;  ///this makes #tracks in noISRFSR in sig is <= 3

    if( !sel_midTracks ) continue;
    
    //if( Photon_EhadOverEem[itrk] > 0.05 ) continue;
    //if( track->PT < 0.5 ) continue;
    if( track->PT < 1 ) continue;


    ////now check if the track is a component of any of the jets
    ///jet1
    bool selTrk = true;
    for(int icon=0; icon<consti_jet1.size(); icon++){
    
    Track* trkconst = (Track*)consti_jet1[icon];

    double dEta = fabs(trkconst->Eta - track->Eta);
    if(dEta < 0.009){
      selTrk = false;
      break;
    }
    
    }///for(int icon=0; icon<consti_jet1.size(); icon++)

    if(!selTrk) continue;

    ///jet 2
    selTrk = true;
    for(int icon=0; icon<consti_jet2.size(); icon++){
    
    Track* trkconst = (Track*)consti_jet2[icon];

    double dEta = fabs(trkconst->Eta - track->Eta);
    if(dEta < 0.009){
      selTrk = false;
      break;
    }
    
    }///for(int icon=0; icon<consti_jet1.size(); icon++)

    if(!selTrk) continue;
    
    double dEta1 = track->Eta - jet1->Eta;
    double dEta2 = track->Eta - jet2->Eta;
    
    //std::cout<<" Selected track eta - track PT : jet1Eta : jet2Eta : dEta1 : dEta2 : "<<track->Eta<<" "<<track->PT<<" "<<jet1->Eta<<" "<<jet2->Eta<<" "<<dEta1<<" "<<dEta2<<std::endl;

    nTrks++;
    sumTrkPt += track->PT;
    
    //cout<<"track eta : track pt "<<track->Eta<<" "<<track->PT<<endl;
    
    hmap["trkEta"]->Fill(track->Eta);
    hmap["trkPt"]->Fill(track->PT);
    
    
    centrality_trk += track->PT * exp( -(track->Eta - avJetEta) );
    //double centrality_trk = 
  }

  
  return nTrks;
}


void Analyse::getJetConstituents(Jet *jet, vector<Track*> &consti, vector<Tower*> &towerconsti){
  
  bool foundMin1Const = false;
  
  for(int jc = 0; jc < jet->Constituents.GetEntriesFast(); ++jc)
    {
      TObject *object = jet->Constituents.At(jc);
      
      // Check if the constituent is accessible
      if(object == 0) continue;
      
      
      if(object->IsA() == Track::Class())
	{
	  Track* trk = (Track*) object;
	  //momentum += trk->P4();
	  consti.push_back(trk);
	}
     
      else if(object->IsA() == Tower::Class())
	{
	  Tower *tower = (Tower*) object;
	  //momentum += trk->P4();
	  towerconsti.push_back(tower);
	}
      
      
    }//for(int jc = 0; jc < jet1b->Constituents.GetEntriesFast(); 
  
}


///check if there is a third jet
bool Analyse::foundThirdJet(int jet1ind, int jet2ind, int &nMoreJets){

  nMoreJets = 0;
  bool found3jet = false;
  ///1st jet selection
  for(int ijet=0; ijet<branchJet->GetEntriesFast(); ijet++){

    //std::cout<<"iJet : jet1Ind : jet2Ind "<<ijet<<" "<<jet1ind<<" "<<jet2ind<<std::endl;
    
    if(ijet==jet1ind || ijet==jet2ind) continue;
    
    
    jet = (Jet*) branchJet->At(ijet);
    if( fabs(jet->Eta) > 4.7 ) continue;
    //if( jet->PT < 30 ) continue;
    //if( jet->PT < 10 ) continue;
    if( jet->PT < 0 ) continue;
    
    nMoreJets += 1;

    found3jet = true;
  }///for(int ijet=0; ijet<Jet_; ijet++)

  //std::cout<<"nAdditional Jets in the even "<<nMoreJets<<std::endl;
	
  return found3jet;
}


double Analyse::legendrePolynomial(int l, double x){

  double P0 = 1;
  double P1 = x;

  if(l==0) return P0;
  if(l==1) return P1;

  if(l>1){
    double f = ( (2*l-1) * x * legendrePolynomial(l-1, x) - (l-1) * legendrePolynomial(l-2, x) )/l;
    return f;
  }

}


double Analyse::weightFactorFWM(string option, TLorentzVector *p1, TLorentzVector *p2){

  TLorentzVector p = (*p1) + (*p2);
  
  double weight = 1;
  
  ///https://arxiv.org/pdf/1212.4436.pdf - equation 4
  switch(str2int(option.c_str())){
    
  case str2int("s") :
    weight = (p1->P() * p2->P());
    return weight;
    
  case str2int("p"):
    weight = (p1->P() * p2->P());
    return weight;
    
  case str2int("T"):
    weight = (p1->Pt() * p2->Pt());
    return weight;
    
  case str2int("z"):
    weight = (p1->Pz() * p2->Pz());
    return weight;

  case str2int("y"):{
    double eta1 = p1->Eta();
    double eta2 = p2->Eta();
    double y_avg = (eta1 + eta2)/2.;
     weight = pow(eta1-y_avg, -1) * pow(eta2-y_avg, -1);
    return weight;
  }

  case str2int("1"):
    return 1;

  }

}

double Analyse::denominatorFWM(string option, TLorentzVector *p1, TLorentzVector *p2){

  TLorentzVector p = (*p1) + (*p2);

  double weight = 1;
  
  ///https://arxiv.org/pdf/1212.4436.pdf - equation 4
  switch(str2int(option.c_str())){
    
  case str2int("s") :
    weight = pow(p.P(),2);
    return weight;
    
  case str2int("p"):
    weight = pow(p1->P() + p2->P(),2);
    return weight;
    
  case str2int("T"):
    weight = pow(p.Pt(),2);
    return weight;
    
  case str2int("z"):
    weight = pow(p.Pz(),2);
    return weight;
  

  case str2int("y"):{
    double eta1 = p1->Eta();
    double eta2 = p2->Eta();
    double y_avg = (eta1 + eta2)/2.;
    weight = pow( pow(eta1-y_avg, -1) + pow(eta2-y_avg, -1), 2 );
    return weight;
  }

  case str2int("1"):
    return 1;

  }


}


double Analyse::calculateFWM(int l, TLorentzVector *p1, TLorentzVector *p2, string weightOption){

  double theta1 = p1->Theta();
  double theta2 = p2->Theta();
  double phi1 = p1->Phi();
  double phi2 = p2->Phi();
  double dPhi = deltaPhi(phi1,phi2);
  

  double cosOmega = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(dPhi);

  ///eqn 36 from here for 2 jet case: https://www.thphys.uni-heidelberg.de/~plehn/includes/theses/butter_b.pdf 
  double fwm11 = weightFactorFWM(weightOption, p1, p1) * legendrePolynomial(l,1); ///same particle so angle cosOmega = 1
  double fwm12 = weightFactorFWM(weightOption, p1, p2) * legendrePolynomial(l,cosOmega);
  double fwm22 = weightFactorFWM(weightOption, p2, p2) * legendrePolynomial(l,1); ///same particle so angle cosOmega = 1
  double fwm = (fwm11 + 2*fwm12 + fwm22)/denominatorFWM(weightOption,p1,p2);

  //cout<<"p1 : p2 "<<p1->P()<<" "<<p2->P()<<endl;
  //cout<<"fwm11 : fwm12 : fwm22 : deno "<<fwm11<<" "<<fwm12<<" "<<fwm22<<" deno "<<denominatorFWM(weightOption,p1,p2)<<endl;
  //cout<<"legendre(l,1) : weight11 : deno : "<<legendrePolynomial(l,1)<<" "<<weightFactorFWM(weightOption, p1, p1)<<" "<<denominatorFWM(weightOption,p1,p2)<<endl;
  
  //cout<<"fwm "<<fwm<<endl;

  return fwm;

}


TChain* Analyse::list_files(const char *dirname, const char *ext){
  
  TChain *t = new TChain("Delphes");
  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 
  if (files) 
    { 
      TSystemFile *file; 
      TString fname; 
      TIter next(files); 
      while ((file=(TSystemFile*)next())) { 
	fname = file->GetName(); 
	if (!file->IsDirectory() && fname.EndsWith(ext)) 
	  { 
	    cout << fname.Data() << endl; 
	    t->Add(Form("%s/%s",dirname,fname.Data()));
	  } 
      } 
    } 
  return t;
}


