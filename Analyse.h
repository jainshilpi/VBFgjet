//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 21 00:28:07 2021 by ROOT version 6.22/06
// from TTree Delphes/Analysis tree
// found on file: signal1.root
//////////////////////////////////////////////////////////

#ifndef Analyse_h
#define Analyse_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TRef.h"
// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TRefArray.h"

//#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"

//#else
class ExRootTreeReader;
class ExRootResult;
//#endif

#include "TH1.h"

#include <iostream>

class Analyse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TClonesArray *branchParticle;
   TClonesArray *branchElectron;
   TClonesArray *branchPhoton;
   TClonesArray *branchMuon;
   TClonesArray *branchEFlowTrack;
   TClonesArray *branchEFlowPhoton;
   TClonesArray *branchEFlowNeutralHadron;
   TClonesArray *branchJet;
   TClonesArray *branchTrack;

  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;
 

   Analyse(TTree *tree=0);
   virtual ~Analyse();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(ExRootTreeReader *treeReader);
   virtual void     Loop(string fname, string fnameout);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual TChain*  list_files(const char *dirname="C:/root/folder/", const char *ext=".root");
   
   virtual int selPho();
   virtual void     selVBFJets(int &jet1ind, int &jet2ind); 
   virtual int      selTracks(int jet1Ind, int jet2Ind, vector<Track*> consti_jet1, vector<Track*> consti_jet2, vector<Tower*> towerconsti_jet1, vector<Tower*> towerconsti_jet2, double &sumTrkPt,map<string, TH1F*> hmap, double &centrality_trk );
   virtual void     getJetConstituents(Jet *jet, vector<Track*> &consti, vector<Tower*> &towerconsti);
   
   virtual bool     foundThirdJet(int jet1ind, int jet2ind, int &nMoreJets);

   virtual double calculateFWM(int l, TLorentzVector *p1, TLorentzVector *p2, string weightOption);
   
   virtual double denominatorFWM(string option, TLorentzVector *p1, TLorentzVector *p2);

   virtual double weightFactorFWM(string option, TLorentzVector *p1, TLorentzVector *p2);

   virtual double legendrePolynomial(int l, double x);

};

#endif

#ifdef Analyse_cxx
Analyse::Analyse(TTree *tree) : fChain(0) 
{
std::cout<<"initialized the Analyse class"<<std::endl;
  /*
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("signal1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("signal1.root");
      }
      f->GetObject("Delphes",tree);

   }
   Init(tree);
  */
}

Analyse::~Analyse()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyse::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyse::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyse::Init(ExRootTreeReader *treeReader)
{

  branchParticle = treeReader->UseBranch("Particle");
  branchElectron = treeReader->UseBranch("Electron");
  branchPhoton = treeReader->UseBranch("Photon");
  branchMuon = treeReader->UseBranch("Muon");
  branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  branchJet = treeReader->UseBranch("Jet");
  branchTrack = treeReader->UseBranch("Track");

   Notify();
}

Bool_t Analyse::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyse::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyse::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyse_cxx
