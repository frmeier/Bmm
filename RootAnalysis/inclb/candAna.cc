#include "candAna.hh"

#include "common/HFMasses.hh"


using namespace std;

// ----------------------------------------------------------------------
candAna::candAna(inclbReader *pReader, string name, string cutsFile) {
  fpReader = pReader; 
  fVerbose = fpReader->fVerbose;	 
  fYear    = fpReader->fYear; 
  fName    = name; 
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", reading cutsfile " << cutsFile << " setup for year " << fYear << endl;
  fHistDir = gFile->mkdir(fName.c_str());
  readCuts(cutsFile, 1); 
  cout << "======================================================================" << endl;	 
}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}

// ----------------------------------------------------------------------
void candAna::endAnalysis() {
  TH1D *h1 = ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())));
  if (h1) {
    cout << Form("==> mon%s: events seen    = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(2))) << endl;
    cout << Form("==> mon%s: cands analysed = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(11))) << endl;
    cout << Form("==> mon%s: cands passed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(12))) << endl;
    cout << Form("==> mon%s: cands failed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(21))) << endl;
    if (h1->GetBinContent(2) < 1) {
      cout << Form("==> mon%s: error, no events seen!", fName.c_str()) << endl; 
    }
  } else {
    cout << Form("==> mon%s: error, histogram not found!", fName.c_str()) << endl; 
  }    
}



// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {
  
  fpEvt = evt; 

  TAnaTrack *pSigTrack(0);
  if (fVerbose == -66) { 
    cout << "---- Evt: " << fEvt << " n(sig tracks) = " << fpEvt->nSigTracks() << endl;
  }
  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(1);
  fpSigTrack = 0; 
  fpSigJet = 0; 
  for (int iC = 0; iC < fpEvt->nSigTracks(); ++iC) {
    pSigTrack = fpEvt->getSigTrack(iC);
    
    if (TYPE != pSigTrack->fInt1) {
      if (fVerbose > 39) cout << "  skipping sig track at " << iC << " which is of type " << pSigTrack->fInt1 <<endl;
      continue;
    }
    
    if (pSigTrack->fInt2 < 0) continue;

    fpSigTrack = pSigTrack; 
    fpSigJet   = fpEvt->getTrackJet(pSigTrack->fInt2);

    if (fVerbose > 99) {
      cout << "Analyzing sig track at " << iC << " which is of type " << fpSigTrack->fInt1 << " pt " << fpSigTrack->fPlab.Perp()
	   << " ptrel = " << fpSigTrack->fDouble1
	   << " and jet index: " << pSigTrack->fInt2
	   << " jet pt/eta = " << fpSigJet->fPlab.Perp() 
	   << "/" <<  fpSigJet->fPlab.Eta() 
	   << endl;
    }

    candAnalysis();

    fTree->Fill(); 
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(31);
  } 

}

// ----------------------------------------------------------------------
void candAna::candAnalysis() {

  if (!fpSigTrack) return;
  if (!fpSigJet) return;

  //((TH1D*)fHistDir->Get("../monEvents"))->Fill(1); 

  fType = fpSigTrack->fInt1; 

  fMuPtRel = fpSigTrack->fDouble1; 

  fMuPt    = fpSigTrack->fPlab.Perp(); 
  fMuEta   = fpSigTrack->fPlab.Eta(); 
  fMuPhi   = fpSigTrack->fPlab.Phi(); 

  fJetPt   = fpSigJet->fPlab.Perp(); 
  fJetEta  = fpSigJet->fPlab.Eta(); 
  fJetPhi  = fpSigJet->fPlab.Phi(); 

  fillRedTreeData();

}


// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  TH1D *h11(0); 
  (void)h11; 
  h11 = new TH1D(Form("mon%s", fName.c_str()), Form("mon%s", fName.c_str()), 50, 0., 50.); 


  TH2D *h22(0); 
  (void)h22; 

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupRedTree(fTree); 
  
  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.); 
  (void)h;

}


// ----------------------------------------------------------------------
void candAna::setupRedTree(TTree *t) {

  t->Branch("run",     &fRTD.run,          "run/L");
  t->Branch("evt",     &fRTD.evt,          "evt/L");
  t->Branch("ls",      &fRTD.ls,           "ls/I");

  t->Branch("muid",    &fRTD.muid,         "muid/O");
  t->Branch("hlt",     &fRTD.hlt,          "hlt/O");
  t->Branch("hltmatch",&fRTD.hltmatch,     "hltmatch/O");
  t->Branch("json",    &fRTD.json,         "json/O");

  t->Branch("type",    &fRTD.type,         "type/I");

  t->Branch("pt",      &fRTD.pt,           "pt/F");
  t->Branch("eta",     &fRTD.eta,          "eta/F");
  t->Branch("phi",     &fRTD.phi,          "phi/F");
  t->Branch("ptrel",   &fRTD.ptrel,        "ptrel/F");

  t->Branch("jpt",     &fRTD.jpt,          "jpt/F");
  t->Branch("jeta",    &fRTD.jeta,         "jeta/F");
  t->Branch("jphi",    &fRTD.jphi,         "jphi/F");

}


// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  // -- define default values for some cuts
  NOPRESELECTION = 0; 
  IGNORETRIGGER  = 0; 

  // -- set up cut sequence for analysis
  fCutFile = fileName;

  if (dump) cout << "==> candAna: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[1000];
  fHistDir->cd();
  if (dump) cout << "gDirectory: "; fHistDir->pwd();
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); 
      if (dump) cout << "TYPE:           " << TYPE << endl;
      if (100100 == TYPE) cstring = "TrackJets";
      ibin = 1;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "TRIGRANGE")) {
      char triggerlist[1000]; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      int r1(0), r2(0); 
      string hlt = splitTrigRange(tl, r1, r2); 
      HLTRANGE.insert(make_pair(hlt, make_pair(r1, r2))); 
      if (dump) {
	cout << "HLTRANGE:       " << hlt << " from " << r1 << " to " << r2 << endl; 
      }
      ibin = 3; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }


    if (!strcmp(CutName, "IGNORETRIGGER")) {
      IGNORETRIGGER = int(CutValue); 
      if (dump) cout << "IGNORETRIGGER      " << IGNORETRIGGER << endl;
      ibin = 5;
      hcuts->SetBinContent(ibin, IGNORETRIGGER);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore trigger :: %i", CutName, IGNORETRIGGER));
    }

    if (!strcmp(CutName, "NOPRESELECTION")) {
      NOPRESELECTION = int(CutValue); 
      if (dump) cout << "NOPRESELECTION     " << NOPRESELECTION << endl;
      ibin = 6;
      hcuts->SetBinContent(ibin, NOPRESELECTION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore preselection :: %i", CutName, NOPRESELECTION));
    }


    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; 
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 202;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(#mu) :: %3.1f", CutName, MUPTLO));
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; 
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 203;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(#mu) :: %3.1f", CutName, MUPTHI));
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; 
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 204;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(#mu) :: %3.1f", CutName, MUETALO));
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; 
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 205;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(#mu) :: %3.1f", CutName, MUETAHI));
    }

  }

  if (dump)  cout << "------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void candAna::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    exit(1);
  }
  char input[1000]; 
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] != '+') {
      lines.push_back(string(buffer));
    } else {
      sscanf(buffer, "+input %s", input);
      readFile(input, lines); 
    }
  }

}


// ----------------------------------------------------------------------
void candAna::getSigTracks(vector<int> &v, TAnaCand *pC) {
  TAnaCand *pD; 
  TAnaTrack *pT; 
  vector<int> bla; 

  // -- loop over daughters
  if (pC->fDau1 > -1) {
    for (int j = pC->fDau1; j <= pC->fDau2; ++j) {
      pD = fpEvt->getCand(j); 
      getSigTracks(bla, pD); 
    }

    for (unsigned j = 0; j < bla.size(); ++j) v.push_back(bla[j]);
  }

  // -- add direct sigtracks
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    if (v.end() == find(v.begin(), v.end(), pT->fIndex)) {
      v.push_back(pT->fIndex); 
    }
  }

}


// ----------------------------------------------------------------------
void candAna::fillRedTreeData() {
  fRTD.run       = fRun;
  fRTD.evt       = fEvt;
  fRTD.ls        = fLS;

  fRTD.muid      = fMuId;
  fRTD.hlt       = fHLT;
  fRTD.hltmatch  = fHLTmatch;
  fRTD.json      = fJSON;

  fRTD.pt        = fMuPt; 
  fRTD.eta       = fMuEta; 
  fRTD.phi       = fMuPhi; 
  fRTD.ptrel     = fMuPtRel; 

}


// ----------------------------------------------------------------------
string candAna::splitTrigRange(string tl, int &r1, int &r2) {

  string::size_type id1 = tl.find_first_of("("); 
  string::size_type id2 = tl.find_first_of(":"); 
  string::size_type id3 = tl.find_first_of(")"); 

  //cout << "tl: " << tl << endl;
  string hlt = tl.substr(0, id1);
  //cout << "hlt: " << hlt << endl;
  string a   = tl.substr(id1+1, id2-id1-1);
  r1 = atoi(a.c_str());
  //cout << "1st a: " << a << " -> r1 = " << r1 << endl;
  a  = tl.substr(id2+1, id3-id2-1); 
  r2 = atoi(a.c_str());
  //cout << "2nd a: " << a << " -> r2 = " << r2 << endl;

  return hlt; 

}
