#include "HFLb2JpsiL0.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"

#include <iostream>

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFLb2JpsiL0::HFLb2JpsiL0(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
  fL0Window(iConfig.getUntrackedParameter<double>("L0Window", 0.2)),
  fLbWindow(iConfig.getUntrackedParameter<double>("LbWindow", 0.8)) {
  dumpConfiguration();
} // HFLb2JpsiL0()

// ----------------------------------------------------------------------
void HFLb2JpsiL0::dumpConfiguration() {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFLb2JpsiL0 configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  L0Window:                " << fL0Window << endl;
  cout << "---  LbWindow:                 " << fLbWindow << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()


// ----------------------------------------------------------------------
void HFLb2JpsiL0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  }
  catch (HFSetupException e) {
    cout << "==>HFLb2JpsiL0> " << e.fMsg << endl;
    return;
  }
	
  fListBuilder->setMinPt(fMuonPt); // work with muon pt and not with track pt
  vector<int> muonList = fListBuilder->getMuonList();
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();
	
  if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList),MMUON,
						     (fPsiMuons < 2 ? trkList : muonList),MMUON,
						     MJPSI-fPsiWindow,MJPSI+fPsiWindow, 1);
  HFTwoParticleCombinatoricsSet L0List = a.combine(trkList,MPROTON,trkList,MPION,MLAMBDA_0-fL0Window,MLAMBDA_0+fL0Window, 0);
	
  if (fVerbose > 0) cout << "==>HFLb2JpsiL0> J/psi list size: " << psiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFLb2JpsiL0> L0 list size: " << L0List.size() << endl;
	
  // -- Build J/psi + L0
  TLorentzVector psi, L0, m1, m2, pr, pi, lb;
  for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
    unsigned int iMuon1 = psiIt->first;
    unsigned int iMuon2 = psiIt->second;
		
    reco::TrackBaseRef mu1TrackView(fTracksHandle, iMuon1);
    reco::Track tMuon1(*mu1TrackView);
    if (tMuon1.pt() < fMuonPt)  continue;
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);
		
    reco::TrackBaseRef mu2TrackView(fTracksHandle, iMuon2);
    reco::Track tMuon2(*mu2TrackView);
    if (tMuon2.pt() < fMuonPt)  continue;
    m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON);
		
    psi = m1 + m2;
    if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
    for (HFTwoParticleCombinatoricsNew::iterator L0It = L0List.begin(); L0It != L0List.end(); ++L0It) {
      unsigned int iProton = L0It->first;
      unsigned int iPion = L0It->second;
			
      if (iProton == iMuon1 || iProton == iMuon2) continue;
      if (iPion == iMuon1 || iPion == iMuon2) continue;

      reco::TrackBaseRef rTrackView1(fTracksHandle, iProton);
      reco::Track tProton(*rTrackView1);
      if (tProton.pt() < fTrackPt) continue;
      pr.SetXYZM(tProton.px(), tProton.py(), tProton.pz(), MKAON); 
      if (psi.DeltaR(pr) > fDeltaR) continue;
			
      reco::TrackBaseRef rTrackView2(fTracksHandle, iPion);
      reco::Track tPion(*rTrackView2);
      if (tPion.pt() < fTrackPt) continue;
      pi.SetXYZM(tPion.px(), tPion.py(), tPion.pz(), MKAON); 
      if (psi.DeltaR(pi) > fDeltaR) continue;
			
      L0 = pr + pi; 
      if ((TMath::Abs(L0.M() - MLAMBDA_0) > fL0Window)) continue;
			
      lb = psi + L0; 
      if (TMath::Abs(lb.M() - MLAMBDA_B) > fLbWindow) continue;
			
      // -- sequential fit: J/Psi kaons
      HFDecayTree theTree(305122, true, MLAMBDA_B, false, -1.0, true);
			
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443, true, MJPSI, false); // Don't use kinematic particle for the Psi
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      iterator = theTree.addDecayTree(303122, true, MLAMBDA_0, false);
      iterator->addTrack(iProton,321);
      iterator->addTrack(iPion,321);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      fSequentialFitter->doFit(&theTree);
			
      // -- sequential fit: J/Psi (constraint) L0 (unconstraint)
      theTree.clear(405122, true, MLAMBDA_B, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      iterator = theTree.addDecayTree(403122, true, MLAMBDA_0, false); // Note: this differs from the lifetime analysis, where we constrained the L0 mass
      iterator->addTrack(iProton,321);
      iterator->addTrack(iPion,321);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      fSequentialFitter->doFit(&theTree);
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFLb2JpsiL0);
