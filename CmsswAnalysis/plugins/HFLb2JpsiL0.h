#ifndef _HFLB2JPSIL0_h_
#define _HFLB2JPSIL0_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFLb2JpsiL0 : public HFVirtualDecay {
	public:
		explicit HFLb2JpsiL0(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();
		
		int           fPsiMuons;
		double        fPsiWindow, fL0Window, fLbWindow;
};

#endif
