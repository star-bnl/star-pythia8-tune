// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
//! c++ headers
#include <iostream>     
#include <algorithm>    
#include <vector>

namespace Rivet {

  class STAR_2012_I930463: public Analysis {
    public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I930463);

      /// Book histograms and initialise projections before the run
      void init() {

        // Initialise and register projections
        const ChargedFinalState cfs(Cuts::abseta < 1.5 && Cuts::abscharge > 0);
        declare(cfs, "CFS");
        const FinalState fs(Cuts::abseta < 1.5);
        declare(fs, "FS");

        // Book histograms
        book(_h_pT_piplus ,1, 1, 1);
        book(_h_pT_piminus,1, 1, 2);
        book(_h_pT_kplus  ,1, 1, 3);
        book(_h_pT_kminus ,1, 1, 4);
        book(_h_pT_p      ,1, 1, 5);
        book(_h_pT_pbar   ,1, 1, 6);
        // book(_h_pT_Kshort ,2, 1, 1);
        // book(_h_pT_pho    ,3, 1, 1);
        book(_h_ratio_pimpip ,7,1,1);
        book(_h_ratio_pbarp ,7,1,2);
        book(_h_ratio_KmKp ,7,1,3);
        book(_h_ratio_ppi ,7,1,4);
        book(_h_ratio_pbarpibar ,7,1,5);
        book(_h_ratio_Kpi ,7,1,6);
        book(_h_K,"TMP/K+K-" ,refData(7,1,6));
        book(_h_pi,"TMP/pi+pi-" ,refData(7,1,6));

        book(_sumWeightSelected, "TMP/sumWselected");

      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const FinalState& finalstate = applyProjection<FinalState>(event, "FS");

        for(Particle p: finalstate.particles()){

          if (p.absrap() < 0.5 && p.pT()/GeV>3) {
            const PdgId pid = p.pid();
            const double pT = p.pT() / GeV;
            switch (abs(pid)) {
              case PID::PIPLUS:
                _h_pi->fill(pT, (1./(2.0 * M_PI))/pT);
                pid > 0 ? _h_pT_piplus->fill(pT, (1./(2.0 * M_PI))/pT) : _h_pT_piminus->fill(pT, (1./(2.0 * M_PI))/pT);
                break;
              case PID::PROTON:
                pid > 0 ? _h_pT_p->fill(pT, (1./(2.0 * M_PI))/pT) : _h_pT_pbar->fill(pT, (1./(2.0 * M_PI))/pT);
                break;
              case PID::KPLUS:
                _h_K->fill(pT, 1.0/pT);
                pid > 0 ? _h_pT_kplus->fill(pT, (1./(2.0 * M_PI))/pT) : _h_pT_kminus->fill(pT, (1./(2.0 * M_PI))/pT);
                break;
            }
          }

        } //end of loop FS particle 

        _sumWeightSelected->fill();
      } // end event ana

      /// Normalise histograms etc., after the run
      void finalize() {
	
        double ppCrossSection = 30*1e3*microbarn; // STAR pp MB 

        const YODA::Scatter1D factor = 1. / *_sumWeightSelected;

        scale(_h_pT_piplus,       factor);
        scale(_h_pT_piminus,       factor);
        scale(_h_pT_kminus,       factor);
        scale(_h_pT_kplus,       factor);
        scale(_h_pT_p,       factor);
        scale(_h_pT_pbar,       factor);

        divide(_h_pT_piminus, _h_pT_piplus, _h_ratio_pimpip);
        divide(_h_pT_pbar, _h_pT_p, _h_ratio_pbarp);
        divide(_h_pT_kminus, _h_pT_kplus, _h_ratio_KmKp);
        divide(_h_pT_p, _h_pT_piplus, _h_ratio_ppi);
        divide(_h_pT_pbar, _h_pT_piminus, _h_ratio_pbarpibar);
        divide(_h_K, _h_pi, _h_ratio_Kpi);

      }

    private:
      CounterPtr _sumWeightSelected;
      Histo1DPtr _h_pT_piplus, _h_pT_piminus, _h_pT_kminus, _h_pT_kplus, _h_pT_p, _h_pT_pbar, _h_K,_h_pi;
      Scatter2DPtr _h_ratio_pimpip, _h_ratio_pbarp, _h_ratio_ppi, _h_ratio_KmKp, _h_ratio_pbarpibar, _h_ratio_Kpi;
  };


  DECLARE_RIVET_PLUGIN(STAR_2012_I930463);

}
