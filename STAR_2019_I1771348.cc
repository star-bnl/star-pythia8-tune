// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2019_I1771348 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2019_I1771348);


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 1.0 && Cuts::pT > 0.2*GeV);
      declare(fs, "fs");

      FastJets jetfs(fs, FastJets::ANTIKT, 0.6);
      declare(jetfs, "jets");

      // Book histograms
      book(_p["_p_Table1"], 1, 1, 1); //! figure 2 <dNch/detadphi> - toward region pT > 0.2
      book(_p["_p_Table2"], 1, 1, 2); //! figure 2 <dNch/detadphi> - away region pT > 0.2
      book(_p["_p_Table3"], 1, 1, 3); //! figure 2 <dNch/detadphi> - transverse region pT > 0.2

      book(_p["_p_Table5"], 2, 1, 2); //! figure 2 <dNch/detadphi> - transverse, pT > 0.5
      
      book(_p["_p_Table6"], 3, 1, 1); //! figure 3 <pTch> - toward region pT > 0.2
      book(_p["_p_Table7"], 3, 1, 2); //! figure 3 <pTch> - away region pT > 0.2
      book(_p["_p_Table8"], 3, 1, 3); //! figure 3 <pTch> - transverse region pT > 0.2

      book(_p["_p_Table10"], 4, 1, 2); //! figure 4 <pTch> - transverse, pT > 0.5
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 5*GeV &&
							  Cuts::abseta < 0.4 &&
							  Cuts::pT < 45*GeV);
      if(jets.size() == 0)
	vetoEvent;

      //! get the final state particles and start clustering jets -
      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
      PseudoJets ch_parts;
      for(const Particle& p : fsParticles){
	if(p.charge() != 0)
	  ch_parts.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
      }

      double table1_counts = 0;
      double table2_counts = 0;
      double table3_counts = 0;
      double table5_counts = 0;
      for (const PseudoJet& ch : ch_parts){
	if(fabs(ch.delta_phi_to(jets.at(0))) < M_PI/3){ //! toward - delta phi < pi/3
	  if(ch.pt() > 0.2){
	    table1_counts+=1.0;
	    _p["_p_Table6"]->fill(jets.at(0).pt(), ch.pt());
	  }
	}
	if(M_PI - fabs(ch.delta_phi_to(jets.at(0)))< M_PI/3){ //! away - |delta phi - pi| < pi/3 
	  if(ch.pt() > 0.2){
	    table2_counts+=1.0;
	    _p["_p_Table7"]->fill(jets.at(0).pt(), ch.pt());
	  }	  
	}
	if(fabs(ch.delta_phi_to(jets.at(0))) > M_PI/3 &&
	   fabs(ch.delta_phi_to(jets.at(0))) < 2*M_PI/3){ //! transverse - pi/3 < |delta phi| < 2pi/3
	  if(ch.pt() > 0.2){
	    table3_counts+=1.0;
	    _p["_p_Table8"]->fill(jets.at(0).pt(), ch.pt());
	  }	  
	  if(ch.pt() > 0.5){
	    table5_counts+=1.0;
	    _p["_p_Table10"]->fill(jets.at(0).pt(), ch.pt());
	  }	  
	}	
      }

      table1_counts = (double)table1_counts / (2 * 2 * M_PI/3);
      table2_counts = (double)table2_counts / (2 * 2 * M_PI/3);
      table3_counts = (double)table3_counts / (2 * 2 * M_PI/3);
      table5_counts = (double)table5_counts / (2 * 2 * M_PI/3);
      
      _p["_p_Table1"]->fill(jets.at(0).pt(), table1_counts);
      _p["_p_Table2"]->fill(jets.at(0).pt(), table2_counts);
      _p["_p_Table3"]->fill(jets.at(0).pt(), table3_counts);
      _p["_p_Table5"]->fill(jets.at(0).pt(), table5_counts);
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }


    /// @name Histograms
    //@{
    map<string, Profile1DPtr> _p;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2019_I1771348);

}
