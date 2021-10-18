// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2021_I1853218 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2021_I1853218);

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 1.0);
      declare(fs, "fs");

      //! parameters for splitting function
      z_cut=0.1;
      beta = 0;      
      
      // Book histograms
      //! Figure-4  R = 0.4, M for various pT bins 
      book(_h["_h_Table1"], 3, 1, 1);
      book(_h["_h_Table2"], 4, 1, 1);
      book(_h["_h_Table3"], 5, 1, 1);
      //Mg
      book(_h["_h_Table4"], 6, 1, 1);
      book(_h["_h_Table5"], 7, 1, 1);
      book(_h["_h_Table6"], 8, 1, 1);

      //! Figure-5  M for various radii  
      book(_h["_h_Table7"], 9, 1, 1);
      book(_h["_h_Table8"], 11, 1, 1);
      //Mg
      book(_h["_h_Table9"], 12, 1, 1);
      book(_h["_h_Table10"], 14, 1, 1);

      
    }//! init function 


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      //! get the final state particles and start clustering jets -
      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

      //! necessary jet selectors 
      fastjet::Selector select_eta_R02  = fastjet::SelectorAbsEtaMax(0.8);
      fastjet::Selector select_eta_R04  = fastjet::SelectorAbsEtaMax(0.6);
      fastjet::Selector select_eta_R06  = fastjet::SelectorAbsEtaMax(0.4);
      fastjet::Selector select_pt   = fastjet::SelectorPtMin(5.0);
      fastjet::Selector select_both_R02 = select_pt && select_eta_R02;
      fastjet::Selector select_both_R04 = select_pt && select_eta_R04;
      fastjet::Selector select_both_R06 = select_pt && select_eta_R06;

      fastjet::JetDefinition jet_def_R02(fastjet::antikt_algorithm, 0.2);
      fastjet::JetDefinition jet_def_R04(fastjet::antikt_algorithm, 0.4);
      fastjet::JetDefinition jet_def_R06(fastjet::antikt_algorithm, 0.6);
      
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd_R02(beta, z_cut, symmetry_measure, 0.2);
      fastjet::contrib::SoftDrop sd_R04(beta, z_cut, symmetry_measure, 0.4);
      fastjet::contrib::SoftDrop sd_R06(beta, z_cut, symmetry_measure, 0.6);
      
      PseudoJets parts;
      for(const Particle& p : fsParticles){
	parts.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
      }

      
       //! Analyze R = 0.2 Jets 
      fastjet::ClusterSequence cs_R02(parts, jet_def_R02);
      PseudoJets jets_R02 = select_both_R02(cs_R02.inclusive_jets(20.0));
      if(jets_R02.size()!=0){
	for(const PseudoJet& jet : jets_R02){
	  if(jet.pt() > 30 && jet.pt() < 40){
	    _h["_h_Table7"]->fill(jet.m());	      
	  }
	  PseudoJet sd_jet = sd_R02(jet);	  
	  if(sd_jet != 0){	    
	    if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table9"]->fill(sd_jet.m());	      
	    }
	  }	  
	}       
      }//! 0.2 jets
      
      //! Analyze R = 0.4 Jets 
      fastjet::ClusterSequence cs_R04(parts, jet_def_R04);
      PseudoJets jets_R04 = select_both_R04(cs_R04.inclusive_jets(20.0));
      if(jets_R04.size()!=0){
	for(const PseudoJet& jet : jets_R04){
	  if(jet.pt() > 20 && jet.pt() < 25){
	    _h["_h_Table1"]->fill(jet.m());
	  }else if(jet.pt() > 25 && jet.pt() < 30){
	    _h["_h_Table2"]->fill(jet.m());	      
	  }else if(jet.pt() > 30 && jet.pt() < 40){
	    _h["_h_Table3"]->fill(jet.m());
	    }
	  PseudoJet sd_jet = sd_R04(jet);	  
	  if(sd_jet != 0){	    
	    if(jet.pt() > 20 && jet.pt() < 25){
	      _h["_h_Table4"]->fill(sd_jet.m());
	    }else if(jet.pt() > 25 && jet.pt() < 30){
	      _h["_h_Table5"]->fill(sd_jet.m());	      
	    }else if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table6"]->fill(sd_jet.m());
	    }
	    }	  
	}       
      }//! 0.4 jets 
      
      //! Analyze R = 0.6 Jets 
      fastjet::ClusterSequence cs_R06(parts, jet_def_R06);
      PseudoJets jets_R06 = select_both_R06(cs_R06.inclusive_jets(20.0));
      if(jets_R06.size()!=0){
	for(const PseudoJet& jet : jets_R06){
	  if(jet.pt() > 30 && jet.pt() < 40){
	    _h["_h_Table8"]->fill(jet.m());	      
	  }
	  PseudoJet sd_jet = sd_R06(jet);	  
	  if(sd_jet != 0){	    
	    if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table10"]->fill(sd_jet.m());	      
	    }
	  }	  
	}       
      }//! 0.6 jets 
      
    }//! analyze function 


    /// Normalise histograms etc., after the run
    void finalize() {

      //! normalize all histograms 
      normalize(_h["_h_Table1"]);
      normalize(_h["_h_Table2"]);
      normalize(_h["_h_Table3"]);

      normalize(_h["_h_Table4"]);      
      normalize(_h["_h_Table5"]);      
      normalize(_h["_h_Table6"]);
      
      normalize(_h["_h_Table7"]);
      normalize(_h["_h_Table8"]);      

      normalize(_h["_h_Table9"]);		    
      normalize(_h["_h_Table10"]);
      
    }//! finalize functiomn 

    double z_cut;
    double beta;
    /// @name Histograms
    map<string, Histo1DPtr> _h;


  };


  DECLARE_RIVET_PLUGIN(STAR_2021_I1853218);

}
