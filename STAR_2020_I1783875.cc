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
  class STAR_2020_I1783875 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2020_I1783875);

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 1.0);
      declare(fs, "fs");

      //! parameters for splitting function
      z_cut=0.1;
      beta = 0;      
      
      // Book histograms
      //! Figure-4  R = 0.4, z{g} for various pT bins 
      book(_h["_h_Table4"], 4, 1, 1);
      book(_h["_h_Table5"], 5, 1, 1);
      book(_h["_h_Table6"], 6, 1, 1);
      book(_h["_h_Table7"], 7, 1, 1);
      book(_h["_h_Table8"], 8, 1, 1);

      //! Figure-5  R = 0.4, R_{g} for various pT bins 
      book(_h["_h_Table9"], 9, 1, 1);
      book(_h["_h_Table10"], 10, 1, 1);
      book(_h["_h_Table11"], 11, 1, 1);
      book(_h["_h_Table12"], 12, 1, 1);
      book(_h["_h_Table13"], 13, 1, 1);

      //! Figure-6  z_{g} radial scans 
      book(_h["_h_Table14"], 14, 1, 1);
      book(_h["_h_Table16"], 16, 1, 1);
      book(_h["_h_Table17"], 17, 1, 1);
      book(_h["_h_Table19"], 19, 1, 1);

      //! Figure-7  R_{g} radial scans 
      book(_h["_h_Table20"], 20, 1, 1);
      book(_h["_h_Table22"], 22, 1, 1);
      book(_h["_h_Table23"], 23, 1, 1);
      book(_h["_h_Table25"], 25, 1, 1);

      //! Figure-9 supplementary material z_{g} radial scans 
      book(_h["_h_Table26"], 26, 1, 1);
      book(_h["_h_Table27"], 27, 1, 1);
      book(_h["_h_Table28"], 28, 1, 1);
      book(_h["_h_Table29"], 29, 1, 1);
      book(_h["_h_Table30"], 30, 1, 1);
      book(_h["_h_Table31"], 31, 1, 1);

      //! Figure-10 supplementary material R_{g} radial scans 
      book(_h["_h_Table32"], 32, 1, 1);
      book(_h["_h_Table33"], 33, 1, 1);
      book(_h["_h_Table34"], 34, 1, 1);
      book(_h["_h_Table35"], 35, 1, 1);
      book(_h["_h_Table36"], 36, 1, 1);
      book(_h["_h_Table37"], 37, 1, 1);
      
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
      PseudoJets jets_R02 = select_both_R02(cs_R02.inclusive_jets(15.0));
      if(jets_R02.size()!=0){
	for(const PseudoJet& jet : jets_R02){
	  PseudoJet sd_jet = sd_R02(jet);	  
	  if(sd_jet != 0){	    
	    double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	    double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();	    
	    if(jet.pt() > 15 && jet.pt() < 20){
	      _h["_h_Table14"]->fill(z);
	      _h["_h_Table20"]->fill(r);	    
	    }else if(jet.pt() > 20 && jet.pt() < 25){
	      _h["_h_Table26"]->fill(z);	     
	      _h["_h_Table32"]->fill(r);	    
	    }else if(jet.pt() > 25 && jet.pt() < 30){
	      _h["_h_Table27"]->fill(z);
	      _h["_h_Table33"]->fill(r);	    
	    }else if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table17"]->fill(z);
	      _h["_h_Table23"]->fill(r);	    
	    }else if(jet.pt() > 40 && jet.pt() < 60){
	      _h["_h_Table28"]->fill(z);	      
	      _h["_h_Table34"]->fill(r);
	    }
	  }	  
	}       
      }//! 0.2 jets 

      //! Analyze R = 0.4 Jets 
      fastjet::ClusterSequence cs_R04(parts, jet_def_R04);
      PseudoJets jets_R04 = select_both_R04(cs_R04.inclusive_jets(15.0));
      if(jets_R04.size()!=0){
	for(const PseudoJet& jet : jets_R04){
	  PseudoJet sd_jet = sd_R04(jet);	  
	  if(sd_jet != 0){	    
	    double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	    double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();	    
	    if(jet.pt() > 15 && jet.pt() < 20){
	      _h["_h_Table4"]->fill(z);
	      _h["_h_Table9"]->fill(r);
	    }else if(jet.pt() > 20 && jet.pt() < 25){
	      _h["_h_Table5"]->fill(z);	      
	      _h["_h_Table10"]->fill(r);
	    }else if(jet.pt() > 25 && jet.pt() < 30){
	      _h["_h_Table6"]->fill(z);
	      _h["_h_Table11"]->fill(r);	    
	    }else if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table7"]->fill(z);
	      _h["_h_Table12"]->fill(r);	    
	    }else if(jet.pt() > 40 && jet.pt() < 60){
	      _h["_h_Table8"]->fill(z);
	      _h["_h_Table13"]->fill(r);
	    }
	  }	  
	}       
      }//! 0.4 jets 

      //! Analyze R = 0.6 Jets 
      fastjet::ClusterSequence cs_R06(parts, jet_def_R06);
      PseudoJets jets_R06 = select_both_R06(cs_R06.inclusive_jets(15.0));
      if(jets_R06.size()!=0){
	for(const PseudoJet& jet : jets_R06){
	  PseudoJet sd_jet = sd_R06(jet);	  
	  if(sd_jet != 0){	    
	    double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	    double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();	    
	    if(jet.pt() > 15 && jet.pt() < 20){
	      _h["_h_Table16"]->fill(z);	      
	      _h["_h_Table22"]->fill(r);	    
	    }else if(jet.pt() > 20 && jet.pt() < 25){
	      _h["_h_Table29"]->fill(z);
	      _h["_h_Table35"]->fill(r);	    
	    }else if(jet.pt() > 25 && jet.pt() < 30){
	      _h["_h_Table30"]->fill(z);
	      _h["_h_Table36"]->fill(r);	    
	    }else if(jet.pt() > 30 && jet.pt() < 40){
	      _h["_h_Table19"]->fill(z);
	      _h["_h_Table25"]->fill(r);	    
	    }else if(jet.pt() > 40 && jet.pt() < 60){
	      _h["_h_Table31"]->fill(z);
	      _h["_h_Table37"]->fill(r);
	    }
	  }	  
	}       
      }//! 0.6 jets 
      
    }//! analyze function 


    /// Normalise histograms etc., after the run
    void finalize() {

      //! normalize all histograms 
      normalize(_h["_h_Table4"]);
      normalize(_h["_h_Table5"]);
      normalize(_h["_h_Table6"]);
      normalize(_h["_h_Table7"]);
      normalize(_h["_h_Table8"]);
      
      normalize(_h["_h_Table9"]);
      normalize(_h["_h_Table10"]);
      normalize(_h["_h_Table11"]);      
      normalize(_h["_h_Table12"]);
      normalize(_h["_h_Table13"]);
		    
      normalize(_h["_h_Table14"]);
      normalize(_h["_h_Table16"]);
      normalize(_h["_h_Table17"]);
      normalize(_h["_h_Table19"]);
      
      normalize(_h["_h_Table20"]);
      normalize(_h["_h_Table22"]);
      normalize(_h["_h_Table23"]);
      normalize(_h["_h_Table25"]);
      
      normalize(_h["_h_Table26"]);		    
      normalize(_h["_h_Table27"]);
      normalize(_h["_h_Table28"]);
      normalize(_h["_h_Table29"]);
      normalize(_h["_h_Table30"]);
      normalize(_h["_h_Table31"]);
      
      normalize(_h["_h_Table32"]);
      normalize(_h["_h_Table33"]);
      normalize(_h["_h_Table34"]);
      normalize(_h["_h_Table35"]);
      normalize(_h["_h_Table36"]);
      normalize(_h["_h_Table37"]);
      
    }//! finalize functiomn 

    double z_cut;
    double beta;
    /// @name Histograms
    map<string, Histo1DPtr> _h;


  };


  DECLARE_RIVET_PLUGIN(STAR_2020_I1783875);

}
