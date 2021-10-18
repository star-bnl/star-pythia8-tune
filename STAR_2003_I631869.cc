// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "math.h"

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2003_I631869 : public Analysis {
  public:
      
    double pimass = 0;

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2003_I631869);

    /// Book histograms and initialise projections before the run
    void init() {
        
      const UnstableParticles np(Cuts::eta > 3.4 && Cuts::eta < 4. && Cuts::pid == 111);
      declare(np, "np");
        
      book(_h["InclPI0Xsec"], 1, 1, 1);
        
    }
      
      
    //   Perform the per-event analysis
    void analyze(const Event& event)
    {

      Particles npParticles = applyProjection<UnstableParticles>(event,"np").particles();
       
      Particle lead;
        
      if(npParticles.size() > 0) lead = npParticles[0];


      for(const Particle& p : npParticles)
        {
	  double ppT = p.pT()/GeV;
	  double piEner = p.E()/GeV;
            
	  pimass = p.mass()/ GeV;
            
	  if(ppT > 1)
            {                
	      if(p.E()/GeV > lead.E()/GeV) lead = p;
                
	      _h["InclPI0Xsec"]->fill(lead.E()/GeV, piEner);

            }
        }
    }

            
    // Normalise histograms etc., after the run
    void finalize() {
      
      int nselected = 0;
      for(size_t i =0; i<_h["InclPI0Xsec"]->numBinsX(); i++)
	{
          nselected++;
          std::pair<double, double> binedges = _h["InclPI0Xsec"]->bin(i).xEdges();;

          double piP1 = pow(pow(binedges.first,2)-pow(pimass,2),0.5);
          double piP2 = pow(pow(binedges.second,2)-pow(pimass,2),0.5);
          double trigpiP = piP2 - piP1;
          
          double width = _h["InclPI0Xsec"]->bin(i).xMax() - _h["InclPI0Xsec"]->bin(i).xMin();
          _h["InclPI0Xsec"]->bin(i).scaleW(width/pow(fabs(trigpiP),3));

	}
      
      double norm = (crossSection()/microbarn)/(sumOfWeights());

      _h["InclPI0Xsec"]->scaleW(norm);
      
      
    }


    map<string, Histo1DPtr> _h;      
      
  };
  DECLARE_RIVET_PLUGIN(STAR_2003_I631869);

}
