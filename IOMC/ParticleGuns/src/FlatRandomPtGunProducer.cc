/*
 *  \author Julia Yarba
 */

#include <ostream>

#include "IOMC/ParticleGuns/interface/FlatRandomPtGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

FlatRandomPtGunProducer::FlatRandomPtGunProducer(const ParameterSet& pset) : 
   BaseFlatGunProducer(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters") ;
  
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");
  
  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

FlatRandomPtGunProducer::~FlatRandomPtGunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}


bool myIsMuonPassScint(double dVx, double dVy, double dVz, double dPx, double dPy, double dPz) {
  // To test the drop-down of efficiency at edges, we can set the cut looser
  //double ScintilXMin = -1000.0;
  //double ScintilXMax =  1000.0;
  //double ScintilZMin =  -605.6;
  //double ScintilZMax =   950.0;
/*
  double ScintilXMin = -1000.0;
  double ScintilXMax =  1000.0;
  double ScintilZMin = -1000.0;
  double ScintilZMax =  1000.0;
  
  double ScintilLowerY = -114.85;
  double ScintilUpperY = 1540.15;
  
  double dTLower = ( ScintilLowerY - dVy ) / dPy;  
  double dXLower = dVx + dTLower * dPx;
  double dZLower = dVz + dTLower * dPz;
  
  double dTUpper = ( ScintilUpperY - dVy ) / dPy;
  double dXUpper = dVx + dTUpper * dPx;
  double dZUpper = dVz + dTUpper * dPz;
  
  if (( ScintilXMin <= dXLower && dXLower <= ScintilXMax && ScintilZMin <= dZLower && dZLower <= ScintilZMax ) &&
      ( ScintilXMin <= dXUpper && dXUpper <= ScintilXMax && ScintilZMin <= dZUpper && dZUpper <= ScintilZMax ))
  {
    return true;
  }
*/
  double ScintilXMin = -1000.0;
  double ScintilXMax =  1000.0;
  double ScintilYMin = -1000.0;
  double ScintilYMax =  1000.0;
  
  double ScintilLowerZ = -114.85;
  double ScintilUpperZ = 1540.15;
  
  double dTLower = ( ScintilLowerZ - dVz ) / dPz;  
  double dXLower = dVx + dTLower * dPx;
  double dYLower = dVy + dTLower * dPy;
  
  double dTUpper = ( ScintilUpperZ - dVz ) / dPz;
  double dXUpper = dVx + dTUpper * dPx;
  double dYUpper = dVy + dTUpper * dPy;
  
  if (( ScintilXMin <= dXLower && dXLower <= ScintilXMax && ScintilYMin <= dYLower && dYLower <= ScintilYMax ) &&
      ( ScintilXMin <= dXUpper && dXUpper <= ScintilXMax && ScintilYMin <= dYUpper && dYUpper <= ScintilYMax ))
  {
    return true;
  }
  
  else return false;
}

void FlatRandomPtGunProducer::produce(Event &e, const EventSetup& es) 
{
   edm::Service<edm::RandomNumberGenerator> rng;
   CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

   if ( fVerbosity > 0 )
   {
      cout << " FlatRandomPtGunProducer : Begin New Event Generation" << endl ; 
   }
   // event loop (well, another step in it...)
          
   // no need to clean up GenEvent memory - done in HepMCProduct
   // 
   
   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
   
   // now actualy, cook up the event from PDGTable and gun parameters
   //
   // 1st, primary vertex
   //
   //HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));
   /*double dVx = CLHEP::RandFlat::shoot(engine, -1200.0, 1200.0) ;
   double dVy = -100.0;
   double dVz = CLHEP::RandFlat::shoot(engine, -650.0, 750.0) ;
   HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(dVx,dVy,dVz));*/
   double dVx;
   //double dVy = 1540.15; // same Y as the upper scintillator
   //double dVz;
   double dVy;
   double dVz = 1540.15; // same Y as the upper scintillator
   HepMC::GenVertex* Vtx = NULL;

   // loop over particles
   //
   int barcode = 1 ;
   for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
   {
       double px, py, pz, mom;
       double phi, theta;
       int j = 0;
       
       while (j < 10000) // j < 10000 to avoid too long computational time
       {

         //dVx = CLHEP::RandFlat::shoot(engine, -1000.0, 1000.0) ;
         //dVz = CLHEP::RandFlat::shoot(engine, -605.6, 950.0) ;
         dVx = CLHEP::RandFlat::shoot(engine, -1000.0, 1000.0) ;
         dVy = CLHEP::RandFlat::shoot(engine, -1000.0, 1000.0) ;
         
         mom   = CLHEP::RandFlat::shoot(engine, fMinPt, fMaxPt) ;
         phi   = CLHEP::RandFlat::shoot(engine, fMinPhi, fMaxPhi) ;
         theta = 0;

         if (fIsThetaFlat)
         {
             theta  = CLHEP::RandFlat::shoot(engine, fMinTheta, fMaxTheta);
         }
         
         if (!fIsThetaFlat)
         {
             double u = CLHEP::RandFlat::shoot(engine, 0.0, 0.785398); // u = Uniform[0;Pi/4]
             theta = 0;
             while(abs(u-(0.5*theta+0.25*sin(2*theta)))>0.000015)
             {
                 theta+=0.00001;
             }             
         }

         // mom = total momentum in this code
         px     =  mom*sin(theta)*cos(phi) ;
         //pz     =  mom*sin(theta)*sin(phi) ;
         //py     = -mom*cos(theta) ; // with the - sign, the muons are going downwards: falling from the sky
         // z-axis -> -y-azis
         // y-axis -> z-axis

         py     =  mom*sin(theta)*sin(phi) ;
         pz     =  mom*cos(theta) ;
         // If -3.141592<phi<0 or 3.141592<phi<2.141592*2 then the muons are going downwards (-z direction) : falling from the sky

         if ( myIsMuonPassScint(dVx, dVy, dVz, px, py, pz) == true ) break; // muon passing through both the scintillators => valid: the loop can be stopped
         
         j++;
         
       }

       //cout<<"j "<<j<<", mom "<<mom<<", px "<<px<<", py "<<py<<", pz "<<pz<<", theta "<<theta<<", phi "<<phi<<endl;
       
       int PartID = fPartIDs[ip] ;
       const HepPDT::ParticleData* 
          PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
       double mass   = PData->mass().value() ;
       Vtx = new HepMC::GenVertex(HepMC::FourVector(dVx,dVy,dVz));

       double energy2= mom*mom + mass*mass ;
       double energy = sqrt(energy2) ;
       HepMC::FourVector p(px,py,pz,energy) ;
       HepMC::GenParticle* Part = new HepMC::GenParticle(p,PartID,1);
       Part->suggest_barcode( barcode ) ;
       barcode++ ;
       Vtx->add_particle_out(Part);

       if ( fAddAntiParticle )
       {
          HepMC::FourVector ap(-px,-py,-pz,energy) ;
	  int APartID = -PartID ;
	  if ( PartID == 22 || PartID == 23 )
	  {
	     APartID = PartID ;
	  }	  
	  HepMC::GenParticle* APart = new HepMC::GenParticle(ap,APartID,1);
	  APart->suggest_barcode( barcode ) ;
	  barcode++ ;
	  Vtx->add_particle_out(APart) ;
       }

   }

   fEvt->add_vertex(Vtx) ;
   fEvt->set_event_number(e.id().event()) ;
   fEvt->set_signal_process_id(20) ; 
        
   if ( fVerbosity > 0 )
   {
      fEvt->print() ;  
   }

   unique_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
   BProduct->addHepMCData( fEvt );
   e.put(std::move(BProduct), "unsmeared");

   unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
   e.put(std::move(genEventInfo));
    
   if ( fVerbosity > 0 )
   {
      // for testing purpose only
      // fEvt->print() ; // prints empty info after it's made into edm::Event
      cout << " FlatRandomPtGunProducer : Event Generation Done " << endl;
   }
}
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomPtGunProducer);
