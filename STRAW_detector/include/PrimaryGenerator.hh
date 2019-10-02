#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "G4VPrimaryGenerator.hh"

class G4Event;

///////////////////////////////////////////////////////////////////////////////

class PrimaryGenerator : public G4VPrimaryGenerator
{
  public:
    PrimaryGenerator();    
   ~PrimaryGenerator();

  public:
	void SetParticlePosition(G4ThreeVector position) {fpos = position;}
	void SetParticleEnergy(G4double energy) { fenergy = energy; }
	void SetParticleType(G4int type) { ftype = type; }
    virtual void GeneratePrimaryVertex(G4Event*);
	

  private:
	G4ThreeVector fpos;
	G4double fenergy;
	G4int ftype;
};

///////////////////////////////////////////////////////////////////////////////

#endif
