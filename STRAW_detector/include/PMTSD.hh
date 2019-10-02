#ifndef PMTSD_h
#define PMTSD_h 1

#include "G4VSensitiveDetector.hh"

#include "PMTHit.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class PMTSD : public G4VSensitiveDetector
{
  public:
    PMTSD(G4String name);
    virtual ~PMTSD();

	// Called for each step in sensitive volume
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);

	// Called at the beginning/end of the event
	virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    PMTHitsCollection* fHitsCollection;
    G4int fHCID;
	G4bool print = false;
};

///////////////////////////////////////////////////////////////////////////////

#endif

