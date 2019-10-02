#ifndef EventAction_h
#define EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#include "PMTHit.hh"
#include <vector>
#include <array>

class G4Event;

class EventAction : public G4UserEventAction
{
public:
	EventAction();
	virtual ~EventAction();

	virtual void  BeginOfEventAction(const G4Event* event);
	virtual void  EndOfEventAction(const G4Event* event);

	void IncPhotonCount_Ceren() { fPhotonCount_Ceren++; }
	//void IncAbsorption() { fAbsorptionCount++; }

	G4int GetPhotonCount_Ceren()const { return fPhotonCount_Ceren; }
	//G4int GetAbsorptionCount()const { return fAbsorptionCount; }

private:
	PMTHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;    
	G4int  fPMTHCID;

	G4int fPhotonCount_Ceren;
	//G4int fAbsorptionCount;

};

#endif
