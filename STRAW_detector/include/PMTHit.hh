#ifndef PMTHit_h
#define PMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

class G4AttDef;
class G4AttValue;

// PMT hit
// - particle time


class PMTHit : public G4VHit
{
  public:
	PMTHit();
	PMTHit(const PMTHit &right);
    virtual ~PMTHit();

	// operators
	const PMTHit& operator=(const PMTHit &right);
	G4bool operator==(const PMTHit &right) const;

	inline void* operator new(size_t);
	inline void  operator delete(void* aHit);

    // methods from base class
	// virtual void Draw();
  //   virtual void Print();
	virtual const std::map<G4String, G4AttDef>* GetAttDefs() const;
	virtual std::vector<G4AttValue>* CreateAttValues() const;

    // methods to handle data

	void SetTime(G4double t) { fTime = t; }
	G4double GetTime() const { return fTime; }


  private:
    G4double fTime;
};

using PMTHitsCollection = G4THitsCollection<PMTHit>;
extern G4ThreadLocal G4Allocator<PMTHit>* PMTHitAllocator;

inline void* PMTHit::operator new(size_t)
{
  if (!PMTHitAllocator) {
	  PMTHitAllocator = new G4Allocator<PMTHit>;
  }
  return (void *) PMTHitAllocator->MallocSingle();
}

inline void PMTHit::operator delete(void *aHit)
{
  if (!PMTHitAllocator) {
	  PMTHitAllocator = new G4Allocator<PMTHit>;
  }
  PMTHitAllocator->FreeSingle((PMTHit*) aHit);
}

#endif
