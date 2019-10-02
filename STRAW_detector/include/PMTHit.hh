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
// - layer ID
// - column ID
// - CellID
// - particle time
// - particle energy
// - world position

class PMTHit : public G4VHit
{
  public:
	PMTHit();
	PMTHit(G4int layerID, G4int columnID, G4int cellID);
	PMTHit(const PMTHit &right);
    virtual ~PMTHit();

	// operators
	const PMTHit& operator=(const PMTHit &right);
	G4bool operator==(const PMTHit &right) const;

	inline void* operator new(size_t);
	inline void  operator delete(void* aHit);

    // methods from base class
	virtual void Draw();
    virtual void Print();
	virtual const std::map<G4String, G4AttDef>* GetAttDefs() const;
	virtual std::vector<G4AttValue>* CreateAttValues() const;

    // methods to handle data
	void SetLayerID(G4int z) { fLayerID = z; }
	G4int GetLayerID() const { return fLayerID; }

	void SetColumnID(G4int z) { fColumnID = z; }
	G4int GetColumnID() const { return fColumnID; }

	void SetCellID(G4int z) { fCellID = z; }
	G4int GetCellID() const { return fCellID; }

	void SetTime(G4double t) { fTime = t; }
	G4double GetTime() const { return fTime; }

	void SetEnergy(G4double e) { fEnergy = e; }
	G4double GetEnergy() const { return fEnergy; }

	void SetWorldPos(G4ThreeVector t) { fWorldPos = t; }
	G4ThreeVector GetWorldPos() const { return fWorldPos; }

      
  private:
	G4int fLayerID;
	G4int fColumnID;
	G4int fCellID;
    G4double fTime;
	G4double fEnergy;
	G4ThreeVector fWorldPos;
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