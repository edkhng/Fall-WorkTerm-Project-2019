#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

// number of PMTs
const G4int kDimX = 2;
const G4int kDimY = 1;
const G4int kDimZ = 3;

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

private:

	void DefineMaterials();
	G4Material* Water;
	G4Material* Seawater;
	G4LogicalVolume* flogicDet;
};


#endif

