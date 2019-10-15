#include "DetectorConstruction.hh"
#include "PMTSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"

#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include <math.h>

///////////////////////////////////////////////////////////////////////////////

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
	DefineMaterials();
}

///////////////////////////////////////////////////////////////////////////////

DetectorConstruction::~DetectorConstruction()
{ }

///////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::DefineMaterials() {

	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	G4double density;
	G4double temp;
	G4int ncomponents;

	// Water
	Water = nist->FindOrBuildMaterial("G4_WATER");

	// Seawater

	density = 1.04*g / cm3;
	temp = 283.15*kelvin;
	Seawater = new G4Material("SeaWater", density, ncomponents = 4, kStateLiquid, temp);
	G4Element* Cl = nist->FindOrBuildElement("Cl");
	G4Element* Na = nist->FindOrBuildElement("Na");
	G4Element* Mg = nist->FindOrBuildElement("Mg");
	Seawater->AddMaterial(Water, 96.88 * perCent); //fractional mass
	Seawater->AddElement(Cl, 1.92 * perCent);
	Seawater->AddElement(Na, 1.07 * perCent);
	Seawater->AddElement(Mg, 0.13 * perCent);

	// Material properties tables
	// Seawater Data

	G4double photonEnergy[] =
	{ 2.666*eV, 3.061*eV, 3.396*eV};
	const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

	G4double refractiveIndex[] =
	{ 1.333, 1.333, 1.333};
	assert(sizeof(refractiveIndex) == sizeof(photonEnergy));

	G4double absorption[] =
	{ 55.0*m, 35.0*m, 30.5*m};
	assert(sizeof(absorption) == sizeof(photonEnergy));

	G4double mie[] =
	{ 32.3*m, 19.0*m, 10.6*m };
	assert(sizeof(mie) == sizeof(photonEnergy));

	G4double MIE_water_const[3] = { 0.9204, 0.1491, 0.8831 };

	// Mar-19 STRAW data by Matthew Man


	G4MaterialPropertiesTable* MPT_Seawater = new G4MaterialPropertiesTable();

	MPT_Seawater->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries)->SetSpline(true);
	MPT_Seawater->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true);


	// sea water only MieHG enabled, comment out the following four lines for pure water simulaiton
	MPT_Seawater->AddProperty("MIEHG", photonEnergy, mie, nEntries)->SetSpline(true);
	MPT_Seawater->AddConstProperty("MIEHG_FORWARD",MIE_water_const[0]);
	MPT_Seawater->AddConstProperty("MIEHG_BACKWARD",MIE_water_const[1]);
	MPT_Seawater->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_water_const[2]);


	Seawater->SetMaterialPropertiesTable(MPT_Seawater);
}

///////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Parameters
  // Detector size
  G4double world_radius = 100 * m;
  // PMT size, PMT number can be adjusted in DetectorConstruction.hh
  G4double det_radius = 355 * mm / 2;

  // Check volumes overlaps
  G4bool checkOverlaps = true;

  // World

  G4Orb* solidWorld = new G4Orb("World", world_radius);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Seawater, "WorldLV");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,
	  "World", 0, false, 0, checkOverlaps);


  // sensitive detector

  for (int i = 0; i < 21; i++) {
	  for (int j = 0; j < 21; j++) {

		  G4Orb* solidDet = new G4Orb("Det", det_radius);

		  flogicDet =
			  new G4LogicalVolume(solidDet, Seawater, "DetLV");

		  G4double phi = 2 * std::acos(-1) * j / 20;
		  G4double theta = 2 * std::acos(-1) * (i - 0.5) / 20;

		  G4double z = sqrt(2969) * std::cos(theta) * m;
		  G4double x = sqrt(2969) * std::sin(theta) * std::cos(phi) * m;
		  G4double y = sqrt(2969) * std::sin(theta) * std::sin(phi) * m;

		  new G4PVPlacement(0, G4ThreeVector(x, y, z), flogicDet, "Det", logicWorld, false, 0, checkOverlaps);

	  }
  }

  // visualization attributes

  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  flogicDet->SetVisAttributes(visAttributes);

  return physWorld;

}

///////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::ConstructSDandField()
{
	//
	// Sensitive detectors
	//
	auto SDManager = G4SDManager::GetSDMpointer();
	G4String SDname;

	auto PMTSD1 = new PMTSD(SDname = "/PMT");
	SDManager->AddNewDetector(PMTSD1);
	flogicDet->SetSensitiveDetector(PMTSD1);

	// Add optical photon filter
	G4SDParticleFilter* OpPhotonFilter =
		new G4SDParticleFilter("OpPhotonFilter", "opticalphoton");
	PMTSD1->SetFilter(OpPhotonFilter);

}
///////////////////////////////////////////////////////////////////////////////
