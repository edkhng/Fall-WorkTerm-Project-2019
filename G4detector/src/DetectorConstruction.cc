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
	//Seawater Data

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
	{ 32.3*m, 19.0*m, 10.6 * m };
	assert(sizeof(mie) == sizeof(photonEnergy));

	G4double MIE_water_const[3] = { 0.9204, -0.1491, 0.8831 };

	// Mar-19 STRAW data by Matthew Man
	


	G4MaterialPropertiesTable* MPT_Seawater = new G4MaterialPropertiesTable();

	MPT_Seawater->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries)->SetSpline(true);
	MPT_Seawater->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries)->SetSpline(true);


	/* sea waterm only MieHG enabled, comment out the following four lines for pure water simulaiton*/
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
  G4double world_sizeX = 150 * m, world_sizeY = 150 * m, world_sizeZ = 150 * m;
  G4double box_sizeX = 150 * m, box_sizeY = 150 * m, box_sizeZ = 150 * m;
  G4ThreeVector box_pos = G4ThreeVector(0 * m, 0 * m, 0 * m);
  // PMT size, PMT number can be adjusted in DetectorConstruction.hh
  G4double det_radius = 355 * mm / 2;

  // Check volumes overlaps
  G4bool checkOverlaps = true;

  // World
  G4Box* solidWorld =    
    new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld, Seawater, "WorldLV"); //solid, material, name
                                 
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  // detector box  
  G4Box* solidbox = 
	new G4Box("Box", 0.5*box_sizeX, 0.5*box_sizeY, 0.5*box_sizeZ);
      
  G4LogicalVolume* logicbox =                         
    new G4LogicalVolume(solidbox, Seawater, "BoxLV"); //solid, material, name
               
  new G4PVPlacement(0,                       //rotation
					box_pos,                 //position
                    logicbox,                //its logical volume
                    "Box",                   //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  // Layer X
  G4Box* solidLayer =
	  new G4Box("Layer", 0.5*box_sizeX/kDimX, 0.5*box_sizeY, 0.5*box_sizeZ);

  G4LogicalVolume* logicLayer =
	  new G4LogicalVolume(solidLayer, Seawater, "LayerLV");
 
  new G4PVReplica(
	  "Layer",             //its name
	  logicLayer,          //its logical volume
	  logicbox,            //its mother
	  kXAxis,              //axis of replication
	  kDimX,               //number of replica
	  box_sizeX / kDimX);  // witdth of replica

  // Column Y
  G4Box* solidColumn =
	  new G4Box("Column", 0.5*box_sizeX/kDimX, 0.5*box_sizeY/kDimY, 0.5*box_sizeZ);

  G4LogicalVolume* logicColumn =
	  new G4LogicalVolume(solidColumn, Seawater, "ColumnLV");

  new G4PVReplica(
	  "Column",             //its name
	  logicColumn,          //its logical volume
	  logicLayer,            //its mother
	  kYAxis,              //axis of replication
	  kDimY,               //number of replica
	  box_sizeY / kDimY);  // witdth of replica

  // Cell Z
  G4Box* solidCell =
	  new G4Box("Cell", 0.5*box_sizeX/kDimX, 0.5*box_sizeY/kDimY, 0.5*box_sizeZ/kDimY);

  G4LogicalVolume* logicCell =
	  new G4LogicalVolume(solidCell, Seawater, "CellLV");

  new G4PVReplica(
	  "Cell",                  //its name
	  logicCell,               //its logical volume
	  logicColumn,             //its mother
	  kZAxis,                  //axis of replication
	  kDimZ,                   //number of replica
	  box_sizeZ / kDimZ);      // witdth of replica


  // sensitive detector
  G4Orb* solidDet =
	  new G4Orb("Det", det_radius);

  flogicDet =
	  new G4LogicalVolume(solidDet, Seawater, "DetLV");

  new G4PVPlacement(0,            //no rotation
	  G4ThreeVector(),            //at (0,0,0)
	  flogicDet,                  //its logical volume
	  "Det",                      //its name
	  logicCell,                  //its mother  volume
	  false,                      //no boolean operation
	  0,                          //copy number
	  checkOverlaps);             //overlaps checking



  // visualization attributes 

  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes->SetVisibility(false);
  logicCell->SetVisAttributes(visAttributes);
  logicColumn->SetVisAttributes(visAttributes);
  logicLayer->SetVisAttributes(visAttributes);
  visAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // LightGray
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
