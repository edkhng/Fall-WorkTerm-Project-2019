/* created 2019-10-07 by Andreas Gaertner
based on the 2015 POCAM simulation by Kai Krings for the PINGU collaboration */

#include "FTFP_BERT.hh"
#include "G4MuonMinus.hh"
#include "G4NistManager.hh"
#include "G4OpticalPhysics.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4SingleParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"
using namespace std;

map<string, double> var = {
	{"WorldRadius", 100},
	{"DetectorRadius", 0.1775},
	{"SourceZ", 0},
	{"OpeningAngle", 180},
	{"RefractiveIndex", 1.333},
	{"AbsorptionLength", 50},
	{"ScatteringLength", 50},
	{"n", 0},
	{"ui", 1}
};

class PitchforkDetectorConstruction : public G4VUserDetectorConstruction{
public:
	PitchforkDetectorConstruction(): G4VUserDetectorConstruction(){}
    virtual G4VPhysicalVolume* Construct(){
		// Get nist material manager
		G4NistManager* nist = G4NistManager::Instance();
		G4double density;
		G4double temp;
		G4int ncomponents;
		G4Material* Water = nist->FindOrBuildMaterial("G4_WATER");

		// Seawater
		density = 1.04*g / cm3;
		temp = 283.15*kelvin;
		G4Material* Seawater = new G4Material("SeaWater", density, ncomponents = 4, kStateLiquid, temp);
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

		// Parameters
	  // Detector size
	  G4double world_radius = 100 * m;
	  G4double det_radius = 355 * mm / 2;

		// Check volumes overlaps
		G4bool checkOverlaps = true;

		G4Orb* solidWorld = new G4Orb("World", world_radius);
		G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Seawater, "WorldLV");
		G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,
		  	"World", 0, false, 0, checkOverlaps);

		for (int i = 0; i < 21; i++) {
			for (int j = 0; j < 21; j++) {
				G4String name = "d_" + std::to_string(i) + "_" + std::to_string(j);
				G4Orb* solidDet = new G4Orb(name, det_radius);

				G4LogicalVolume* flogicDet =
					new G4LogicalVolume(solidDet, Seawater, name);

				G4double phi = 2 * std::acos(-1) * j / 20;
				G4double theta = 2 * std::acos(-1) * (i - 0.5) / 20;

				G4double z = sqrt(2969) * std::cos(theta) * m;
				G4double x = sqrt(2969) * std::sin(theta) * std::cos(phi) * m;
				G4double y = sqrt(2969) * std::sin(theta) * std::sin(phi) * m;

				new G4PVPlacement(0, G4ThreeVector(x, y, z), flogicDet, name, logicWorld, false, 0, checkOverlaps);

			}
		}
		return physWorld;
	}
};

class PitchforkPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
public:
	PitchforkPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(){
		source = new G4SingleParticleSource();
		source->SetParticleDefinition(G4OpticalPhoton::Definition());
		G4SPSPosDistribution* position = source->GetPosDist();
		position->SetPosDisType("Point");
		position->SetCentreCoords(G4ThreeVector(0,0,var["SourceZ"]*m));
		G4SPSAngDistribution* angle = source->GetAngDist();
		angle->SetAngDistType("iso");
		angle->DefineAngRefAxes("angref1", G4ThreeVector(1., 0., 0.));
		angle->DefineAngRefAxes("angref2", G4ThreeVector(0., 1., 0.));
		angle->SetMinTheta(0.);
		angle->SetMaxTheta(var["OpeningAngle"]*deg);
		angle->SetMinPhi(0.);
		angle->SetMaxPhi(360.*deg);
		G4SPSEneDistribution* energy = source->GetEneDist();
		energy->SetEnergyDisType("Gauss");
		energy->SetMonoEnergy(3.061*eV);    // 405 nm
		energy->SetBeamSigmaInE(7.56e-2*eV);
	}
	virtual ~PitchforkPrimaryGeneratorAction(){
		delete source;
	}
	virtual void GeneratePrimaries(G4Event* event){
    	source->SetParticlePolarization(G4RandomDirection());
    	source->SetParticleTime(0);
    	source->GeneratePrimaryVertex(event);
	}
private:
	G4SingleParticleSource* source;
};

class PitchforkSteppingAction : public G4UserSteppingAction{
public:
	PitchforkSteppingAction() : G4UserSteppingAction(){}
	virtual void UserSteppingAction(const G4Step* step){
		G4StepPoint* startpoint = step->GetPreStepPoint();
		G4TouchableHandle geometry = startpoint->GetTouchableHandle();
		G4String volumeName = geometry->GetVolume()->GetName();
		static bool first = true;
		if (volumeName.substr(0,2) == "d_"){
			G4double time = startpoint->GetGlobalTime()/ns;
			G4ThreeVector position = startpoint->GetPosition();
			G4ThreeVector direction = startpoint->GetMomentumDirection();
			if (first){ G4cout << "#volume\ttime[ns]\t theta[deg]\tphi[deg]" << G4endl; first = false; }
			G4cout << volumeName << "\t" << time << "\t" << direction.theta()/deg
					<< "\t" << direction.phi()/deg << G4endl;
			step->GetTrack()->SetTrackStatus(fStopAndKill);
		}
	}
};


int main(int argc,char** argv){
	for (int i = 1; i < argc; i += 2){
		if (var.find(string(argv[i])) == var.end()){
			throw invalid_argument("Error parsing args. Unknown parameter.");}
		try {
			var[string(argv[i])] = stod(argv[i+1]);
		}catch (const std::exception& e){
			throw invalid_argument("Error parsing args. Value missing.");}}
	G4cout << "--- Using the following variables ---" << G4endl;
	for (map<string, double>::iterator it = var.begin(); it != var.end(); it++)
		G4cout << it->first << " : " << it->second << G4endl;
	G4cout << "---" << G4endl;

	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(0);
	G4RunManager* runManager = new G4RunManager;
	runManager->SetUserInitialization(new PitchforkDetectorConstruction());
	G4VModularPhysicsList* p = new FTFP_BERT;
	G4OpticalPhysics* o = new G4OpticalPhysics();
	o->SetMaxNumPhotonsPerStep(100);
	o->SetMaxBetaChangePerStep(0.1);
	o->SetTrackSecondariesFirst(kCerenkov, true);
	p->RegisterPhysics(o);
	runManager->SetUserInitialization(p);
	runManager->SetUserAction(new PitchforkPrimaryGeneratorAction());
	runManager->SetUserAction(new PitchforkSteppingAction());

	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
	G4UImanager* u = G4UImanager::GetUIpointer();
	u->ApplyCommand("/run/initialize");

	if (var["ui"] != 0.0){
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		u->ApplyCommand("/vis/open OGL 600x600-0+0");
		u->ApplyCommand("/vis/drawVolume");
		u->ApplyCommand("/vis/viewer/set/lineSegmentsPerCircle 15");
		u->ApplyCommand("/vis/scene/add/trajectories smooth");
		u->ApplyCommand("/vis/scene/endOfEventAction accumulate");
		u->ApplyCommand("/vis/viewer/set/style wireframe");
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
		ui->SessionStart();
		delete ui;
	}else{
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
	}
	delete visManager;
	delete runManager;
}
