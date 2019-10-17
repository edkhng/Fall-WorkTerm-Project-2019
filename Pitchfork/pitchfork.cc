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
	{"WorldRadius", 10},
	{"DetectorRadius", 0.5},
	{"SourceZ", 0},
	{"OpeningAngle", 180},
	{"RefractiveIndex", 1.3},
	{"AbsorptionLength", 50},
	{"ScatteringLength", 50},
	{"n", 0},
	{"ui", 1}
};

class PitchforkDetectorConstruction : public G4VUserDetectorConstruction{
public:
	PitchforkDetectorConstruction(): G4VUserDetectorConstruction(){}
    virtual G4VPhysicalVolume* Construct(){  
		G4NistManager* nist = G4NistManager::Instance();

		G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
		const size_t entries = 2;
		G4double PhotonEnergy [entries] = {1*eV, 10*eV};
		G4double RefractiveIndex [entries] = {var["RefractiveIndex"], var["RefractiveIndex"]};
		G4double AbsorptionLength [entries] = {var["AbsorptionLength"]*m, var["AbsorptionLength"]*m};
		G4double ScatteringLength [entries] = {var["ScatteringLength"]*m, var["ScatteringLength"]*m};
		G4MaterialPropertiesTable* wp = new G4MaterialPropertiesTable();
		wp->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, entries);
		wp->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLength, entries);
		wp->AddProperty("RAYLEIGH", PhotonEnergy, ScatteringLength, entries);
		water->SetMaterialPropertiesTable(wp);

		G4Orb* world = new G4Orb("World", var["WorldRadius"]*m);
		G4LogicalVolume* world_logic = new G4LogicalVolume(world, water, "World");
		G4VPhysicalVolume* world_phys = new G4PVPlacement(0, G4ThreeVector(), world_logic, 
				"World", 0, false, 0, true);

		for (int i=-2; i<3; i++){
			G4String name = "d_" + to_string(i);
			G4Orb* d = new G4Orb(name, var["DetectorRadius"]*m);
			G4LogicalVolume* d_logic = new G4LogicalVolume(d, water, name);
			new G4PVPlacement(0, G4ThreeVector(5*m,2*i*m,0), d_logic, name, world_logic, false, 0, true);
		}
		return world_phys;
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
		energy->SetMonoEnergy(3.06*eV);    // 405 nm
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
