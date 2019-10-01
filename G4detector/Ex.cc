#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4TrajectoryDrawByParticleID.hh"

//#ifdef G4MULTITHREADED
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif
#include "G4Threading.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronicProcessStore.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

///////////////////////////////////////////////////////////////////////////////

namespace {
	void PrintUsage() {
		G4cerr << " Usage: " << G4endl;
		G4cerr << " Ex [-m macro ] [-e eventid] [-f filename]"
			<< G4endl;
		G4cerr << "    default filemane: default" << G4endl;
		G4cerr << "    default eventid: 0" << G4endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc,char** argv)
{
  if (argc > 7) {
	PrintUsage();
	return 1;
  }
  G4String macro;
  G4String eventid = "0";
  G4String filename = "default";
  for (G4int i = 1; i < argc; i = i + 2) {
	  if (G4String(argv[i]) == "-m") macro = argv[i + 1];
	  else if (G4String(argv[i]) == "-e") eventid = argv[i + 1];
	  else if (G4String(argv[i]) == "-f") filename = argv[i + 1];
  }

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if (macro.size() == 0) {
	  ui = new G4UIExecutive(argc, argv);
  }

  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new DetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  opticalPhysics->SetWLSTimeProfile("delta");

  // the total photons generated by the track is unchanged.
  opticalPhysics->SetMaxNumPhotonsPerStep(20); // Step size
  opticalPhysics->SetMaxBetaChangePerStep(40.0); // Allowed change in percent

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);

  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization());
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;

  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  // set initial verbosities
  UImanager->ApplyCommand("/control/verbose 0");
  UImanager->ApplyCommand("/run/verbose 0");
  UImanager->ApplyCommand("/event/verbose 0");
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/hits/verbose 0");
  UImanager->ApplyCommand("/material/verbose 0");
  UImanager->ApplyCommand("/process/setVerbose 0 all");
  UImanager->ApplyCommand("/process/verbose 0");
  UImanager->ApplyCommand("/process/eLoss/verbose 0");
  UImanager->ApplyCommand("/process/had/verbose 0");
  UImanager->ApplyCommand("/process/em/verbose 0");
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  
  
  // Process macro or start UI session
  //

  if (macro.size()) {
	  // Batch mode
	  UImanager->ApplyCommand("/run/initialize");
	  G4String command = "/Ex/eventID ";
	  UImanager->ApplyCommand(command + eventid);

	  command = "/analysis/setFileName ";
	  UImanager->ApplyCommand(command + filename + eventid);

	  command = "/control/execute ";
	  UImanager->ApplyCommand(command + macro);

  }
  else // Define UI session for interactive mode
  {
	  UImanager->ApplyCommand("/control/execute init_vis.mac");
	  ui->SessionStart();
	  delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}