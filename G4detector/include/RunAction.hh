#ifndef B5RunAction_h
#define B5RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class EventAction;
class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
};


#endif
