#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>

class MySteppingAction : public G4UserSteppingAction {
public:
    MySteppingAction() {
        outFile.open("results.csv");
        outFile << "Energy_keV,Angle_deg\n";
    }
    ~MySteppingAction() {
        outFile.close();
    }

    void UserSteppingAction(const G4Step* step) override {
        
        auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

        if (volume->GetName() == "logicTarget") {
            G4double edep = step->GetTotalEnergyDeposit();

            if (edep > 0) {
                auto preDir  = step->GetPreStepPoint()->GetMomentumDirection();
                auto postDir = step->GetPostStepPoint()->GetMomentumDirection();

                G4double angle = preDir.angle(postDir);

                outFile << edep/keV << "," << angle/deg << "\n";

                outFile.flush();
            }
        }
    }

private:
    std::ofstream outFile;
};


#endif

