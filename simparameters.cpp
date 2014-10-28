#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;
    timeStep = 0.001;
    NewtonMaxIters = 2000;
    NewtonTolerance = 1e-8;

    activeForces = F_GRAVITY;
    gravityG = -9.8;
    penaltyStiffness = 1000000;

    bodyDensity = 1.0;
    launchBody = R_SPHERE;
    launchVel = 10;
    randomLaunchOrientation = false;
    randomLaunchAngVel = false;
    randomLaunchVelMagnitude = 5.0;

    coeffRestitution = 0.1;
    coeffFriction = 1;
}
