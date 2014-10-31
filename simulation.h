#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "simparameters.h"
#include <QGLWidget>

class RigidBodyTemplate;
class RigidBodyInstance;

typedef Eigen::Triplet<double> Tr;

class SimParameters;

struct Plane
{
    Eigen::Vector3d pos;
    Eigen::Vector3d normal;
};

struct Impulse
{
    Impulse(int body1, int body2, Eigen::Vector3d impactPoint, Eigen::Vector3d contactNormal, Eigen::Vector3d tangent1, Eigen::Vector3d tangent2, double impulseMagnitude)
    {
        this->body1 = body1;
        this->body2 = body2;
        this->impactPoint = impactPoint;
        this->contactNormal = contactNormal;
        this->tangent1 = tangent1;
        this->tangent2 = tangent2;
        this->impulseMagnitude = impulseMagnitude;
        this->planeContact = false;
    }

    Impulse(int body1, Eigen::Vector3d impactPoint, Eigen::Vector3d contactNormal, Eigen::Vector3d tangent1, Eigen::Vector3d tangent2, double impulseMagnitude)
    {
        this->body1 = body1;
        this->body2 = NULL;
        this->impactPoint = impactPoint;
        this->contactNormal = contactNormal;
        this->tangent1 = tangent1;
        this->tangent2 = tangent2;
        this->impulseMagnitude = impulseMagnitude;
        this->planeContact = true;
    }

    int body1;
    int body2;
    Eigen::Vector3d impactPoint;
    Eigen::Vector3d contactNormal;
    Eigen::Vector3d tangent1;
    Eigen::Vector3d tangent2;
    double impulseMagnitude;
    bool planeContact;
};

class Simulation
{
public:
    Simulation(const SimParameters &params);
    ~Simulation();

    void takeSimulationStep();
    void initializeGL();

    void renderPlanes(bool transparent);
    void renderObjects();
    void clearScene();
    void addRigidBody(Eigen::Vector3d pos, Eigen::Vector3d lookdir);
    double computeSignedDistancePointToBody(Eigen::Vector3d point, RigidBodyInstance body);
    Eigen::Vector3d signedDistanceGrad(Eigen::Vector3d point, RigidBodyInstance body);
    void computeFrictionForces();

private:
    void loadFloorTexture();
    void loadWallTexture();
    void loadRigidBodies();

    void renderPlane(const Plane &p, bool isFloor);

    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);

    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;
    GLuint wallTex_;

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    std::vector<Plane> planes_;
    std::vector<Impulse *> impulses_;
};

#endif // SIMULATION_H
