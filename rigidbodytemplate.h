#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <map>

class Mesh;
class SignedDistanceField;

class RigidBodyTemplate
{
public:
    RigidBodyTemplate(std::string &meshFilename);
    ~RigidBodyTemplate();

    const Mesh &getMesh() const {return *m_;}
    double getVolume() const {return volume_;}
    const Eigen::Vector3d getPrincipleAxis(int axis) const {return principAxes_[axis];}

    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}
    const double getNoOfCubes() const {return noOfCubes_;}
    const double getCubeSideLength() const {return cubeSideLength_;}
    const double getSubCubeSideLength() const {return subCubeSideLength_;}
    const std::map<double, double> getSignedDistanceField() const {return signedDistMap_;}

private:
    RigidBodyTemplate(const RigidBodyTemplate &other);
    RigidBodyTemplate &operator=(const RigidBodyTemplate &other);

    void computeVolume();
    Eigen::Vector3d computeCenterOfMass();
    double boundingSphereRadius();
    void computeInertiaTensor();

    bool insideTemplate(Eigen::Vector3d point);
    double shortestDistanceFromPointToM(Eigen::Vector3d point);
    void computeSignedDistanceField(std::string meshFilename);

    void initializeSignedDistanceField(std::string fileName);

    Mesh *m_;

    double volume_;
    Eigen::Matrix3d inertiaTensor_;
    Eigen::Vector3d principAxes_[3];

    std::map<double, double> signedDistMap_;
    double noOfCubes_;
    double cubeSideLength_;
    double subCubeSideLength_;
};

#endif // RIGIDBODYTEMPLATE_H
