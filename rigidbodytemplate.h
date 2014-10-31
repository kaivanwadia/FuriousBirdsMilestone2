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
    double getSignedDistance(double i, double j, double k) const;

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
    double computeMapIndex(double i, double j, double k) const;

    Mesh *m_;

    double volume_;
    Eigen::Matrix3d inertiaTensor_;
    Eigen::Vector3d principAxes_[3];

    Eigen::VectorXd signedDistMap_;
    double noOfCubes_;
    double cubeSideLength_;
    double subCubeSideLength_;
};

#endif // RIGIDBODYTEMPLATE_H
