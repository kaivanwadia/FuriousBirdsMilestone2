﻿#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include "mesh.h"
#include "quadprog/eiquadprog.hpp"
const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0), wallTex_(0)
{
    gamewon = false;
    loadRigidBodies();
}

Simulation::~Simulation()
{
    clearScene();
    for(vector<RigidBodyTemplate *>::iterator it = templates_.begin(); it != templates_.end(); ++it)
    {
        delete *it;
    }
}

void Simulation::initializeGL()
{
    loadFloorTexture();
    loadWallTexture();
}

void Simulation::loadRigidBodies()
{
    const int numobjs = 5;
    string objNames[numobjs] = {"resources/sphere.obj", "resources/2by4.obj", "resources/bunny.obj", "resources/custom.obj", "resources/pyramid.obj"};
    for(int i=0; i<numobjs; i++)
    {
        RigidBodyTemplate *rbt = new RigidBodyTemplate(objNames[i]);
        templates_.push_back(rbt);
    }
}

void Simulation::loadFloorTexture()
{
    floorTex_ = SOIL_load_OGL_texture("resources/grid.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(floorTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}

void Simulation::loadWallTexture()
{
    wallTex_ = SOIL_load_OGL_texture("resources/wall.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(wallTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}

void Simulation::renderPlanes(bool transparent)
{
    renderLock_.lock();

    glEnable(GL_CULL_FACE);
    if(transparent)
    {
        glCullFace(GL_FRONT);
        glColor4f(1.0, 1.0, 1.0, 0.5);
    }
    else
    {
        glCullFace(GL_BACK);
        glColor4f(1.0, 1.0, 1.0, 1.0);
    }

    for(vector<Plane>::iterator it = planes_.begin(); it != planes_.end(); ++it)
        renderPlane(*it, it == planes_.begin());

    glDisable(GL_CULL_FACE);

    renderLock_.unlock();
}

void Simulation::renderPlane(const Plane &p, bool isFloor)
{    
    if(isFloor && floorTex_)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else if(!isFloor && wallTex_)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else
        glColor3f(0.5, 0.5, 0.5);

    double texsize = 5.0;
    double gridsize = 1000.0;

    double texmax = gridsize/texsize;

    Vector3d tangent1 = VectorMath::perpToAxis(p.normal);
    Vector3d tangent2 = tangent1.cross(p.normal);

    Vector3d corner;

    glBegin(GL_QUADS);
    {
        glTexCoord2f(texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);
    }
    glDisable(GL_TEXTURE_2D);
    glEnd();
}

void Simulation::renderObjects()
{
    renderLock_.lock();
    {
        for(vector<RigidBodyInstance *>::iterator it = bodies_.begin(); it != bodies_.end(); ++it)
        {
            (*it)->render();
        }
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    time_ += params_.timeStep;
    if(params_.gameRunning){

        if(gamewon)
        {
            cout<<"You WON !!!"<<endl;
            setupGameMode();
        }
        for(vector<RigidBodyInstance *>::iterator it = bodies_.end()-1; it != bodies_.begin(); --it)
        {
            if(((*it)->c).norm()>100){
                RigidBodyInstance *b = *it;
                bodies_.erase(it);
                delete b;
            }
        }
    }

    vector<Matrix3d> Roldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c += params_.timeStep*body.cvel;
        Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*body.w);
        Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);
        Roldthetas.push_back(Rtheta);
        body.theta = VectorMath::axisAngle(Rhw*Rtheta);

    }

    VectorXd cForce;
    VectorXd thetaForce;
    computeForces(cForce, thetaForce);

    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_.timeStep*cForce.segment<3>(3*bodyidx)/body.density/body.getTemplate().getVolume();

        Vector3d newwguess(body.w);
        Matrix3d &Roldtheta = Roldthetas[bodyidx];

        int iter = 0;
        for(iter=0; iter<params_.NewtonMaxIters; iter++)
        {
            Matrix3d Dw1 = -VectorMath::TMatrix(params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(-body.theta);
            Matrix3d Dw2 = VectorMath::TMatrix(-params_.timeStep*body.w).inverse() * VectorMath::TMatrix(-body.theta);
            Matrix3d DRw = VectorMath::DrotVector(-body.theta, newwguess);
            Matrix3d Rnewtheta = VectorMath::rotationMatrix(body.theta);
            Vector3d fval = body.density * Dw1.transpose()*Rnewtheta*Mi*Rnewtheta.transpose()*newwguess;
            fval += -params_.timeStep*body.density * DRw.transpose() * Mi * Rnewtheta.transpose() * newwguess;
            fval += body.density * Dw2.transpose() * Roldtheta * Mi * Roldtheta.transpose() * body.w;
            fval += params_.timeStep*thetaForce.segment<3>(3*bodyidx);

            if(fval.norm() < params_.NewtonTolerance)
                break;

            Matrix3d Df = body.density * Dw1.transpose()*Rnewtheta*Mi*Rnewtheta.transpose();// + -params_.timeStep*(*it)->density * DRw.transpose() * Mi * Rnewtheta.transpose();

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        body.w = newwguess;
        if(impulses_.size()>0 && params_.activeForces&SimParameters::F_FRICTION)
        {

            computeFrictionForces();
        }
    }
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        for(vector<RigidBodyInstance *>::iterator it = bodies_.begin(); it != bodies_.end(); ++it)
            delete *it;
        bodies_.clear();
        planes_.clear();
        impulses_.clear();
        Plane groundPlane;
        groundPlane.pos << 0,0,0;
        groundPlane.normal << 0,0,1;
        planes_.push_back(groundPlane);
    }
    renderLock_.unlock();
}

void Simulation::addRigidBody(Vector3d pos, Vector3d lookdir)
{
    renderLock_.lock();
    {
        if (params_.launchBody == SimParameters::R_PLANE)
        {
            Vector3d base = pos + 10*lookdir;
            Vector3d normal = -1 * lookdir;
            Plane plane;
            plane.pos = base;
            plane.normal = normal;
            planes_.push_back(plane);
        }
        else
        {
            Vector3d orient(0,0,0);
            if(params_.randomLaunchOrientation)
            {
                for(int i=0; i<3; i++)
                    orient[i] = 2.0*VectorMath::randomUnitIntervalReal()-1.0;
            }

            RigidBodyInstance *newbody = new RigidBodyInstance(*templates_[params_.launchBody], pos, orient, params_.bodyDensity);
            newbody->cvel = lookdir;
            newbody->cvel.normalize();
            newbody->cvel *= params_.launchVel;
            Vector3d orientvel(0,0,0);
            if(params_.randomLaunchAngVel)
            {
                for(int i=0; i<3; i++)
                    orientvel[i] = (2.0*VectorMath::randomUnitIntervalReal()-1.0);
                orientvel.normalize();
                orientvel *= params_.randomLaunchVelMagnitude;
            }
            newbody->w = VectorMath::rotationMatrix(orient)*orientvel;

            bodies_.push_back(newbody);
        }
    }
    renderLock_.unlock();
}

void Simulation::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc.resize(3*bodies_.size());
    Ftheta.resize(3*bodies_.size());
    Fc.setZero();
    Ftheta.setZero();
    impulses_.clear();

    Vector3d z(0,0,1);

    for(int bodyidx=0; bodyidx<(int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d bodyRotMatrix = VectorMath::rotationMatrix(body.theta);
        if(params_.activeForces & SimParameters::F_GRAVITY)
            Fc[3*bodyidx+2] += params_.gravityG*body.density*body.getTemplate().getVolume();

        Matrix3d rot = VectorMath::rotationMatrix(body.theta);


        int noOfVerts = body.getTemplate().getMesh().getNumVerts();
        for(int i=0; i<noOfVerts; i++)
        {
            double epsilon = params_.coeffRestitution;
            const Vector3d pt = body.getTemplate().getMesh().getVert(i);
            Vector3d embpt = body.c + rot*pt;

            for (int planeID = 0; planeID < planes_.size(); planeID++)
            {
                double bodyPlaneDist = (body.c - planes_[planeID].pos).dot(planes_[planeID].normal);
                if (bodyPlaneDist > 1)
                {
                    continue;
                }
                double checkValue = (embpt - planes_[planeID].pos).dot(planes_[planeID].normal);
                if (checkValue >= 0)
                {
                    continue;
                }
                double velocity = (body.cvel + body.w.cross(VectorMath::rotationMatrix(body.theta)*pt)).dot(planes_[planeID].normal);
                if (velocity<=0)
                {
                    epsilon = 1;
                }
                else
                {
                    epsilon = params_.coeffRestitution;
                }
                double planeStiffness = params_.penaltyStiffness;
                Fc.segment<3>(3*bodyidx) += -epsilon*planeStiffness*checkValue*planes_[planeID].normal/noOfVerts;
                Ftheta.segment<3>(3*bodyidx) += -epsilon*planeStiffness*checkValue/noOfVerts * VectorMath::DrotVector(body.theta, pt).transpose()*planes_[planeID].normal;

                Vector3d tangent1 = VectorMath::perpToAxis(planes_[planeID].normal);
                Vector3d tangent2 = planes_[planeID].normal.cross(tangent1);
                tangent1 = tangent1/tangent1.norm();
                tangent2 = tangent2/tangent2.norm();
                Impulse *impulse = new Impulse(bodyidx, embpt, planes_[planeID].normal, tangent1, tangent2, params_.timeStep*(-epsilon*planeStiffness*checkValue*planes_[planeID].normal/noOfVerts).norm());
                impulses_.push_back(impulse);
            }
        }
        for (int m2BodyIdx = 0; m2BodyIdx < (int)bodies_.size(); m2BodyIdx++)
        {
            if (bodyidx == m2BodyIdx)
            {
                continue;
            }
            RigidBodyInstance &m2Body = *bodies_[m2BodyIdx];
            Matrix3d m2BodyRotMatrix = VectorMath::rotationMatrix(m2Body.theta);
            if ((m2Body.c - body.c).norm() > 2)
            {
                continue;
            }
            Vector3d forceC1(0,0,0);
            Vector3d forceT1(0,0,0);
            Vector3d forceC2(0,0,0);
            Vector3d forceT2(0,0,0);
            for (int pID = 0; pID < noOfVerts; pID++)
            {
                Vector3d templatePoint = body.getTemplate().getMesh().getVert(pID);
                Vector3d embPoint = body.c + bodyRotMatrix*templatePoint;
                double checkValueSignedDistance = computeSignedDistancePointToBody(embPoint, m2Body);
                if (checkValueSignedDistance >= 0)
                {
                    continue;
                }
                if(params_.gameRunning)
                {
                    gamewon = true;
                }
                double epsilon = params_.coeffRestitution;
                Vector3d relativeVelocity = (body.cvel + body.w.cross(bodyRotMatrix*templatePoint))
                                            - (m2Body.cvel + m2Body.w.cross(body.c + bodyRotMatrix*templatePoint - m2Body.c));
                Vector3d gradDwrtQ = signedDistanceGrad(embPoint, m2Body);
                if (relativeVelocity.dot(gradDwrtQ) < 0)
                {
                    epsilon = 1;
                }
                Vector3d gradDwrtC1 = m2BodyRotMatrix*gradDwrtQ;
                Vector3d gradDwrtT1 = VectorMath::DrotVector(body.theta, templatePoint).transpose()*m2BodyRotMatrix*gradDwrtQ;
                Vector3d gradDwrtC2 = -m2BodyRotMatrix*gradDwrtQ;
                Vector3d gradDwrtT2 = VectorMath::DrotVector(-1*m2Body.theta, bodyRotMatrix*templatePoint + body.c - m2Body.c).transpose()*gradDwrtQ;
                forceC1 = forceC1 + epsilon * checkValueSignedDistance * (-gradDwrtC1);
                forceT1 = forceT1 + epsilon * checkValueSignedDistance * (-gradDwrtT1);
                forceC2 = forceC2 + epsilon * checkValueSignedDistance * (-gradDwrtC2);
                forceT2 = forceT2 + epsilon * checkValueSignedDistance * (-gradDwrtT2);

                Vector3d tangent1 = VectorMath::perpToAxis(gradDwrtQ);
                Vector3d tangent2 = gradDwrtQ.cross(tangent1);
                tangent1 = tangent1/tangent1.norm();
                tangent2 = tangent2/tangent2.norm();
                Impulse *impulse = new Impulse(bodyidx, m2BodyIdx, embPoint, gradDwrtQ, tangent1, tangent2, params_.timeStep*((epsilon * checkValueSignedDistance * (-gradDwrtC1)).norm()));
                impulses_.push_back(impulse);
            }
            double kByVerts = params_.penaltyStiffness/noOfVerts;
            Fc.segment<3>(3*bodyidx) += forceC1*kByVerts;
            Ftheta.segment<3>(3*bodyidx) += forceT1*kByVerts;
            Fc.segment<3>(3*m2BodyIdx) += forceC2*kByVerts;
            Ftheta.segment<3>(3*m2BodyIdx) += forceT2*kByVerts;
        }
    }
}

void Simulation::computeFrictionForces()
{
    VectorXd J(bodies_.size()*6);
    VectorXd x(impulses_.size()*2);
    MatrixXd Q(6 * bodies_.size(), 6 * bodies_.size());
    VectorXd v(6 * bodies_.size());
    MatrixXd B(3 * impulses_.size(), 2* impulses_.size());
    MatrixXd P(6 * bodies_.size(), 3 * impulses_.size());
    B.setZero();
    v.setZero();
    Q.setZero();
    J.setZero();
    x.setZero();
    P.setZero();
    Matrix3d ident;
    ident.setIdentity();
    for(int i = 0; i < (int) bodies_.size(); i++)
    {
        RigidBodyInstance body = *bodies_[i];
        Q.block<3,3>(6*i, 6*i) = ident* params_.bodyDensity * body.getTemplate().getVolume();
        Q.block<3,3>(6*i+3, 6*i+3) = VectorMath::rotationMatrix(body.theta).transpose() * body.getTemplate().getInertiaTensor() * VectorMath::rotationMatrix(body.theta);
        v.segment<3>(6*i) = body.cvel;
        v.segment<3>(6*i + 3) = body.w;
    }
    for(int i = 0; i < (int) impulses_.size(); i++)
    {
        B.block<3,1>(3*i, 2*i) = impulses_[i]->tangent1;
        B.block<3,1>(3*i, 2*i+1) = impulses_[i]->tangent2;
        RigidBodyInstance r1 = *bodies_[impulses_[i]->body1];
        RigidBodyInstance r2 = *bodies_[impulses_[i]->body2];
        if(impulses_[i]->body1 == impulses_[i]->body2)
        {
            Matrix3d s = (1/params_.bodyDensity) * (VectorMath::crossProductMatrix(impulses_[i]->impactPoint)*VectorMath::rotationMatrix(-r1.theta) * r1.getTemplate().getInertiaTensor().inverse()*VectorMath::rotationMatrix(-r1.theta).inverse());
            P.block<3,3>(6*impulses_[i]->body1, 3*i) = ident * (1.0/(r1.getTemplate().getVolume()*params_.bodyDensity));
            P.block<3,3>(6*impulses_[i]->body1 + 3, 3*i) = s;
        }
        else
        {
            Matrix3d s1 = (1/params_.bodyDensity) * (VectorMath::crossProductMatrix(impulses_[i]->impactPoint)*VectorMath::rotationMatrix(-r1.theta) * r1.getTemplate().getInertiaTensor().inverse()*VectorMath::rotationMatrix(-r1.theta).inverse());
            Matrix3d s2 = (1/params_.bodyDensity) * (VectorMath::crossProductMatrix(impulses_[i]->impactPoint)*VectorMath::rotationMatrix(-r2.theta) * r2.getTemplate().getInertiaTensor().inverse()*VectorMath::rotationMatrix(-r2.theta).inverse());
            P.block<3,3>(6*impulses_[i]->body1, 3*i) = ident * (1.0/(r1.getTemplate().getVolume()*params_.bodyDensity));
            P.block<3,3>(6*impulses_[i]->body1 + 3, 3*i) = s1;
            P.block<3,3>(6*impulses_[i]->body2, 3*i) = -ident * (1.0/(r2.getTemplate().getVolume()*params_.bodyDensity));
            P.block<3,3>(6*impulses_[i]->body2 + 3, 3*i) = -s2;
        }
    }
    VectorXd ce0;
    VectorXd ci0(4 * impulses_.size());
    MatrixXd CE;
    MatrixXd CI(4*impulses_.size(), 2 * impulses_.size());
    CI.setZero();
    MatrixXd constraintMatrix(4,2);
    constraintMatrix.setZero();
    constraintMatrix.coeffRef(0,0) = 1;
    constraintMatrix.coeffRef(0,1) = 1;
    constraintMatrix.coeffRef(1,0) = 1;
    constraintMatrix.coeffRef(1,1) = -1;
    constraintMatrix.coeffRef(2,0) = -1;
    constraintMatrix.coeffRef(2,1) = 1;
    constraintMatrix.coeffRef(3,0) = -1;
    constraintMatrix.coeffRef(3,1) = -1;
    for(int i=0; i<impulses_.size(); i++)
    {
        CI.block<4,2>(4*i, 2*i) = constraintMatrix;
        VectorXd mu(4);
        mu.setOnes();
        mu*=params_.coeffFriction*impulses_[i]->impulseMagnitude;
        ci0.segment<4>(4*i) = mu;
    }

    MatrixXd G(2 * impulses_.size(), 2 * impulses_.size());
    G.setZero();
    G = 2* ( B.transpose() * P.transpose() * Q * P * B );
    MatrixXd idG(2* impulses_.size(), 2 * impulses_.size());
    idG.setIdentity();
    idG *= 1e-6;
    G+= idG;
    VectorXd g0 = v.transpose() * Q * P * B;
    Eigen::solve_quadprog(G, g0, CE, ce0, -CI.transpose(), ci0, x);
    J = B * x;
    VectorXd deltaV(6*bodies_.size());
    deltaV.setZero();
    deltaV = P * J;
    for(int i = 0; i < (int) bodies_.size(); i++)
    {
        RigidBodyInstance &body = *bodies_[i];
        body.cvel += deltaV.segment<3>(6*i);
        body.w += deltaV.segment<3>(6*i+3);
    }
}

double Simulation::computeSignedDistancePointToBody(Vector3d point, RigidBodyInstance body)
{
    Vector3d q = VectorMath::rotationMatrix(-1*body.theta) * (point - body.c);
    if (q.norm() > 1.0)
    {
        return numeric_limits<double>::max();
    }
    double subCubeLength = body.getTemplate().getSubCubeSideLength();
    double cubeLength = body.getTemplate().getCubeSideLength();
    Vector3d temp(-cubeLength/2,-cubeLength/2,-cubeLength/2);
    Vector3d dXYZ = q - temp;
    double i = floor(dXYZ[0]/subCubeLength);
    double j = floor(dXYZ[1]/subCubeLength);
    double k = floor(dXYZ[2]/subCubeLength);
    if (i == 29)
    {
        i--;
    }
    if (j == 29)
    {
        j--;
    }
    if (k == 29)
    {
        k--;
    }
    double SbyL = (dXYZ[0] - i*subCubeLength)/subCubeLength;
    double TbyL = (dXYZ[1] - j*subCubeLength)/subCubeLength;
    double UbyL = (dXYZ[2] - k*subCubeLength)/subCubeLength;

    double signedDistance = UbyL*(TbyL*( SbyL*body.getTemplate().getSignedDistance(i+1, j+1, k+1) + (1-SbyL)*body.getTemplate().getSignedDistance(i, j+1, k+1) )
                                  + (1-TbyL)*( SbyL*body.getTemplate().getSignedDistance(i+1, j, k+1) + (1-SbyL)*body.getTemplate().getSignedDistance(i, j, k+1) ))
                            + (1-UbyL)*(TbyL*( SbyL*body.getTemplate().getSignedDistance(i+1, j+1, k) + (1-SbyL)*body.getTemplate().getSignedDistance(i, j+1, k) )
                                  + (1-TbyL)*( SbyL*body.getTemplate().getSignedDistance(i+1, j, k) + (1-SbyL)*body.getTemplate().getSignedDistance(i, j, k) ));

    return signedDistance;

}

Vector3d Simulation::signedDistanceGrad(Vector3d point, RigidBodyInstance body)
{
    Vector3d q = VectorMath::rotationMatrix(-1*body.theta) * (point - body.c);
    if (q.norm() > 1.1)
    {
        return Vector3d(numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max());
    }
    double L = body.getTemplate().getSubCubeSideLength();
    double cubeLength = body.getTemplate().getCubeSideLength();
    Vector3d dXYZ = q - Vector3d(-cubeLength/2,-cubeLength/2,-cubeLength/2);
    double i = floor(dXYZ[0]/L);
    double j = floor(dXYZ[1]/L);
    double k = floor(dXYZ[2]/L);
    double S = dXYZ[0] - i*L;
    double T = dXYZ[1] - j*L;
    double U = dXYZ[2] - k*L;
    Vector3d gradQ;
    gradQ.setZero();
    gradQ[0] = (U/(L*L*L)) * (T*body.getTemplate().getSignedDistance(i+1, j+1, k+1) - T*body.getTemplate().getSignedDistance(i, j+1, k+1) + (L-T)*body.getTemplate().getSignedDistance(i+1, j, k+1) + (T-L)*body.getTemplate().getSignedDistance(i, j, k+1))
            + ((L-U)/(L*L*L)) * (T*body.getTemplate().getSignedDistance(i+1, j+1, k) - T*body.getTemplate().getSignedDistance(i, j+1, k) + (L-T)*body.getTemplate().getSignedDistance(i+1, j, k) + (T-L)*body.getTemplate().getSignedDistance(i, j, k));
    gradQ[1] = (U/(L*L*L)) * (S*body.getTemplate().getSignedDistance(i+1, j+1, k+1) + (L-S)*body.getTemplate().getSignedDistance(i, j+1, k+1) - S*body.getTemplate().getSignedDistance(i+1, j, k+1) + (S-L)*body.getTemplate().getSignedDistance(i, j, k+1))
            + ((L-U)/(L*L*L)) * (S*body.getTemplate().getSignedDistance(i+1, j+1, k) + (L-S)*body.getTemplate().getSignedDistance(i, j+1, k) - S*body.getTemplate().getSignedDistance(i+1, j, k) + (S-L)*body.getTemplate().getSignedDistance(i, j, k));
    gradQ[2] = (1/L)*((T/L)*( (S/L)*body.getTemplate().getSignedDistance(i+1, j+1, k+1) + ((L-S)/L)*body.getTemplate().getSignedDistance(i, j+1, k+1) )
                    + ((L-T)/L)*( (S/L)*body.getTemplate().getSignedDistance(i+1, j, k+1) + ((L-S)/L)*body.getTemplate().getSignedDistance(i, j, k+1) ))
              - (1/L)*((T/L)*((S/L)*body.getTemplate().getSignedDistance(i+1, j+1, k) + ((L-S)/L)*body.getTemplate().getSignedDistance(i, j+1, k) )
                    + ((L-T)/L)*( (S/L)*body.getTemplate().getSignedDistance(i+1, j, k) + ((L-S)/L)*body.getTemplate().getSignedDistance(i, j, k) ));
    return gradQ;
}

void Simulation::setupGameMode()
{
        clearScene();
        gamewon = false;
        double x = rand()%70 + 20;
        double y = rand()%70 + 20;
        double z = rand()%40 + 5;
        Vector3d pos(x, y, z);
        Vector3d orient(0,0,0);
        RigidBodyInstance *newbody = new RigidBodyInstance(*templates_[SimParameters::R_SPHERE], pos, orient, params_.bodyDensity);
        bodies_.push_back(newbody);
}
