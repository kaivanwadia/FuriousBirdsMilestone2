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

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0), wallTex_(0)
{
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
    const int numobjs = 4;
    string objNames[numobjs] = {"resources/sphere.obj", "resources/2by4.obj", "resources/bunny.obj", "resources/custom.obj"};
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

    Vector3d z(0,0,1);

    for(int bodyidx=0; bodyidx<(int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];

        if(params_.activeForces & SimParameters::F_GRAVITY)
            Fc[3*bodyidx+2] += params_.gravityG*body.density*body.getTemplate().getVolume();

        Matrix3d rot = VectorMath::rotationMatrix(body.theta);


        int nverts = body.getTemplate().getMesh().getNumVerts();
        for(int i=0; i<nverts; i++)
        {
            double epsilon = params_.coeffRestitution;
            const Vector3d pt = body.getTemplate().getMesh().getVert(i);
            Vector3d embpt = body.c + rot*pt;

            for (int planeID = 0; planeID < planes_.size(); planeID++)
            {
                double velocity = (body.cvel + body.w.cross(VectorMath::rotationMatrix(body.theta)*pt)).dot(planes_[planeID].normal);
                if (velocity<=0)
                {
                    epsilon = 1;
                }
                else
                {
                    epsilon = params_.coeffRestitution;
                }
                double bodyPlaneDist = (body.c - planes_[planeID].pos).dot(planes_[planeID].normal);
                if (bodyPlaneDist > 1)
                {
                    continue;
                }
                double checkValue = (embpt - planes_[planeID].pos).dot(planes_[planeID].normal);
                if (checkValue < 0)
                {
                    double planeStiffness = params_.penaltyStiffness;
                    Fc.segment<3>(3*bodyidx) += -epsilon*planeStiffness*checkValue*planes_[planeID].normal/nverts;
                    Ftheta.segment<3>(3*bodyidx) += -epsilon*planeStiffness*checkValue/nverts * VectorMath::DrotVector(body.theta, pt).transpose()*planes_[planeID].normal;
                }
            }
        }
    }
}

