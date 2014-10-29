#include "rigidbodytemplate.h"
#include "mesh.h"
#include <iostream>
#include <Eigen/Dense>
#include "vectormath.h"
#include "collisions/CTCD.h"
#include "collisions/Distance.h"
#include <fstream>
#include <map>

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(std::string &meshFilename) : volume_(0)
{
    m_ = new Mesh(meshFilename);
    inertiaTensor_.setZero();
    for(int i=0; i<3; i++)
        principAxes_[i].setZero();

    computeVolume();
    Vector3d cm = computeCenterOfMass();
    m_->translate(-cm);
    double r = boundingSphereRadius();
    m_->scale(1.0/r);
    volume_ /= r*r*r;
    computeInertiaTensor();
    if (meshFilename != "resources/bunny.obj")
    {
        computeSignedDistanceField(meshFilename);
    }
}

RigidBodyTemplate::~RigidBodyTemplate()
{
    delete m_;
}

void RigidBodyTemplate::computeVolume()
{
    volume_ = 0;
    int nfaces = m_->getNumFaces();
    for(int i=0; i<nfaces; i++)
    {
        Vector3d sum(0,0,0);
        for(int j=0; j<3; j++)
            sum += m_->getVert(m_->getFace(i)[j]);
        volume_ += m_->getFaceArea(i)/9.0 * sum.dot(m_->getFaceNormal(i));
    }
}

Vector3d RigidBodyTemplate::computeCenterOfMass()
{
    Vector3d cm(0,0,0);
    int nfaces = m_->getNumFaces();
    for(int i=0; i<nfaces; i++)
    {
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = m_->getVert(m_->getFace(i)[j]);
        double area = m_->getFaceArea(i);
        Vector3d normal = m_->getFaceNormal(i);

        Vector3d term(0,0,0);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                for(int l=k; l<3; l++)
                    term[j] += pts[k][j]*pts[l][j];
            }
            term[j] *= area*normal[j]/12.0;
        }

        cm += term;
    }
    return cm/volume_;
}

double RigidBodyTemplate::boundingSphereRadius()
{
    double maxradius = 0;
    for(int i=0; i<m_->getNumVerts(); i++)
    {
        maxradius = max(maxradius, m_->getVert(i).norm());
    }
    return maxradius;
}

void RigidBodyTemplate::computeInertiaTensor()
{
    Vector3d quads(0,0,0);
    Vector3d mixed(0,0,0);
    int nfaces = m_->getNumFaces();
    for(int i=0; i<nfaces; i++)
    {
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = m_->getVert(m_->getFace(i)[j]);
        double area = m_->getFaceArea(i);
        Vector3d normal = m_->getFaceNormal(i);

        Vector3d term(0,0,0);
        Vector3d mixterm(0,0,0);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                for(int l=k; l<3; l++)
                {
                    for(int m=l; m<3; m++)
                        term[j] += pts[k][j]*pts[l][j]*pts[m][j];
                }
            }            
            term[j] *= area*normal[j]/30.0;
        }
        double mix = 0;
        for(int j=0; j<3; j++)
        {
            mix += 6.0*pts[j][0]*pts[j][1]*pts[j][2];
            for(int k=0; k<3; k++)
            {
                mix += 2.0*pts[j][k]*pts[j][(k+1)%3]*pts[(j+1)%3][(k+2)%3];
                mix += 2.0*pts[j][k]*pts[j][(k+1)%3]*pts[(j+2)%3][(k+2)%3];
            }
            mix += pts[j][0]*pts[(j+1)%3][1]*pts[(j+2)%3][2];
            mix += pts[j][2]*pts[(j+1)%3][1]*pts[(j+2)%3][0];
        }
        for(int j=0; j<3; j++)
            mixterm[j] = mix*area*normal[j]/60.0;

        quads += term;
        mixed += mixterm;
    }

    inertiaTensor_ << quads[1]+quads[2], -mixed[2], -mixed[1],
            -mixed[2], quads[0]+quads[2], -mixed[0],
            -mixed[1], -mixed[0], quads[0]+quads[1];

    SelfAdjointEigenSolver<Matrix3d> es(inertiaTensor_);
    for(int i=0; i<3; i++)
    {
        principAxes_[i] = es.eigenvectors().col(i);
        principAxes_[i].normalize();
    }
}

bool RigidBodyTemplate::insideTemplate(Vector3d point)
{
    Vector3d pointEnd;
    int evenCount = 0;
    int oddCount = 0;
    bool inside = false;
    double time = 0;
    for (int rayNo = 1; rayNo <= 10; rayNo++)
    {
        pointEnd.setZero();
        pointEnd<<VectorMath::randomUnitIntervalReal()*2 - 1, VectorMath::randomUnitIntervalReal()*2 - 1, VectorMath::randomUnitIntervalReal()*2 - 1;
        pointEnd = pointEnd * 2.1;
        inside = false;
        for (int faceNo = 0; faceNo < m_->getNumFaces(); faceNo++)
        {
            if (CTCD::vertexFaceCTCD(point, m_->getVert(m_->getFace(faceNo)[0]), m_->getVert(m_->getFace(faceNo)[1]), m_->getVert(m_->getFace(faceNo)[2]),
                                     pointEnd, m_->getVert(m_->getFace(faceNo)[0]), m_->getVert(m_->getFace(faceNo)[1]), m_->getVert(m_->getFace(faceNo)[2]), 10e-8, time))
            {
                inside = !inside;
            }
        }
        if (inside)
            oddCount++;
        else
            evenCount++;
    }
    if (evenCount > oddCount)
        return false;
    else
        return true;
}

double RigidBodyTemplate::shortestDistanceFromPointToM(Vector3d point)
{
    double minimumDistance = numeric_limits<double>::max();
    double baryIgnore1 = 0;
    double baryIgnore2 = 0;
    double baryIgnore3 = 0;
    for (int faceNo = 0; faceNo < m_->getNumFaces(); faceNo++)
    {
        Vector3d dist = Distance::vertexFaceDistance(point, m_->getVert(m_->getFace(faceNo)[0]), m_->getVert(m_->getFace(faceNo)[1]), m_->getVert(m_->getFace(faceNo)[2]),
                                                   baryIgnore1, baryIgnore2, baryIgnore3);
        if (dist.norm() < minimumDistance)
        {
            minimumDistance = dist.norm();
        }
    }
    return minimumDistance;
}

void RigidBodyTemplate::computeSignedDistanceField(string meshFileName)
{
    meshFileName = meshFileName.substr(0, meshFileName.length() - 3) + "sdf";
    ifstream inputFileStream(meshFileName.c_str());
    if (inputFileStream.good())
    {
        inputFileStream.close();
        initializeSignedDistanceField(meshFileName);
        return;
    }
    ofstream outputFileStream(meshFileName.c_str(), ofstream::out);
    outputFileStream<<"30"<<endl;
    for (int xIndex = 0; xIndex < 30; xIndex++)
    {
        for (int yIndex = 0; yIndex < 30; yIndex++)
        {
            for (int zIndex = 0; zIndex < 30; zIndex++)
            {
                Vector3d point;
                point<<(xIndex/29.0)*2-1,(yIndex/29.0)*2-1,(zIndex/29.0)*2-1;
                double distance = shortestDistanceFromPointToM(point);
                if (insideTemplate(point))
                {
                    distance = distance * -1;
                }
                outputFileStream<<distance<<" ";
//                cout<<"Distance : "<<distance;
            }
        }
    }
    outputFileStream.close();
    initializeSignedDistanceField(meshFileName);
}

void RigidBodyTemplate::initializeSignedDistanceField(string fileName)
{
    ifstream inputFileStream(fileName.c_str());
    inputFileStream >> noOfCubes_;
    int x = 0;
    int y = 0;
    int z = 0;
    while (inputFileStream)
    {
        inputFileStream >> signedDistMap_[x][y][z];
        z++;
        if (z == 30)
        {
            z = 0;
            y++;
        }
        if (y == 30)
        {
            y = 0;
            x++;
        }
        if (x == 30)
        {
            break;
            x++;
        }
    }
    cout<<signedDistMap_[0][0][0]<<endl;
    cout<<signedDistMap_[29][29][29]<<endl;
    inputFileStream.close();
}
