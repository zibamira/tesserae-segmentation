#include "HxPlaneDistanceMap.h"

#include <hxcore/HxObjectPool.h>
#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxcore/HxSettingsMgr.h>

#include <hxfield/HxUniformLabelField3.h>

#include <Inventor/SbMatrix.h>
#include <Inventor/SbRotation.h>

#include <mclib/internal/McDMatrix.h>
#include <mclib/internal/McWatch.h>
#include <mclib/internal/McFitPlane.h>

#ifdef _OPENMP
#include <omp.h>
#endif

HX_INIT_CLASS(HxPlaneDistanceMap, HxCompModule)

HxUniformScalarField3* m_backgroundSegmentation;
HxUniformVectorField3* m_planeNormalsField;
HxSurface* m_sphere;
HxLandmarkSet* m_visualizationInputLandmarks;

HxUniformScalarField3* m_outputDistanceField;

HxLandmarkSet* m_rayCastingVisualizationLandmarks;

HxPlaneDistanceMap::HxPlaneDistanceMap()
    : HxCompModule(HxUniformScalarField3::getClassTypeId())
    , portSphere(this, "sphere", tr("Sphere"), HxSurface::getClassTypeId())
    , portPlaneNormals(this, "planeNormals", tr("Plane Normals"), HxUniformVectorField3::getClassTypeId())
    , portVisualizationLandmarks(this, "visualizationlandmarks", tr("Visualization Landmarks"), HxLandmarkSet::getClassTypeId())
    , portMask(this, "Mask", tr("Mask"), HxUniformLabelField3::getClassTypeId())
    , portType(this, "type", tr("Type"), 3)
    , portDistance(this, "distance", tr("Distance"), 3)
    , portNumberOfRays(this, "numberOfRays", tr("Number Of Rays"), 1)
    , portSizeBox(this, "sizebox", tr("Sizebox"), 1)
    , portDoIt(this, "doIt", tr("Action"))
{
    m_backgroundSegmentation = 0;
    m_planeNormalsField = 0;
    m_sphere = 0;
    m_visualizationInputLandmarks = 0;
    m_outputDistanceField = 0;
    m_rayCastingVisualizationLandmarks = 0;

    portNumberOfRays.setLabel(0, "In plane");
    portNumberOfRays.setValue(0, 361);

    portType.setLabel(0, "No sphere");
    portType.setLabel(1, "Sphere");
    portType.setLabel(2, "Landmark sphere visualization");
    portType.setValue(1);

    portDistance.setLabel(0, "one-directional");
    portDistance.setLabel(1, "two-directional");
    portDistance.setLabel(2, "average");

    portSizeBox.setValue(0, 40);
}

HxPlaneDistanceMap::~HxPlaneDistanceMap()
{
}

void
HxPlaneDistanceMap::update()
{
    if (portData.isNew(HxData::NEW_SOURCE))
    {
        m_backgroundSegmentation = hxconnection_cast<HxUniformScalarField3>(portData);
    }
    if (portVisualizationLandmarks.isNew(HxData::NEW_SOURCE))
    {
        m_visualizationInputLandmarks = hxconnection_cast<HxLandmarkSet>(portVisualizationLandmarks);
    }
    if (portPlaneNormals.isNew(HxData::NEW_SOURCE))
    {
        m_planeNormalsField = hxconnection_cast<HxUniformVectorField3>(portPlaneNormals);
    }
    if (portSphere.isNew(HxData::NEW_SOURCE))
    {
        m_sphere = hxconnection_cast<HxSurface>(portSphere);
    }
}

void
HxPlaneDistanceMap::compute()
{
    if (!checkPreconditions())
    {
        return;
    }

    if (portDoIt.wasHit())
    {
        if (portType.getValue() == 2)
        {
            visualizeRayShootingUsingSphere();
        }
        else
        {
            initOutputDistanceField();
            createDistanceField();
        }
    }
}

void
HxPlaneDistanceMap::initOutputDistanceField()
{
    theWorkArea->startWorking("Initializing data ...");

    const McDim3l& dims = m_backgroundSegmentation->lattice().getDims();
    const McBox3f& bbox = m_backgroundSegmentation->getBoundingBox();

    const McPrimType outputFieldPrimType = McPrimType::MC_FLOAT;

    m_outputDistanceField = dynamic_cast<HxUniformScalarField3*>(getResult());

    bool initializeWithZero = false;

    if (m_outputDistanceField == 0)
    {
        m_outputDistanceField = new HxUniformScalarField3(dims, outputFieldPrimType);
        initializeWithZero = true;
    }
    else
    {
        const McDim3l& dimsOutputField = m_outputDistanceField->lattice().getDims();
        if ((m_outputDistanceField->lattice().primType() != outputFieldPrimType) ||
            (dims[0] != dimsOutputField[0]) || (dims[1] != dimsOutputField[1]) ||
            (dims[2] != dimsOutputField[2]))
        {
            m_outputDistanceField->lattice().setPrimType(outputFieldPrimType);
            m_outputDistanceField->lattice().resize(dims);
            initializeWithZero = true;
        }
    }

    m_outputDistanceField->lattice().setBoundingBox(bbox);

    if (initializeWithZero)
    {
        const mculong size = (mculong)dims[0] * (mculong)dims[1] * (mculong)dims[2];

        float* data = (float*)m_outputDistanceField->lattice().dataPtr();

        const float background = 0.0;
        for (mculong i = 0; i < size; ++i)
            data[i] = background;
    }

    m_outputDistanceField->composeLabel(m_backgroundSegmentation->getLabel(), "distance");

    theWorkArea->stopWorking();
}

bool
HxPlaneDistanceMap::checkPreconditions()
{
    if (!m_backgroundSegmentation)
    {
        return false;
    }

    if (portType.getValue() == 0)
    {
        return true;
    }
    else if (portType.getValue() == 1)
    {
        return m_sphere;
    }
    else
    {
        return (m_visualizationInputLandmarks);
    }
}

// Calculate distance of every foreground voxel to the background in the PCA plane
void
HxPlaneDistanceMap::createDistanceField()
{
    theWorkArea->startWorking("Create distance field ...");
    McWatch watch;

    const McDim3l& dims = m_backgroundSegmentation->lattice().getDims();
    const McBox3f& bbox = m_backgroundSegmentation->getBoundingBox();

    // Debug variables
    const bool createNormals = false;
    const bool createVarianceField = false;

    HxUniformVectorField3* vectorField = 0;
    if (createNormals && !m_planeNormalsField)
    {
        vectorField = new HxUniformVectorField3(dims, McPrimType::MC_FLOAT);
        vectorField->lattice().setBoundingBox(bbox);
    }

    HxUniformScalarField3* varianceField = 0;
    if ((createVarianceField) && (!m_planeNormalsField) && (portType.getValue() == 1))
    {
        varianceField = new HxUniformScalarField3(dims, McPrimType::MC_FLOAT);
        varianceField->lattice().setBoundingBox(bbox);
    }

    HxUniformLabelField3* maskField = hxconnection_cast<HxUniformLabelField3>(portMask);

    int numIteration = 0;

// Iterate over all foreground voxels
#ifdef _OPENMP
    int numThreads = theSettingsMgr->getPreferences().maxNumberOfComputeThreads;
    if (numThreads < 1)
        numThreads = 1;
    omp_set_num_threads(numThreads);
#pragma omp parallel for
#endif
    for (int i = 0; i < dims[0]; i++)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
        numIteration++;
#ifdef _OPENMP
        if (omp_get_thread_num() == 0) // Main thread.
#endif
        {
            theWorkArea->setProgressValue(numIteration / (float)(dims[0]));
        }

        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float label;
                float maskValue = 1.0;
                m_backgroundSegmentation->lattice().eval(i, j, k, &label);
                if (maskField)
                {
                    maskField->lattice().eval(i, j, k, &maskValue);
                }

                if ((label > 0.0) && (maskValue > 0.0))
                {
                    SbVec3f planeVector(0, 0, 0);
                    SbVec3f orthogonalVector(0, 0, 0);
                    getPlaneVector(i, j, k, vectorField, varianceField, createNormals, createVarianceField, orthogonalVector, planeVector);

                    // Set distance value
                    McVec3f globalPositionVector = m_backgroundSegmentation->lattice().coords()->pos(McVec3i(i, j, k));
                    float distance = distanceToBackgroundInPlane(globalPositionVector, planeVector, orthogonalVector, false);
                    m_outputDistanceField->lattice().set(i, j, k, &distance);
                }
            }
        }
    }

    setResult(m_outputDistanceField);
    if (vectorField)
    {
        vectorField->setLabel("VectorField");
        theObjectPool->addObject(vectorField);
    }
    if (varianceField)
    {
        varianceField->setLabel("VarianceField");
        theObjectPool->addObject(varianceField);
    }

    theMsg->printf("time: %f seconds", watch.stop());
    theWorkArea->stopWorking();
}

void
HxPlaneDistanceMap::getPlaneVector(const int i,
                                   const int j,
                                   const int k,
                                   HxUniformVectorField3* vectorField,
                                   HxUniformScalarField3* varianceField,
                                   bool createNormals,
                                   bool createVarianceField,
                                   SbVec3f& orthogonalVector,
                                   SbVec3f& planeVector)
{
    if (!m_planeNormalsField)
    {
        if (createNormals)
        {
            float zeroValues[3];
            zeroValues[0] = 0;
            zeroValues[1] = 0;
            zeroValues[2] = 0;
            vectorField->lattice().set(i, j, k, zeroValues);
        }
        if (portType.getValue() == 1)
        {
            // Use sphere
            const int numberOfRays = m_sphere->points().sizeInt();
            McDMatrix<float> m(3, numberOfRays);

            McVec3f globalPositionVector = m_backgroundSegmentation->lattice().coords()->pos(McVec3i(i, j, k));
            McVec3f endPoints[2];
            endPoints[0] = globalPositionVector;

            float averageDist = 0.0;

            int numberOfSuccessfulRays = 0;

            McFitPlane plane;

            for (int ii = 0; ii < numberOfRays; ii++)
            {
                McVec3f spherePoint = m_sphere->points()[ii];
                endPoints[1] = globalPositionVector + spherePoint;

                McVec3f backgroundIntersectionPoint;
                float dist = bresenham(endPoints, backgroundIntersectionPoint);
                if (dist != 100000)
                {
                    // Add as ii-th column
                    m[0][numberOfSuccessfulRays] = backgroundIntersectionPoint[0];
                    m[1][numberOfSuccessfulRays] = backgroundIntersectionPoint[1];
                    m[2][numberOfSuccessfulRays] = backgroundIntersectionPoint[2];

                    plane.addPoint(backgroundIntersectionPoint);

                    averageDist += dist;
                    numberOfSuccessfulRays++;
                }
            }
            averageDist /= (float)numberOfSuccessfulRays;
            //m_outputDistanceField->lattice().set(i, j, k, &averageDist);

            m.resizeKeep(3, numberOfSuccessfulRays);

            // Fit plane to intersection points
            plane.fit();
            McVec3f orthogonalVectorMc = plane.getNormal();
            orthogonalVector[0] = orthogonalVectorMc[0];
            orthogonalVector[1] = orthogonalVectorMc[1];
            orthogonalVector[2] = orthogonalVectorMc[2];
            planeVector = getOrthogonalVector(orthogonalVector);

            if (createVarianceField)
            {
                float totalSquaredDiff = 0;
                float variance = 0;
                for (int ii = 0; ii < numberOfSuccessfulRays; ii++)
                {
                    McVec3f point;
                    point[0] = m[0][ii];
                    point[1] = m[1][ii];
                    point[2] = m[2][ii];

                    // Calculate squared difference to plane
                    McVec3f planeCenter = plane.getCentroid();
                    orthogonalVectorMc.normalize();
                    float diff = abs(planeCenter.dot(orthogonalVectorMc) - point.dot(orthogonalVectorMc));
                    float squaredDiff = diff * diff;
                    totalSquaredDiff += squaredDiff;
                }
                totalSquaredDiff /= numberOfSuccessfulRays;
                for (int ii = 0; ii < numberOfSuccessfulRays; ii++)
                {
                    McVec3f point;
                    point[0] = m[0][ii];
                    point[1] = m[1][ii];
                    point[2] = m[2][ii];

                    // Calculate variance
                    McVec3f planeCenter = plane.getCentroid();
                    orthogonalVectorMc.normalize();
                    float diff = abs(planeCenter.dot(orthogonalVectorMc) - point.dot(orthogonalVectorMc));
                    float squaredDiff = diff * diff;
                    variance += (squaredDiff - totalSquaredDiff) * (squaredDiff - totalSquaredDiff);
                }

                variance /= numberOfSuccessfulRays;
                varianceField->lattice().set(i, j, k, &variance);
            }
        }
        else
        {
            // No sphere
            SbVec3d center;

            // Compute principal axes
            McMat3d axesModel = computePrincipalAxes(i, j, k, m_backgroundSegmentation, center);

            // Set plane vectors
            planeVector = SbVec3f(axesModel[0][0], axesModel[1][0], axesModel[2][0]);
            orthogonalVector = SbVec3f(axesModel[0][2], axesModel[1][2], axesModel[2][2]);
        }

        // Visualize orthogonal vector
        if (createNormals)
        {
            orthogonalVector.normalize();
            float orthogonalValues[3];
            orthogonalValues[0] = orthogonalVector[0];
            orthogonalValues[1] = orthogonalVector[1];
            orthogonalValues[2] = orthogonalVector[2];
            vectorField->lattice().set(i, j, k, orthogonalValues);
        }
    }
    else
    {
        // There are normals
        float orthogonalVectorValues[3];
        m_planeNormalsField->lattice().eval(i, j, k, orthogonalVectorValues);
        orthogonalVector[0] = orthogonalVectorValues[0];
        orthogonalVector[1] = orthogonalVectorValues[1];
        orthogonalVector[2] = orthogonalVectorValues[2];
        planeVector = getOrthogonalVector(orthogonalVector);
    }
}

// Cast rays from position in given plane
// Calculate shortest distance of a ray to a background voxel
// Input coordinates are global coordinates (not mesh coords)
// If setLandmarks is true, landmarkSet gets a landmark for each ray
float
HxPlaneDistanceMap::distanceToBackgroundInPlane(SbVec3f position, SbVec3f planeDirectionVector1, SbVec3f planeNormalVector, bool setLandmarks)
{
    float eps = 0.0001;
    if (planeNormalVector.length() < eps)
    {
        return 100000.0;
    }

    int steps = portNumberOfRays.getValue(0); // number of rays
    float stepsize = (2 * M_PI) / (steps - 1);

    SbRotation rotation(planeNormalVector, stepsize);

    SbVec3f directionVector(planeDirectionVector1);

    if (portDistance.getValue() == 0)
    {
        return shortestOneDirectionalDistance(position, directionVector, rotation, setLandmarks, stepsize);
    }
    else if (portDistance.getValue() == 1)
    {
        return shortestTwoDirectionalDistance(position, directionVector, rotation, setLandmarks, stepsize);
    }
    else if (portDistance.getValue() == 2)
    {
        return averageDistance(position, directionVector, rotation, setLandmarks, stepsize);
    }

    return 0.f;
}

// Returns 100000 if ray hits boundary first
float
HxPlaneDistanceMap::bresenham(McVec3f bondEndPoints[2], McVec3f& intersectionPoint)
{
    const McVec3f voxelSizeVector = m_backgroundSegmentation->getVoxelSize();
    const float voxelSize = voxelSizeVector[0];

    const McDim3l& dims = m_backgroundSegmentation->lattice().getDims();

    const McBox3f& bbox = m_backgroundSegmentation->getBoundingBox();

    // find major axis
    int major = 0, minor[2] = { 0, 0 };
    float maxDist = fabs(bondEndPoints[0][0] - bondEndPoints[1][0]) / voxelSize;
    float tmpDist = fabs(bondEndPoints[0][1] - bondEndPoints[1][1]) / voxelSize;
    if (tmpDist > maxDist)
    {
        major = 1;
        maxDist = tmpDist;
    }
    tmpDist = fabs(bondEndPoints[0][2] - bondEndPoints[1][2]) / voxelSize;
    if (tmpDist > maxDist)
    {
        major = 2;
        maxDist = tmpDist;
    }

    // assign minor axes
    switch (major)
    {
        case 0: // x-axis
            minor[0] = 1;
            minor[1] = 2;
            break;
        case 1: // y-axis
            minor[0] = 0;
            minor[1] = 2;
            break;
        case 2: // z-axis
            minor[0] = 0;
            minor[1] = 1;
            break;
    }

    float majorDist = bondEndPoints[1][major] - bondEndPoints[0][major];

    int majorStepsize;

    int id;
    McVec3f ep[2];
    if (majorDist > 0.0)
    {
        ep[0] = bondEndPoints[0];
        ep[1] = bondEndPoints[1];
        majorStepsize = 1;
    }
    else
    {
        ep[0] = bondEndPoints[0];
        ep[1] = bondEndPoints[1];
        majorStepsize = -1;
        majorDist *= -1.0;
    }

    // majorDist always positive

    // get lengths of projections
    float minorDist[2];
    minorDist[0] = ep[1][minor[0]] - ep[0][minor[0]];
    minorDist[1] = ep[1][minor[1]] - ep[0][minor[1]];

    // compute slopes in projection planes
    float m[2];
    m[0] = minorDist[0] / majorDist;
    m[1] = minorDist[1] / majorDist;

    // extend cylinder line to grid plane
    id = (int)((ep[0][major] - bbox[2 * major]) / voxelSize);
    float factor = ((bbox[2 * major] + id * voxelSize) - ep[0][major]);
    ep[0][major] = bbox[2 * major] + id * voxelSize;
    if (majorStepsize == 1)
    {
        ep[0][minor[0]] += factor * m[0];
        ep[0][minor[1]] += factor * m[1];
    }
    else
    {
        ep[0][minor[0]] -= factor * m[0];
        ep[0][minor[1]] -= factor * m[1];
    }

    McVec3f tmp;
    tmp[major] = ep[0][major] + 0.5 * voxelSize;
    tmp[minor[0]] = ep[0][minor[0]] + 0.5 * voxelSize;
    tmp[minor[1]] = ep[0][minor[1]] + 0.5 * voxelSize;

    // compute start node of grid
    McVec4i curId;
    curId[major] = (int)((tmp[major] - bbox[2 * major]) / voxelSize);
    curId[minor[0]] = (int)((tmp[minor[0]] - bbox[2 * minor[0]]) / voxelSize);
    curId[minor[1]] = (int)((tmp[minor[1]] - bbox[2 * minor[1]]) / voxelSize);

    // compute start error
    float error[2];
    error[0] = (ep[0][minor[0]] -
                (bbox[2 * minor[0]] + curId[minor[0]] * voxelSize)) /
               voxelSize;
    error[1] = (ep[0][minor[1]] -
                (bbox[2 * minor[1]] + curId[minor[1]] * voxelSize)) /
               voxelSize;

    McVec3f curCoord = McVec3f(bbox[0], bbox[2], bbox[4]);
    curCoord[0] += ((float)curId[0]) * voxelSize;
    curCoord[1] += ((float)curId[1]) * voxelSize;
    curCoord[2] += ((float)curId[2]) * voxelSize;
    // store first node
    if (!(curId[0] < 0 || curId[1] < 0 || curId[2] < 0 ||
          curId[0] >= dims[0] || curId[1] >= dims[1] || curId[2] >= dims[2]))
    {
        curId[3] = curId[0] + dims[0] * (curId[1] + curId[2] * dims[1]);
    }

    // get points along cylinder line
    float exact;
    for (int i = 0;; i++)
    {
        curId[major] += majorStepsize;

        exact = curId[minor[0]] + m[0] + error[0];
        curId[minor[0]] = (int)(exact + 0.5);
        error[0] = exact - (float)curId[minor[0]];

        exact = curId[minor[1]] + m[1] + error[1];
        curId[minor[1]] = (int)(exact + 0.5);
        error[1] = exact - (float)curId[minor[1]];

        curId[3] = curId[0] + dims[0] * (curId[1] + curId[2] * dims[1]);

        curCoord = McVec3f(bbox[0], bbox[2], bbox[4]);
        curCoord[0] += curId[0] * voxelSize;
        curCoord[1] += curId[1] * voxelSize;
        curCoord[2] += curId[2] * voxelSize;

        if (curId[0] < 0 || curId[1] < 0 || curId[2] < 0 || curId[0] >= dims[0] || curId[1] >= dims[1] || curId[2] >= dims[2])
        {
            intersectionPoint = McVec3f(0, 0, 0);
            return 100000;
        }

        float value;
        m_backgroundSegmentation->lattice().eval(curId[0], curId[1], curId[2], &value);
        if (value == 0.0)
        {
            intersectionPoint = curCoord;
            McVec3f distanceVector = curCoord - ep[0];
            return distanceVector.length();
        }
    }
}

// Compute shortest one-directional distance
float
HxPlaneDistanceMap::shortestOneDirectionalDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize)
{
    float distance = 0.f;

    // Cast rays
    for (float angle = 0; angle < 2 * M_PI; angle += stepsize)
    {
        // Construct ray: rotate planeDirectionVector1 with angle around planeNormalVector
        if (angle > 0)
        {
            rotation.multVec(directionVector, directionVector);
        }

        // Calculate distance to background
        McVec3f positionMc(position[0], position[1], position[2]);
        McVec3f directionVectorMc(directionVector[0], directionVector[1], directionVector[2]);
        directionVectorMc.normalize();
        McVec3f positionMc2 = positionMc + directionVectorMc;
        McVec3f twoPositions[2];
        twoPositions[0] = positionMc;
        twoPositions[1] = positionMc2;

        if (setLandmarks)
        {
            McVec3f positionLandmark = positionMc + 100 * directionVectorMc;
            m_rayCastingVisualizationLandmarks->appendLandmark(positionLandmark);
        }

        float currentDistance = bresenham(twoPositions, positionMc2);
        if ((distance > currentDistance) || (angle == 0))
        {
            distance = currentDistance;
        }
    }

    return distance;
}

// Compute shortest two-directional distance
float
HxPlaneDistanceMap::shortestTwoDirectionalDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize)
{
    float distance = 0.f;

    // Cast rays
    for (float angle = 0; angle < M_PI; angle += stepsize)
    {
        // Construct ray: rotate planeDirectionVector1 with angle around planeNormalVector
        if (angle > 0)
        {
            rotation.multVec(directionVector, directionVector);
        }

        // Calculate distance to background
        McVec3f positionMc(position[0], position[1], position[2]);
        McVec3f directionVectorMc(directionVector[0], directionVector[1], directionVector[2]);
        directionVectorMc.normalize();
        McVec3f positionMc2 = positionMc + directionVectorMc;
        McVec3f positionMc3 = positionMc - directionVectorMc;

        if (setLandmarks)
        {
            McVec3f positionLandmark = positionMc + 100 * directionVectorMc;
            m_rayCastingVisualizationLandmarks->appendLandmark(positionLandmark);
        }

        McVec3f twoPositions[2];
        twoPositions[0] = positionMc;
        twoPositions[1] = positionMc2;
        float currentDistance = bresenham(twoPositions, positionMc2);

        twoPositions[1] = positionMc3;
        currentDistance += bresenham(twoPositions, positionMc3);

        if ((distance > currentDistance) || (angle == 0))
        {
            distance = currentDistance;
        }
    }

    return distance;
}

// Compute average distance
float
HxPlaneDistanceMap::averageDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize)
{
    float distance = 0.f;
    float numDistances = 0.f;

    // Cast rays
    for (float angle = 0; angle < 2 * M_PI; angle += stepsize)
    {
        // Construct ray: rotate planeDirectionVector1 with angle around planeNormalVector
        if (angle > 0)
        {
            rotation.multVec(directionVector, directionVector);
        }

        // Calculate distance to background
        McVec3f positionMc(position[0], position[1], position[2]);
        McVec3f directionVectorMc(directionVector[0], directionVector[1], directionVector[2]);
        directionVectorMc.normalize();
        McVec3f positionMc2 = positionMc + directionVectorMc;
        McVec3f twoPositions[2];
        twoPositions[0] = positionMc;
        twoPositions[1] = positionMc2;

        if (setLandmarks)
        {
            McVec3f positionLandmark = positionMc + 100 * directionVectorMc;
            m_rayCastingVisualizationLandmarks->appendLandmark(positionLandmark);
        }

        float currentDistance = bresenham(twoPositions, positionMc2);
        if (currentDistance != 100000)
        {
            distance += currentDistance;
            numDistances += 1.f;
        }
    }

    if (numDistances > 0.f)
        distance = distance / numDistances;

    return distance;
}

// Shows rays in plane
void
HxPlaneDistanceMap::visualizeRayShootingInPlane()
{
    m_rayCastingVisualizationLandmarks = HxLandmarkSet::createInstance();

    for (int i = 0; i < m_visualizationInputLandmarks->getNumMarkers(); i++)
    {
        McVec3f pos = m_visualizationInputLandmarks->getPosition(i);
        HxLoc3Regular* location = m_planeNormalsField->lattice().coords()->createLocation();
        location->set(pos);
        const int ii = location->getIx();
        const int jj = location->getIy();
        const int kk = location->getIz();

        visualizeRayShootingInPlaneForPoint(ii, jj, kk);
    }

    theObjectPool->addObject(m_rayCastingVisualizationLandmarks);
}

void
HxPlaneDistanceMap::visualizeRayShootingInPlaneForPoint(const int i, const int j, const int k)
{
    SbVec3f globalPositionVector = m_backgroundSegmentation->lattice().coords()->pos(McVec3i(i, j, k));
    float orthogonalVectorValues[3];
    m_planeNormalsField->lattice().eval(i, j, k, orthogonalVectorValues);
    SbVec3f orthogonalVector(orthogonalVectorValues);
    SbVec3f planeVector1 = getOrthogonalVector(orthogonalVector);
    distanceToBackgroundInPlane(globalPositionVector, planeVector1, orthogonalVector, true);
}

// Shows background intersection points
void
HxPlaneDistanceMap::visualizeRayShootingUsingSphere()
{
    m_rayCastingVisualizationLandmarks = HxLandmarkSet::createInstance();

    for (int i = 0; i < m_visualizationInputLandmarks->getNumMarkers(); i++)
    {
        McVec3f pos = m_visualizationInputLandmarks->getPosition(i);
        HxLoc3Regular* location = m_backgroundSegmentation->lattice().coords()->createLocation();
        location->set(pos);
        const int ii = location->getIx();
        const int jj = location->getIy();
        const int kk = location->getIz();

        visualizeRayShootingUsingSphereForPoint(ii, jj, kk);
    }

    theObjectPool->addObject(m_rayCastingVisualizationLandmarks);
}

void
HxPlaneDistanceMap::visualizeRayShootingUsingSphereForPoint(const int i, const int j, const int k)
{
    // Use sphere
    const int numberOfRays = m_sphere->points().sizeInt();
    McDMatrix<float> m(3, numberOfRays);

    McVec3f globalPositionVector = m_backgroundSegmentation->lattice().coords()->pos(McVec3i(i, j, k));
    McVec3f endPoints[2];
    endPoints[0] = globalPositionVector;

    for (int ii = 0; ii < numberOfRays; ii++)
    {
        McVec3f spherePoint = m_sphere->points()[ii];
        endPoints[1] = globalPositionVector + spherePoint;

        McVec3f backgroundIntersectionPoint;
        bresenham(endPoints, backgroundIntersectionPoint);

        // Show background intersection point
        m_rayCastingVisualizationLandmarks->appendLandmark(backgroundIntersectionPoint);
    }
}

// Taken from HxAlign:
SbVec3d
HxPlaneDistanceMap::computeCenter(int x, int y, int z, HxUniformScalarField3* input)
{
    int i, j, k;
    const McDim3l& dims = input->lattice().getDims();
    HxCoord3* coords = input->lattice().coords();

    SbVec3d center(0.0, 0.0, 0.0);

    int s = 20;

    mclong numberOfPoints = 0;
    for (k = z - s; k <= z + s; k++)
    {
        for (j = y - s; j <= y + s; j++)
        {
            for (i = x - s; i <= x + s; i++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]) && (k >= 0) && (k < dims[2]))
                {
                    double a = input->evalReg(i, j, k);
                    if (a > 0.0)
                    {
                        SbVec3f tmp = coords->pos(McVec3i(i, j, k));
                        SbVec3d tmpd(tmp);
                        center += tmpd;
                        numberOfPoints++;
                    }
                }
            }
        }
    }

    if (numberOfPoints)
        center /= numberOfPoints;
    return center;
}

// Taken from HxAlign:
McMat3d
HxPlaneDistanceMap::computePrincipalAxes(int x, int y, int z, HxUniformScalarField3* input, SbVec3d& center)
{
    int i, j, k;

    int s = portSizeBox.getValue(0); //40;

    const McDim3l& dims = input->lattice().getDims();
    HxCoord3* coords = input->lattice().coords();

    McDMatrix<double> inertia(3, 3);
    inertia.fill(0.0);

    center = computeCenter(x, y, z, input);

    for (k = z - s; k <= z + s; k++)
    {
        for (j = y - s; j <= y + s; j++)
        {
            for (i = x - s; i <= x + s; i++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]) && (k >= 0) && (k < dims[2]))
                {
                    double a = input->evalReg(i, j, k);
                    if (a > 0.0)
                    {
                        SbVec3f tmpf = coords->pos(McVec3i(i, j, k));
                        SbVec3d tmp(tmpf);
                        tmp -= center;

                        inertia[0][0] += tmp[0] * tmp[0]; //a*(tmp[1]*tmp[1] + tmp[2]*tmp[2]);
                        inertia[1][1] += tmp[1] * tmp[1]; // a*(tmp[0]*tmp[0] + tmp[2]*tmp[2]);
                        inertia[2][2] += tmp[2] * tmp[2]; // a*(tmp[1]*tmp[1] + tmp[0]*tmp[0]);

                        inertia[0][1] += tmp[0] * tmp[1]; //-= a*tmp[0]*tmp[1];
                        inertia[1][0] = inertia[0][1];

                        inertia[0][2] += tmp[0] * tmp[2]; //-= a*tmp[0]*tmp[2];
                        inertia[2][0] = inertia[0][2];

                        inertia[1][2] += tmp[1] * tmp[2]; //-= a*tmp[1]*tmp[2];
                        inertia[2][1] = inertia[1][2];
                    }
                }
            }
        }
    }

    // Compute eigenvectors and values
    McDMatrix<double> pa(3, 3);
    double v[3];
    inertia.eigenSystem(pa, v);

    // Sort eigenvectors after (descending) magnitude of eigenvalue
    int mo[3];
    for (i = 0; i < 3; i++)
    {
        int count = 0;
        for (j = 0; j < 3; j++)
        {
            if (j != i && v[j] > v[i])
                count++;
        }
        mo[i] = count;
    }
    // Account for equal eigenvalues
    for (i = 0; i < 3; i++)
    {
        for (j = i + 1; j < 3; j++)
        {
            if (v[i] == v[j])
                mo[j]++;
        }
    }

    // Setup matrix of sorted eigenvectors
    McMat3d m = McMat3d::IDENTITY;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            m[i][mo[j]] = pa(i, j);
        }
    }

    return m;
}

SbVec3f
HxPlaneDistanceMap::getOrthogonalVector(SbVec3f v)
{
    SbVec3f orthogonalVector(0, 0, 0);

    float eps = 0.0001;

    if (v.length() < eps)
    {
        //theMsg->printf("Warning: zero vector in orthogonal vector calculation");
        return orthogonalVector;
    }

    int index = 0;
    for (int i = 0; i < 3; i++)
    {
        if (fabs(v[i]) > eps)
        {
            index = i;
            break;
        }
    }

    orthogonalVector[(index + 1) % 3] = 1;
    orthogonalVector[(index + 2) % 3] = 1;
    orthogonalVector[index] = (-v[(index + 1) % 3] - v[(index + 2) % 3]) / (v[index]);

    return orthogonalVector;
}
