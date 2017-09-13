#ifndef HXPLANEDISTANCEMAP_H
#define HXPLANEDISTANCEMAP_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortRadioBox.h>

#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformVectorField3.h>

#include <hxlandmark/internal/HxLandmarkSet.h>

#include <mclib/McMat3d.h>
#include <mclib/McVec4i.h>

#include <hxsurface/HxSurface.h>

#include "api.h"

/**
 * Computes a 2D planar distance map for a given foreground segmentation
 * (foreground voxels with label field values > 0).
 *
 * For each foreground label: Compute the plane that approximates the surrounding label field
 * in the best possible way. Then shoot rays inside this plane to compute the distance to the
 * nearest background voxel.
 *
 * Plane approximation starting from point P: The standard-solution is to shoot rays in all 3D directions and to
 * compute the intersection points of these rays with the background. The plane is then
 * the best fitting plane regarding these points (in least squares sense; and translated to go
 * through P). As a second option, it is possible to take all voxel positions
 * of foreground voxels around P (in a neighborhood) and do a PCA on these positions. The first 2
 * eigenvectors span the computed plane.
 *
 * Distance computation in plane: The distance can be computed one-directional, two-directional or
 * in an average sense.
 *
 * This module was used for distance map computation for the tesserae segmentation project.
 */

class HXTESSERAE_API HxPlaneDistanceMap : public HxCompModule
{
    HX_HEADER(HxPlaneDistanceMap);

public:
    void compute();
    void update();

    HxConnection portSphere;
    HxConnection portPlaneNormals;
    HxConnection portVisualizationLandmarks;
    HxConnection portMask;

    HxPortRadioBox portType;
    HxPortRadioBox portDistance;
    HxPortIntTextN portNumberOfRays;
    HxPortIntTextN portSizeBox;
    HxPortDoIt portDoIt;

private:
    void initOutputDistanceField();
    bool checkPreconditions();

    // Create distance field in plane
    void createDistanceField();

    // Calculate plane
    void getPlaneVector(const int i,
                        const int j,
                        const int k,
                        HxUniformVectorField3* vectorField,
                        HxUniformScalarField3* varianceField,
                        bool createNormals,
                        bool createVarianceField,
                        SbVec3f& orthogonalVector,
                        SbVec3f& planeVector);

    // Calculate distance in given plane
    float distanceToBackgroundInPlane(SbVec3f pos, SbVec3f planeDirectionVector1, SbVec3f planeNormalVector, bool setLandmarks);
    float bresenham(McVec3f bondEndPoints[2], McVec3f& intersectionPoint);

    float shortestOneDirectionalDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize);
    float shortestTwoDirectionalDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize);
    float averageDistance(SbVec3f position, SbVec3f directionVector, SbRotation rotation, bool setLandmarks, float stepsize);

    // Function that helps visualization
    void visualizeRayShootingInPlane();
    void visualizeRayShootingInPlaneForPoint(const int i, const int j, const int k);
    void visualizeRayShootingUsingSphere();
    void visualizeRayShootingUsingSphereForPoint(const int i, const int j, const int k);

    // PCA methods
    SbVec3d computeCenter(int x, int y, int z, HxUniformScalarField3* input);
    McMat3d computePrincipalAxes(int x, int y, int z, HxUniformScalarField3* input, SbVec3d& center);

    // Helper functions
    SbVec3f getOrthogonalVector(SbVec3f v);

    HxUniformScalarField3* m_backgroundSegmentation;
    HxUniformVectorField3* m_planeNormalsField;
    HxSurface* m_sphere;
    HxLandmarkSet* m_visualizationInputLandmarks;

    HxUniformScalarField3* m_outputDistanceField;

    HxLandmarkSet* m_rayCastingVisualizationLandmarks;
};

#endif
