#ifndef HX_MT_LANDMARK_EVALUATION
#define HX_MT_LANDMARK_EVALUATION

#include "api.h"
#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortRadioBox.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxlandmark/internal/HxLandmarkSet.h>

/**
  This module evaluates a segmentation, given as a uniform scalar field, against a set of landmarks. For each landmark,
  a corresponding voxel is computed. It is either the nearest voxel or the nearest foreground voxel (using BFS).
  The label field is evaluated at this voxel position and a counter for the label found there is increased (label was hit
  by the landmark). This counter value and several
  derived measures are output as a spreadsheet with one line per label. For a perfect segmentation, each hit counter (except for the background label)
  equals 1, i.e. each label is hit by exactly one marker. Additionally, there is a second table in the spreadsheet with summarized information
  containing data like the number of correct labels, numer of clusters, ..., and finally the precision and recall value.

*/
class HXINTERACTIVESEGMENTATION_API HxMTLandmarkEvaluation : public HxCompModule
{
    HX_HEADER(HxMTLandmarkEvaluation);

public:
    virtual void compute();

    HxPortRadioBox portNeighborhood;
    HxPortDoIt portDoIt;

private:
    HxConnection portLandmarks;
    void getSurroundingGridPoints(const McVec3f& pos, const HxUniformScalarField3& field, McDArray<McVec3i>& gridPoints);
    float getLabel(const McVec3f& pos, HxUniformScalarField3& field);
    float getNearestForegroundLabelBFS(const McVec3f& pos, HxUniformScalarField3* labels);
};

#endif
