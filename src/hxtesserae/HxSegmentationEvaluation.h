#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortToggleList.h>

#include <hxfield/HxUniformScalarField3.h>

#include "api.h"

/**
 * This module computes the rand index and variation of information to compare a segmentation result
 * with a ground truth segmentation (or an arbitrary second segmentation).
 * Only voxels inside given mask are taken into calculation (mask value > 0).
 * Labelfields are unsigned 16 bit label fields.
 * Background voxels (value = 0) in the masked area are taken into calculation and treated as all other labels.
 * It does not matter if there are labels with no voxels inside the masked region.
 */
class HXTESSERAE_API HxSegmentationEvaluation : public HxCompModule
{
    HX_HEADER(HxSegmentationEvaluation);

public:
    void compute();
    void update();

    /// Tcl command interface.
    virtual int parse(Tcl_Interp* t, int argc, char** argv);

    HxConnection portGroundTruthSegmentation;
    HxConnection portMask;

    HxPortToggleList portUserSelection;
    HxPortDoIt portAction;

private:
    bool checkInputFields(HxUniformScalarField3* segmentationField,
                          HxUniformScalarField3* groundTruthField,
                          HxUniformScalarField3* maskField);

    double diceMeasureOnForeground(HxUniformScalarField3* segmentationField,
                                   HxUniformScalarField3* groundTruthField,
                                   HxUniformScalarField3* maskField);
    double variationOfInformation(HxUniformScalarField3* segmentationField,
                                  HxUniformScalarField3* groundTruthField,
                                  HxUniformScalarField3* maskField);
    double randIndex(HxUniformScalarField3* segmentationField,
                     HxUniformScalarField3* groundTruthField,
                     HxUniformScalarField3* maskField);
};
