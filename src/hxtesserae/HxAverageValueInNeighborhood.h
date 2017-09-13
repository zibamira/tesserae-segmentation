#ifndef HXAVERAGEVALUEINNEIGHBORHOOD_H
#define HXAVERAGEVALUEINNEIGHBORHOOD_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortModuleSwitch.h>
#include <hxcore/HxPortRadioBox.h>
#include <hxcore/HxPortButtonList.h>

#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformLabelField3.h>

#include <hxquant2/internal/HxQuant2GenericModule.h>

#include "api.h"

/**
 * @brief The HxAverageValueInNeighborhood class
 * Given a scalar field and a label field there are 2 ways to compute average values:
 * 1. Compute the average value in a neighborhood for all voxels inside the given label field.
 * The neighborhood contains all voxels that are in the label field and in a cube around the current voxel.
 * The port size determines the size of this cube (sidelength / 2).
 * The output is a label field
 * 2. Compute average values in neighborhoods around the labels (one value per label).
 * The result is written in a spreadheet.
 * The neighborhood is computed via dilation (using ball, size is taken from size port and equals
 * half-kernel size).
 * The average is computed for all voxels inside the dilation and not inside the original label,
 * additionaly only values larger than a given threshold are taken into calculation.
 */
class HXTESSERAE_API HxAverageValueInNeighborhood : public HxCompModule
{
    HX_HEADER(HxAverageValueInNeighborhood);

public:
    void compute();
    void update();

    HxConnection portLabelField;

    HxPortRadioBox portSelection;
    HxPortIntTextN portSize;
    HxPortFloatTextN portThreshold;
    HxPortFloatTextN portVisualizeCoords;
    HxPortButtonList portVisualize;
    HxPortDoIt portDoIt;

private:
    /*
     * Algorithm 1: Available in 2 versions: one slow and one fast
     * Slow: for each voxel iterate through complete neighborhood
     * Fast: iterate through complete neighborhood only for 0,0,0; afterwards each voxel uses previous results from
     * neighbored refererence voxel
     */

    // Slow version
    void computeAverageValuesInsideLabel(HxUniformScalarField3* scalarField, HxUniformScalarField3* labelField);
    float getNeighborhoodValue(HxUniformScalarField3* scalarField,
                               HxUniformScalarField3* labelField,
                               const int x,
                               const int y,
                               const int z,
                               int& numberOfVoxels,
                               float& totalValue);

    // Fast version:
    // we compute the
    // neighborhood average values for the first xy-slice. Then we iterate over increasing
    // z-values in parallel on the CPU, i.e. with one thread for one pair of x and y
    // coordinates.
    void computeAverageValuesInsideLabelFastAlgorithm(HxUniformScalarField3* scalarField, HxUniformScalarField3* labelField);
    float getNeighborhoodValueUsingReferenceX(HxUniformScalarField3* scalarField,
                                              HxUniformScalarField3* labelField,
                                              const int x,
                                              const int y,
                                              const int z,
                                              int numberVoxelsReference,
                                              float totalValueReference,
                                              int& numberOfVoxels,
                                              float& totalValue);
    float getNeighborhoodValueUsingReferenceY(HxUniformScalarField3* scalarField,
                                              HxUniformScalarField3* labelField,
                                              const int x,
                                              const int y,
                                              const int z,
                                              int numberVoxelsReference,
                                              float totalValueReference,
                                              int& numberOfVoxels,
                                              float& totalValue);
    float getNeighborhoodValueUsingReferenceZ(HxUniformScalarField3* scalarField,
                                              HxUniformScalarField3* labelField,
                                              const int x,
                                              const int y,
                                              const int z,
                                              int numberVoxelsReference,
                                              float totalValueReference,
                                              int& numberOfVoxels,
                                              float& totalValue);

    // Only visualization (made for paper)
    void visualizeNeighborhoodInsideLabel(HxUniformLabelField3* labelField,
                                          const int x,
                                          const int y,
                                          const int z);

    /*
     * Algorithm 2:
     * Given a labelfield, dilate each label individually and compute the average value
     * in the new dilated area (without input label area)
     */
    void computeAverageValuesAroundAllLabels(HxUniformScalarField3* scalarField,
                                             HxUniformLabelField3* labelField);
    float computeAverageValueAroundLabel(HxUniformScalarField3* scalarField,
                                         HxUniformLabelField3* labelField,
                                         HxQuant2GenericModule* dilationModule,
                                         const float averageValueInsideLabel,
                                         int& numberOfVoxels);
};

#endif // HXAVERAGEVALUEINNEIGHBORHOOD_H
