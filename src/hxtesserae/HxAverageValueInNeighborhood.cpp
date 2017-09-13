#include "HxAverageValueInNeighborhood.h"

#include <hxcore/HxResource.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxcore/HxObjectPool.h>
#include <hxcore/HxMessage.h>
#include <hxcore/HxSettingsMgr.h>

#include <mclib/internal/McWatch.h>

#include <hxspreadsheet/internal/HxSpreadSheet.h>

#ifdef _OPENMP
#include <omp.h>
#endif

HX_INIT_CLASS(HxAverageValueInNeighborhood, HxCompModule)

HxAverageValueInNeighborhood::HxAverageValueInNeighborhood()
    : HxCompModule(HxUniformScalarField3::getClassTypeId())
    , portLabelField(this, "labelField", tr("Label Field"), HxUniformLabelField3::getClassTypeId())
    , portSelection(this, "selection", tr("Selection"), 2)
    , portSize(this, "size", tr("Size"), 1)
    , portThreshold(this, "threshold", tr("Threshold"), 1)
    , portVisualizeCoords(this, "visualizeCoords", tr("Visualize Coords"), 3)
    , portVisualize(this, "visualize", tr("Visualize"), 1)
    , portDoIt(this, "doIt", tr("Action"))
{
    portSize.setValue(0, 15);

    portSelection.setLabel(0, "Inside Label");
    portSelection.setLabel(1, "Around Labels");

    portVisualize.setLabel(0, "Inside Label");

    portVisualizeCoords.setLabel(0, "x");
    portVisualizeCoords.setLabel(1, "y");
    portVisualizeCoords.setLabel(2, "z");
}

HxAverageValueInNeighborhood::~HxAverageValueInNeighborhood()
{
}

void
HxAverageValueInNeighborhood::update()
{
    if (portSelection.isNew())
    {
        if (portSelection.getValue() == 0)
        {
            portThreshold.hide();
        }
        else
        {
            portThreshold.show();
        }
    }
}

void
HxAverageValueInNeighborhood::compute()
{
    if (portDoIt.wasHit())
    {
        HxUniformScalarField3* scalarField = hxconnection_cast<HxUniformScalarField3>(portData);
        HxUniformLabelField3* labelField = hxconnection_cast<HxUniformLabelField3>(portLabelField);
        if (scalarField && labelField)
        {
            McWatch watch;

            if (portSelection.getValue() == 0)
            {
                computeAverageValuesInsideLabelFastAlgorithm(scalarField, labelField);
            }
            else
            {
                computeAverageValuesAroundAllLabels(scalarField, labelField);
            }

            theMsg->printf("Time for average computation: %f seconds", watch.stop());
        }
    }

    if (portVisualize.wasHit(0))
    {
        HxUniformLabelField3* labelField = hxconnection_cast<HxUniformLabelField3>(portLabelField);
        if (labelField)
        {
            const float x = portVisualizeCoords.getValue(0);
            const float y = portVisualizeCoords.getValue(1);
            const float z = portVisualizeCoords.getValue(2);
            visualizeNeighborhoodInsideLabel(labelField, x, y, z);
        }
    }
}

void
HxAverageValueInNeighborhood::computeAverageValuesInsideLabel(HxUniformScalarField3* scalarField, HxUniformScalarField3* labelField)
{
    const McDim3l& dims = labelField->lattice().getDims();
    HxUniformScalarField3* resultField = new HxUniformScalarField3(dims, McPrimType::MC_FLOAT);

    int numberOfVoxels;
    float totalValue;

#ifdef _OPENMP
    int numThreads = theSettingsMgr->getPreferences().maxNumberOfComputeThreads;
    if (numThreads < 1)
        numThreads = 1;
    omp_set_num_threads(numThreads);
#pragma omp parallel for
#endif
    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float label;
                labelField->lattice().eval(i, j, k, &label);
                int labelInt = (int)label;
                if (labelInt != 0)
                {
                    float valueInNeighborHood = getNeighborhoodValue(scalarField, labelField, i, j, k, numberOfVoxels, totalValue);
                    resultField->lattice().set(i, j, k, &valueInNeighborHood);
                }
                else
                {
                    float valueInNeighborHood = 0.0;
                    resultField->lattice().set(i, j, k, &valueInNeighborHood);
                }
            }
        }
    }

    const McBox3f& bbox = labelField->getBoundingBox();
    resultField->lattice().setBoundingBox(bbox);
    resultField->setLabel("mean");
    setResult(resultField);
}

void
HxAverageValueInNeighborhood::computeAverageValuesInsideLabelFastAlgorithm(HxUniformScalarField3* scalarField, HxUniformScalarField3* labelField)
{
    const McDim3l& dims = labelField->lattice().getDims();
    int numberOfVoxels;
    float totalValue;
    float valueInNeighborHood;
    float label;
    int labelInt;

    HxUniformScalarField3* resultField = new HxUniformScalarField3(dims, McPrimType::MC_FLOAT);

    McDArray<McDArray<float> > totalValuesArray;
    totalValuesArray.resize(dims[0]);
    for (int i = 0; i < dims[0]; i++)
    {
        totalValuesArray[i].resize(dims[1]);
    }

    McDArray<McDArray<int> > numberOfVoxelsArray;
    numberOfVoxelsArray.resize(dims[0]);
    for (int i = 0; i < dims[0]; i++)
    {
        numberOfVoxelsArray[i].resize(dims[1]);
    }

    labelField->lattice().eval(0, 0, 0, &label);
    labelInt = (int)label;
    if (labelInt != 0)
    {
        valueInNeighborHood = getNeighborhoodValue(scalarField, labelField, 0, 0, 0, numberOfVoxels, totalValue);
        numberOfVoxelsArray[0][0] = numberOfVoxels;
        totalValuesArray[0][0] = totalValue;
        resultField->lattice().set(0, 0, 0, &valueInNeighborHood);
    }
    else
    {
        valueInNeighborHood = 0.0;
        resultField->lattice().set(0, 0, 0, &valueInNeighborHood);
    }

    // Compute values for all i,0,0
    for (int i = 1; i < dims[0]; i++)
    {
        labelField->lattice().eval(i, 0, 0, &label);
        labelInt = (int)label;
        if (labelInt != 0)
        {
            labelField->lattice().eval(i - 1, 0, 0, &label);
            labelInt = (int)label;
            if (labelInt != 0)
            {
                // Compute value for i, 0, 0 with i-1, 0, 0 as reference
                valueInNeighborHood = getNeighborhoodValueUsingReferenceX(scalarField,
                                                                          labelField,
                                                                          i,
                                                                          0,
                                                                          0,
                                                                          numberOfVoxelsArray[i - 1][0],
                                                                          totalValuesArray[i - 1][0],
                                                                          numberOfVoxels,
                                                                          totalValue);
            }
            else
            {
                // Not possible to use i-1, 0, 0 as reference
                valueInNeighborHood = getNeighborhoodValue(scalarField, labelField, i, 0, 0, numberOfVoxels, totalValue);
            }
            // Set values
            numberOfVoxelsArray[i][0] = numberOfVoxels;
            totalValuesArray[i][0] = totalValue;
            resultField->lattice().set(i, 0, 0, &valueInNeighborHood);
        }
        else
        {
            valueInNeighborHood = 0.0;
            resultField->lattice().set(i, 0, 0, &valueInNeighborHood);
        }
    }

    // Compute values for all i,j,0 based on existing i,0,0 values
    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 1; j < dims[1]; j++)
        {
            labelField->lattice().eval(i, j, 0, &label);
            labelInt = (int)label;
            if (labelInt != 0)
            {
                labelField->lattice().eval(i, j - 1, 0, &label);
                labelInt = (int)label;
                if (labelInt != 0)
                {
                    // Compute value for i, j, 0 with i, j-1, 0 as reference
                    valueInNeighborHood = getNeighborhoodValueUsingReferenceY(scalarField,
                                                                              labelField,
                                                                              i,
                                                                              j,
                                                                              0,
                                                                              numberOfVoxelsArray[i][j - 1],
                                                                              totalValuesArray[i][j - 1],
                                                                              numberOfVoxels,
                                                                              totalValue);
                }
                else
                {
                    // Not possible to use i, j-1, 0 as reference
                    valueInNeighborHood = getNeighborhoodValue(scalarField, labelField, i, j, 0, numberOfVoxels, totalValue);
                }
                // Set values
                numberOfVoxelsArray[i][j] = numberOfVoxels;
                totalValuesArray[i][j] = totalValue;
                resultField->lattice().set(i, j, 0, &valueInNeighborHood);
            }
            else
            {
                valueInNeighborHood = 0.0;
                resultField->lattice().set(i, j, 0, &valueInNeighborHood);
            }
        }
    }

// Compute values for all i,j,0 based on existing i,j,0 values
#ifdef _OPENMP
    int numThreads = theSettingsMgr->getPreferences().maxNumberOfComputeThreads;
    if (numThreads < 1)
        numThreads = 1;
    omp_set_num_threads(numThreads);
#pragma omp parallel for
#endif
    for (int i = 0; i < dims[0]; i++)
    {
        int numberOfVoxels;
        float totalValue;
        float label;
        int labelInt;
        float valueInNeighborHood;
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 1; k < dims[2]; k++)
            {
                labelField->lattice().eval(i, j, k, &label);
                labelInt = (int)label;
                if (labelInt != 0)
                {
                    labelField->lattice().eval(i, j, k - 1, &label);
                    labelInt = (int)label;
                    if (labelInt != 0)
                    {
                        // Compute value for i, j, k with i, j, k-1 as reference
                        valueInNeighborHood = getNeighborhoodValueUsingReferenceZ(scalarField,
                                                                                  labelField,
                                                                                  i,
                                                                                  j,
                                                                                  k,
                                                                                  numberOfVoxelsArray[i][j],
                                                                                  totalValuesArray[i][j],
                                                                                  numberOfVoxels,
                                                                                  totalValue);
                    }
                    else
                    {
                        // Not possible to use i, j, k-1 as reference
                        valueInNeighborHood = getNeighborhoodValue(scalarField, labelField, i, j, k, numberOfVoxels, totalValue);
                    }
                    // Set values
                    numberOfVoxelsArray[i][j] = numberOfVoxels;
                    totalValuesArray[i][j] = totalValue;
                    resultField->lattice().set(i, j, k, &valueInNeighborHood);
                }
                else
                {
                    valueInNeighborHood = 0.0;
                    resultField->lattice().set(i, j, k, &valueInNeighborHood);
                }
            }
        }
    }

    const McBox3f& bbox = labelField->getBoundingBox();
    resultField->lattice().setBoundingBox(bbox);
    resultField->setLabel("mean");
    setResult(resultField);
}

float
HxAverageValueInNeighborhood::getNeighborhoodValueUsingReferenceZ(HxUniformScalarField3* scalarField,
                                                                  HxUniformScalarField3* labelField,
                                                                  const int x,
                                                                  const int y,
                                                                  const int z,
                                                                  const int numberVoxelsReference,
                                                                  const float totalValueReference,
                                                                  int& numberOfVoxels,
                                                                  float& totalValue)
{
    const int s = portSize.getValue(0);

    float totalValueOldArea = 0.0;
    int numberVoxelsOldArea = 0;
    float totalValueNewArea = 0.0;
    int numberVoxelsNewArea = 0;

    const McDim3l& dims = labelField->lattice().getDims();

    // Old area
    if (z - s - 1 >= 0)
    {
        for (int i = x - s; i <= x + s; i++)
        {
            for (int j = y - s; j <= y + s; j++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]))
                {
                    float labelFloat;
                    labelField->lattice().eval(i, j, z - s - 1, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsOldArea++;
                        float currentValue;
                        scalarField->lattice().eval(i, j, z - s - 1, &currentValue);
                        totalValueOldArea += currentValue;
                    }
                }
            }
        }
    }

    // New area
    if (z + s < dims[2])
    {
        for (int i = x - s; i <= x + s; i++)
        {
            for (int j = y - s; j <= y + s; j++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]))
                {
                    float labelFloat;
                    labelField->lattice().eval(i, j, z + s, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsNewArea++;
                        float currentValue;
                        scalarField->lattice().eval(i, j, z + s, &currentValue);
                        totalValueNewArea += currentValue;
                    }
                }
            }
        }
    }

    numberOfVoxels = numberVoxelsReference - numberVoxelsOldArea + numberVoxelsNewArea;
    totalValue = totalValueReference - totalValueOldArea + totalValueNewArea;

    if (numberOfVoxels == 0)
    {
        return 0;
    }

    return (totalValue) / (numberOfVoxels);
}

float
HxAverageValueInNeighborhood::getNeighborhoodValueUsingReferenceY(HxUniformScalarField3* scalarField,
                                                                  HxUniformScalarField3* labelField,
                                                                  const int x,
                                                                  const int y,
                                                                  const int z,
                                                                  const int numberVoxelsReference,
                                                                  const float totalValueReference,
                                                                  int& numberOfVoxels,
                                                                  float& totalValue)
{
    const int s = portSize.getValue(0);

    float totalValueOldArea = 0.0;
    int numberVoxelsOldArea = 0;
    float totalValueNewArea = 0.0;
    int numberVoxelsNewArea = 0;

    const McDim3l& dims = labelField->lattice().getDims();

    // Old area
    if (y - s - 1 >= 0)
    {
        for (int i = x - s; i <= x + s; i++)
        {
            for (int j = z - s; j <= z + s; j++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[2]))
                {
                    float labelFloat;
                    labelField->lattice().eval(i, y - s - 1, j, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsOldArea++;
                        float currentValue;
                        scalarField->lattice().eval(i, y - s - 1, j, &currentValue);
                        totalValueOldArea += currentValue;
                    }
                }
            }
        }
    }

    // New area
    if (y + s < dims[1])
    {
        for (int i = x - s; i <= x + s; i++)
        {
            for (int j = z - s; j <= z + s; j++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[2]))
                {
                    float labelFloat;
                    labelField->lattice().eval(i, y + s, j, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsNewArea++;
                        float currentValue;
                        scalarField->lattice().eval(i, y + s, j, &currentValue);
                        totalValueNewArea += currentValue;
                    }
                }
            }
        }
    }

    numberOfVoxels = numberVoxelsReference - numberVoxelsOldArea + numberVoxelsNewArea;
    totalValue = totalValueReference - totalValueOldArea + totalValueNewArea;

    if (numberOfVoxels == 0)
    {
        return 0;
    }

    return (totalValue) / (numberOfVoxels);
}

float
HxAverageValueInNeighborhood::getNeighborhoodValueUsingReferenceX(HxUniformScalarField3* scalarField,
                                                                  HxUniformScalarField3* labelField,
                                                                  const int x,
                                                                  const int y,
                                                                  const int z,
                                                                  const int numberVoxelsReference,
                                                                  const float totalValueReference,
                                                                  int& numberOfVoxels,
                                                                  float& totalValue)
{
    const int s = portSize.getValue(0);

    float totalValueOldArea = 0.0;
    int numberVoxelsOldArea = 0;
    float totalValueNewArea = 0.0;
    int numberVoxelsNewArea = 0;

    const McDim3l& dims = labelField->lattice().getDims();

    // Old area
    if (x - s - 1 >= 0)
    {
        for (int i = y - s; i <= y + s; i++)
        {
            for (int j = z - s; j <= z + s; j++)
            {
                if ((i >= 0) && (i < dims[1]) && (j >= 0) && (j < dims[2]))
                {
                    float labelFloat;
                    labelField->lattice().eval(x - s - 1, i, j, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsOldArea++;
                        float currentValue;
                        scalarField->lattice().eval(x - s - 1, i, j, &currentValue);
                        totalValueOldArea += currentValue;
                    }
                }
            }
        }
    }

    // New area
    if (x + s < dims[0])
    {
        for (int i = y - s; i <= y + s; i++)
        {
            for (int j = z - s; j <= z + s; j++)
            {
                if ((i >= 0) && (i < dims[1]) && (j >= 0) && (j < dims[2]))
                {
                    float labelFloat;
                    labelField->lattice().eval(x + s, i, j, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberVoxelsNewArea++;
                        float currentValue;
                        scalarField->lattice().eval(x + s, i, j, &currentValue);
                        totalValueNewArea += currentValue;
                    }
                }
            }
        }
    }

    numberOfVoxels = numberVoxelsReference - numberVoxelsOldArea + numberVoxelsNewArea;
    totalValue = totalValueReference - totalValueOldArea + totalValueNewArea;

    if (numberOfVoxels == 0)
    {
        return 0;
    }

    return (totalValue) / (numberOfVoxels);
}

float
HxAverageValueInNeighborhood::getNeighborhoodValue(HxUniformScalarField3* scalarField,
                                                   HxUniformScalarField3* labelField,
                                                   const int x,
                                                   const int y,
                                                   const int z,
                                                   int& numberOfVoxels,
                                                   float& totalValue)
{
    const int s = portSize.getValue(0);

    numberOfVoxels = 0;
    totalValue = 0;

    const McDim3l& dims = labelField->lattice().getDims();

    for (int k = z - s; k <= z + s; k++)
    {
        for (int j = y - s; j <= y + s; j++)
        {
            for (int i = x - s; i <= x + s; i++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]) && (k >= 0) && (k < dims[2]))
                {
                    float labelFloat;
                    labelField->lattice().eval(i, j, k, &labelFloat);
                    if (labelFloat > 0.0)
                    {
                        numberOfVoxels++;
                        float currentValue;
                        scalarField->lattice().eval(i, j, k, &currentValue);
                        totalValue += currentValue;
                    }
                }
            }
        }
    }

    if (numberOfVoxels == 0)
    {
        return 0;
    }

    return totalValue / numberOfVoxels;
}

void
HxAverageValueInNeighborhood::visualizeNeighborhoodInsideLabel(HxUniformLabelField3* labelField,
                                                               const int x,
                                                               const int y,
                                                               const int z)
{
    HxLoc3Regular* location = (HxLoc3Regular*)labelField->createLocation();
    location->set(McVec3f(x, y, z));
    location->getIx();

    const McDim3l& dims = labelField->labelLattice().getDims();
    const McBox3f& bbox = labelField->getBoundingBox();

    HxUniformLabelField3* outputLabelField = new HxUniformLabelField3(dims);
    outputLabelField->lattice().setBoundingBox(bbox);
    outputLabelField->setLabel("VisualizationLabelField");

    const int s = portSize.getValue(0);

    for (int k = location->getIz() - s; k <= location->getIz() + s; k++)
    {
        for (int j = location->getIy() - s; j <= location->getIy() + s; j++)
        {
            for (int i = location->getIx() - s; i <= location->getIx() + s; i++)
            {
                if ((i >= 0) && (i < dims[0]) && (j >= 0) && (j < dims[1]) && (k >= 0) && (k < dims[2]))
                {
                    float value = 1;
                    outputLabelField->labelLattice().set(i, j, k, &value);
                }
            }
        }
    }

    setResult(outputLabelField);
}

void
HxAverageValueInNeighborhood::computeAverageValuesAroundAllLabels(HxUniformScalarField3* scalarField,
                                                                  HxUniformLabelField3* labelField)
{
    float min = 0;
    float max = 0;
    labelField->getRange(min, max);

    const int size = portSize.getValue(0);

    HxUniformLabelField3* currentLabelField = labelField->duplicate();
    currentLabelField->setLabel("temporaryLabelField");
    theObjectPool->addObject(currentLabelField);

    HxQuant2GenericModule* dilationModule = (HxQuant2GenericModule*)HxResource::createObject("Dilation");
    HxConnection* connection = (HxConnection*)dilationModule->getPort("inputImage");
    connection->connect(currentLabelField);
    dilationModule->fire();
    HxPortIntTextN* portSize = (HxPortIntTextN*)dilationModule->getPort("size");
    portSize->setValue(size);
    dilationModule->fire();
    HxPortModuleSwitch* portType = (HxPortModuleSwitch*)dilationModule->getPort("Type");
    portType->setCurrentType(portType->getType(3));
    dilationModule->fire();

    HxSpreadSheet* spreadsheet = HxSpreadSheet::createInstance();
    spreadsheet->setLabel("AverageValuesAroundLabels-Spreadsheet");
    spreadsheet->addColumn("index", HxSpreadSheet::Column::INT);
    spreadsheet->addColumn("averageValue", HxSpreadSheet::Column::FLOAT);
    spreadsheet->addColumn("numberVoxels", HxSpreadSheet::Column::INT);
    HxSpreadSheet::Column* columnIndex = spreadsheet->column(0);
    HxSpreadSheet::Column* columnAverageValue = spreadsheet->column(1);
    HxSpreadSheet::Column* columnNumberVoxels = spreadsheet->column(2);

    theWorkArea->startWorking("Working..");

    // For all labels
    for (int i = 1; i <= max; i++)
    {
        float averageValueInsideLabel = 0;
        int numberOfVoxelsInsideLabel = 0;

        const McDim3l& dims = labelField->lattice().getDims();
// Create a labelfield that contains only label i (as value 1)
#ifdef _OPENMP
        int numThreads = theSettingsMgr->getPreferences().maxNumberOfComputeThreads;
        if (numThreads < 1)
            numThreads = 1;
        omp_set_num_threads(numThreads);
#pragma omp parallel for
#endif
        for (int l = 0; l < dims[0]; l++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int k = 0; k < dims[2]; k++)
                {
                    float label;
                    float scalarValue;
                    labelField->lattice().eval(l, j, k, &label);
                    scalarField->lattice().eval(l, j, k, &scalarValue);
                    if (label != i)
                    {
                        const float value = 0;
                        currentLabelField->lattice().set(l, j, k, &value);
                    }
                    else
                    {
                        const float value = 1;
#pragma omp critical
                        {
                            averageValueInsideLabel += scalarValue;
                            numberOfVoxelsInsideLabel++;
                        }
                        currentLabelField->lattice().set(l, j, k, &value);
                    }
                }
            }
        }

        // Compute average value around label i and add it to spreadsheet
        averageValueInsideLabel = averageValueInsideLabel / (float)numberOfVoxelsInsideLabel;
        int numberOfVoxelsLargerInsideAverage = 0;
        const float averageValueAroundLabel = computeAverageValueAroundLabel(scalarField, labelField, dilationModule, averageValueInsideLabel, numberOfVoxelsLargerInsideAverage);
        //theMsg->printf("value: %f", averageValueAroundLabel);
        theMsg->printf("number: %i", numberOfVoxelsLargerInsideAverage);
        spreadsheet->addRow();
        columnIndex->setValueInt(i - 1, i);
        columnAverageValue->setValue(i - 1, averageValueAroundLabel);
        columnNumberVoxels->setValue(i - 1, numberOfVoxelsLargerInsideAverage);

        theWorkArea->setProgressValue((float)i / float(max));
    }

    theWorkArea->stopWorking();
    setResult(spreadsheet);
}

float
HxAverageValueInNeighborhood::computeAverageValueAroundLabel(HxUniformScalarField3* scalarField,
                                                             HxUniformLabelField3* labelField,
                                                             HxQuant2GenericModule* dilationModule,
                                                             const float averageValueInsideLabel,
                                                             int& numberOfVoxels)
{
    HxPortDoIt* portDoItDilation = (HxPortDoIt*)dilationModule->getPort("doIt");
    portDoItDilation->setValue(0);
    dilationModule->fire();
    HxUniformLabelField3* dilationResult = dynamic_cast<HxUniformLabelField3*>(dilationModule->getResult());

    const float threshold = portThreshold.getValue();

    float averageValue = 0;
    numberOfVoxels = 0;
    const McDim3l& dims = labelField->lattice().getDims();

    for (int l = 0; l < dims[0]; l++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float labelInputImage;
                float labelDilationImage;
                float scalarValue;
                labelField->lattice().eval(l, j, k, &labelInputImage);
                dilationResult->lattice().eval(l, j, k, &labelDilationImage);
                scalarField->lattice().eval(l, j, k, &scalarValue);
                if ((labelInputImage == 0) && (labelDilationImage > 0) && (scalarValue > averageValueInsideLabel))
                {
                    averageValue += scalarValue;
                    numberOfVoxels++;
                }
            }
        }
    }
    if (numberOfVoxels > 0)
    {
        averageValue = averageValue / (float)numberOfVoxels;
    }
    else
    {
        averageValue = 0;
    }

    return averageValue;
}
