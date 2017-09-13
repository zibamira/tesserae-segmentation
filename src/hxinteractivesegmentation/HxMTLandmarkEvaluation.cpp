#include "HxMTLandmarkEvaluation.h"
#include <hxspreadsheet/internal/HxSpreadSheet.h>
#include <mclib/McVec3i.h>
#include <mclib/McBitfield.h>
#include <queue>
#include <hxcore/HxMessage.h>

HX_INIT_CLASS(HxMTLandmarkEvaluation, HxCompModule);

HxMTLandmarkEvaluation::HxMTLandmarkEvaluation()
    : HxCompModule(HxUniformScalarField3::getClassTypeId())
    , portNeighborhood(this, "neighborhood", tr("Neighborhood"), 2)
    , portDoIt(this, "doIt", tr("Action"))
    , portLandmarks(this, "landmarks", tr("Landmarks"), HxLandmarkSet::getClassTypeId())
{
    portNeighborhood.setLabel(0, "Direct neighbors");
    portNeighborhood.setLabel(1, "BFS");

    portNeighborhood.setValue(1);
}

HxMTLandmarkEvaluation::~HxMTLandmarkEvaluation()
{
}

void
HxMTLandmarkEvaluation::compute()
{
    if (!portDoIt.wasHit())
    {
        return;
    }

    HxUniformScalarField3* labels = hxconnection_cast<HxUniformScalarField3>(portData);
    if (!labels)
        return;
    HxLandmarkSet* landmarks = hxconnection_cast<HxLandmarkSet>(portLandmarks);
    if (!landmarks)
        return;

    // Contains one landmark for each input landmark that leads to background (0) label region
    HxLandmarkSet* backgroundLandmarks = HxLandmarkSet::createInstance();
    // Contains one landmark for each region containing more than one input landmark
    HxLandmarkSet* clusterLandmarks = HxLandmarkSet::createInstance();

    McVec3f VOXEL_SIZE = labels->getVoxelSize();

    McArrayView<McVec3f> coords = landmarks->getCoords();
    const int NUM_LANDMARKS = landmarks->getCoords().size();
    float labelMin, labelMax;
    labels->lattice().computeRange(labelMin, labelMax);
    const int NUM_LABELS = (int)(labelMax - labelMin + 1);

    McDArray<McVec3f> labelToLandmarkPosition;
    labelToLandmarkPosition.resize(NUM_LABELS);

    McDArray<int> labelHit(NUM_LABELS);
    labelHit.fill(0);
    McVec3i gridPos;
    float curLabel;
    for (int i = 0; i < NUM_LANDMARKS; i++)
    {
        if (portNeighborhood.getValue() == 0)
            curLabel = getLabel(coords[i], *labels);
        else
            curLabel = getNearestForegroundLabelBFS(coords[i], labels);

        if (curLabel == 0)
        {
            backgroundLandmarks->appendLandmark(coords[i]);
        }

        labelHit[(int)curLabel]++;
        labelToLandmarkPosition[(int)curLabel] = coords[i];
    }

    HxSpreadSheet* result = HxSpreadSheet::createInstance();
    result->setTableName("Evaluation", 0);
    result->addColumn("Label", 1);
    result->addColumn("#hits", 1);
    result->addColumn("correct", 1);
    result->addColumn("over-segmented", 1);
    result->addColumn("2-cluster", 1);
    result->addColumn("3-cluster", 1);
    result->addColumn("4-cluster", 1);
    result->addColumn("5-cluster", 1);
    result->setNumRows(NUM_LABELS);

    int numberOfCorrectSegments = 0;
    int numberOfOversegmentations = 0;
    int numberOfClusters = 0;
    int numberOfWeightedClusters = 0;
    int numberOfBackgroundHits = labelHit[0];

    for (int i = 0; i < NUM_LABELS; i++)
    {
        result->columns[0].setValue(i, i);
        result->columns[1].setValue(i, labelHit[i]);

        int correct = labelHit[i] == 1;
        result->columns[2].setValue(i, correct);

        int overseg = labelHit[i] == 0;
        result->columns[3].setValue(i, overseg);

        int cluster2 = labelHit[i] == 2;
        result->columns[4].setValue(i, cluster2);
        int cluster3 = labelHit[i] == 3;
        result->columns[5].setValue(i, cluster3);
        int cluster4 = labelHit[i] == 4;
        result->columns[6].setValue(i, cluster4);
        int cluster5 = labelHit[i] == 5;
        result->columns[7].setValue(i, cluster5);

        // Exclude background from computations
        if (i > 0)
        {
            if (correct)
            {
                numberOfCorrectSegments++;
            }
            if (overseg)
            {
                numberOfOversegmentations++;
            }
            if (labelHit[i] > 1)
            {
                numberOfClusters++;
                numberOfWeightedClusters += labelHit[i] - 1;
            }

            if (overseg)
            {
                // output: label
                theMsg->printf("Oversegmentation: label %i", i);
            }
            if (cluster2 || cluster3 || cluster4 || cluster5)
            {
                clusterLandmarks->appendLandmark(labelToLandmarkPosition[i]);
            }
        }
    }

    // Compute precision / recall values
    // Count one true positive for one correct label (exactly one hit) and for one n-cluster
    // Count n-1 false negatives for one n-cluster and m false negatives for m background hits
    // Count one false positive for one oversegmented label
    // We separate between number of clusters and number of false negative
    // Assume you have a 3 cluster, then you would have to split the segment into 3 segments
    // So we have 2 false negatives for one cluster
    const double precision = (double)(numberOfCorrectSegments + numberOfClusters) /
                             (numberOfCorrectSegments + numberOfClusters + numberOfOversegmentations);
    const double recall = (double)(numberOfCorrectSegments + numberOfClusters) /
                          (numberOfCorrectSegments + numberOfClusters + numberOfWeightedClusters + numberOfBackgroundHits);

    // Create second table containing one line of summarized information
    result->addTable("Summary");
    result->addColumn("# correct labels", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# clusters", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# weighted clusters", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# background hits", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# true positives", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# false positives", HxSpreadSheet::Column::INT, 1);
    result->addColumn("# false negatives", HxSpreadSheet::Column::INT, 1);
    result->addColumn("Precision", HxSpreadSheet::Column::FLOAT, 1);
    result->addColumn("Recall", HxSpreadSheet::Column::FLOAT, 1);
    result->setNumRows(1, 1);
    result->column(0, 1)->setValue(0, numberOfCorrectSegments);
    result->column(1, 1)->setValue(0, numberOfClusters);
    result->column(2, 1)->setValue(0, numberOfWeightedClusters);
    result->column(3, 1)->setValue(0, numberOfBackgroundHits);
    result->column(4, 1)->setValue(0, numberOfCorrectSegments + numberOfClusters);
    result->column(5, 1)->setValue(0, numberOfOversegmentations);
    result->column(6, 1)->setValue(0, numberOfWeightedClusters + numberOfBackgroundHits);
    result->column(7, 1)->setValue(0, precision);
    result->column(8, 1)->setValue(0, recall);

    theMsg->printf(
        "All background hits are taken as false negatives, background label does not count as oversegmented"
        ", correct or cluster! Number of correct labels: %i, number of clusters ignoring their size: %i"
        ", weighted clusters (n-cluster counts for n-1): %i"
        ", number of background hits: %i"
        ", number of true positives (number of correct segments + number of clusters): %i"
        ", number of false positives (number of oversegmentations): %i"
        ", number of false negatives (number weighted clusters + number of background hits): %i"
        ", precision: %f, recall: %f",
        numberOfCorrectSegments,
        numberOfClusters,
        numberOfWeightedClusters,
        numberOfBackgroundHits,
        numberOfCorrectSegments + numberOfClusters,
        numberOfOversegmentations,
        numberOfWeightedClusters + numberOfBackgroundHits,
        precision,
        recall);

    result->composeLabel(labels->getLabel(), "landmark-evaluation");

    setResult(0, result);

    backgroundLandmarks->setLabel("Background Landmarks");
    clusterLandmarks->setLabel("Cluster Landmarks");
    setResult(1, backgroundLandmarks);
    setResult(2, clusterLandmarks);
}

// Find nearest non-background label using bfs
float
HxMTLandmarkEvaluation::getNearestForegroundLabelBFS(const McVec3f& pos, HxUniformScalarField3* labels)
{
    const McDim3l& dims = labels->lattice().getDims();
    const int dimsSize = dims[0] * dims[1] * dims[2];

    HxLoc3Regular* location = (HxLoc3Regular*)labels->createLocation();
    location->set(pos);
    McBitfield bitfield(dimsSize);

    const int startingMeshId = location->getIz() * (dims[1] * dims[0]) + location->getIy() * dims[0] + location->getIx();

    std::queue<int> meshIdQueue;
    meshIdQueue.push(startingMeshId);
    bitfield.set(startingMeshId);

    float labelValue;

    while (!meshIdQueue.empty())
    {
        const int meshId = meshIdQueue.front();
        meshIdQueue.pop();
        labels->lattice().eval(meshId, &labelValue);
        if (labelValue != 0.0)
        {
            return labelValue;
        }

        int neighborId = meshId + 1;
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
        neighborId = meshId - 1;
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
        neighborId = meshId + dims[0];
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
        neighborId = meshId - dims[0];
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
        neighborId = meshId + dims[0] * dims[1];
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
        neighborId = meshId - dims[0] * dims[1];
        if ((neighborId < dimsSize) && (neighborId >= 0) && (!bitfield[neighborId]))
        {
            meshIdQueue.push(neighborId);
            bitfield.set(neighborId);
        }
    }

    return 0.0;
}

void
HxMTLandmarkEvaluation::getSurroundingGridPoints(const McVec3f& pos, const HxUniformScalarField3& field, McDArray<McVec3i>& gridPoints)
{
    const McVec3f& voxelSize = field.getVoxelSize();
    const McBox3f& bbox = field.getBoundingBox();

    McVec3i nearest;
    McVec3i low;
    McVec3i high;
    for (int i = 0; i < 3; i++)
    {
        float fractional = (pos[i] - bbox[2 * i] + 0.5f * voxelSize[i]) / voxelSize[i];
        nearest[i] = (int)fractional;
        low[i] = (int)floor(fractional);
        high[i] = (int)ceil(fractional);
    }
    // nearest voxel at pos0, 8 surrounding voxels afterwards
    gridPoints.remax(9);
    gridPoints.append(nearest);
    if (low[0] >= 0 && low[1] >= 0 && low[2] >= 0)
        gridPoints.append(McVec3i(low[0], low[1], low[2]));
    if (low[0] >= 0 && low[1] >= 0 && high[2] < field.lattice().getDims()[2])
        gridPoints.append(McVec3i(low[0], low[1], high[2]));
    if (low[0] >= 0 && high[1] < field.lattice().getDims()[1] && low[2] >= 0)
        gridPoints.append(McVec3i(low[0], high[1], low[2]));
    if (low[0] >= 0 && high[1] < field.lattice().getDims()[1] && high[2] < field.lattice().getDims()[2])
        gridPoints.append(McVec3i(low[0], high[1], high[2]));
    if (high[0] < field.lattice().getDims()[0] && low[1] >= 0 && low[2] >= 0)
        gridPoints.append(McVec3i(high[0], low[1], low[2]));
    if (high[0] < field.lattice().getDims()[0] && low[1] >= 0 && high[2] < field.lattice().getDims()[2])
        gridPoints.append(McVec3i(high[0], low[1], high[2]));
    if (high[0] < field.lattice().getDims()[0] && high[1] < field.lattice().getDims()[1] && low[2] >= 0)
        gridPoints.append(McVec3i(high[0], high[1], low[2]));
    if (high[0] < field.lattice().getDims()[0] && high[1] < field.lattice().getDims()[1] && high[2] < field.lattice().getDims()[2])
        gridPoints.append(McVec3i(high[0], high[1], high[2]));
}

float
HxMTLandmarkEvaluation::getLabel(const McVec3f& pos, HxUniformScalarField3& field)
{
    McDArray<McVec3i> gridPoints;
    getSurroundingGridPoints(pos, field, gridPoints);
    float label = field.evalReg(gridPoints[0][0], gridPoints[0][1], gridPoints[0][2]);

    // if nearest voxel is not background return its label
    if (label > 0)
        return label;

    // get labels of neighborhood
    McDArray<float> neighborhoodLabels;
    for (int i = 1; i < gridPoints.size(); i++)
        neighborhoodLabels.append(field.evalReg(gridPoints[i][0], gridPoints[i][1], gridPoints[i][2]));

    // if there is exactly one non-background label in neighborhood, return it
    for (int i = 0; i < neighborhoodLabels.size(); i++)
        if (label == 0 && neighborhoodLabels[i] > 0)
            label = neighborhoodLabels[i];
        else if (label > 0 && neighborhoodLabels[i] > 0 && label != neighborhoodLabels[i])
            return 0;

    return label;
}
