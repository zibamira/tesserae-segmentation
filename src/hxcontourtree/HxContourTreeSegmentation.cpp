#include <mclib/internal/McSorter.h>
#include <mclib/internal/McComparators.h>
#include <mclib/internal/McWatch.h>
#include <mclib/McTypedPointer.h>
#include <mclib/McRawData.h>
#include <mclib/McHistogram.h>

#include <hxcore/internal/HxWorkArea.h>

#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformLabelField3.h>

#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <hxplot/PzMarkerline.h>
#include <hxplot/PzCurve.h>

#include <hxcontourtree/AugmentedContourTree.h>

#include "HxContourTreeSegmentation.h"
#include "SimplicialMesh3DForHexahedralMesh.h"
#include "SetUnionDataStructure.h"

namespace
{
    template <typename T>
    int
    operator<(const McVec2<T>& a, const McVec2<T>& b)
    {
        if (a.x != b.x)
            return (a.x < b.x);
        return (a.y < b.y);
    }

    template <typename T>
    int
    operator>(const McVec2<T>& a, const McVec2<T>& b)
    {
        if (a.x != b.x)
            return (a.x > b.x);
        return (a.y > b.y);
    }
}

enum PersOptions
{
    PERS_RANGE,
    PERS_HISTO
};

enum PersModes
{
    PERS_ABSOLUTE,
    PERS_RELATIVE,
    PERS_ADAPTIVE
};

HX_INIT_CLASS(HxContourTreeSegmentation, HxCompModule);

HxContourTreeSegmentation::HxContourTreeSegmentation()
    : HxCompModule(HxLattice3::getClassTypeId())
    , portAutoSetValues(this, "auto-initialization", tr("Auto Initialization"))
    , portThreshold(this, "threshold", tr("Threshold"))
    , portPersistenceMode(this, "persistenceMode", tr("Persistence Mode"), 3)
    , portPersistenceOptions(this, "persistenceOptions", tr("Persistence Options"), 2)
    , portPersistenceValue(this, "persistenceValue", tr("Persistence Value"))
    , portMinimalSegmentSize(this, "minimalSegmentSize", tr("Minimal Segment Size"))
    , portSortNeighborsBy(this, "sortNeighborsBy", tr("Sort Neighbors By"), 2)
    , portFastSegmentation(this, "fastSegmentation", tr("Fast Segmentation"), 1)
    , portCalculatePlot(this, "calculatePlot", tr("Calculate Plot"), 1)
    , portDoIt(this, "doIt", tr("Do It"))
{
    m_plot = new PzEasyPlot("Number of segments");

    portAutoSetValues.setValue(0, 1);

    portPersistenceMode.setLabel(PERS_ABSOLUTE, "absolute");
    portPersistenceMode.setLabel(PERS_RELATIVE, "relative");
    portPersistenceMode.setLabel(PERS_ADAPTIVE, "adaptive");

    portPersistenceOptions.setLabel(PERS_RANGE, "range");
    portPersistenceOptions.setLabel(PERS_HISTO, "histogram");

    portSortNeighborsBy.setLabel(0, "value");
    portSortNeighborsBy.setLabel(1, "majority vote");
    portSortNeighborsBy.setValue(1);

    portCalculatePlot.setLabel("Plot");
    portCalculatePlot.setLabel(0, "Show");

    portMinimalSegmentSize.setMinMax(0, 1000);
    portMinimalSegmentSize.setValue(50);

    portFastSegmentation.setLabel("Fast Segmentation");

    m_outputLabelField = 0;
    m_relativePersistenceValue = 0.05f;
    m_recomputeHistogram = true;
    m_plotParametersChanged = true;
    m_needFastSegmentationInitialization = true;

    m_inputField = 0;
    m_contourTree = 0;
    m_mesh = 0;
}

HxContourTreeSegmentation::~HxContourTreeSegmentation()
{
}

void
HxContourTreeSegmentation::update()
{
    if (portData.isNew())
    {
        m_recomputeHistogram = true;
        m_plotParametersChanged = true;

        m_inputField = hxconnection_cast<HxUniformScalarField3>(portData);
        if (m_inputField)
        {
            m_mesh = new SimplicialMesh3DForHexahedralMesh(&(m_inputField->lattice()));
            m_mesh->setNeighborhood(SimplicialMesh3DForHexahedralMesh::NEIGHBORHOOD_26);

            recomputeRange();
            setScalarValueWidth();

            if (portAutoSetValues.getValue(0) == 1)
            {
                initThresholdPort();
                initPersistencePort();
            }
        }
    }

    if (!m_inputField)
        return;

    if (portAutoSetValues.isNew() && portAutoSetValues.getValue(0) == 1)
    {
        initThresholdPort();
        initPersistencePort();
    }

    if (portPersistenceMode.isNew())
    {
        setScalarValueWidth();
        initPersistencePort();
        m_plotParametersChanged = true;
    }

    if (portPersistenceOptions.isNew())
    {
        setScalarValueWidth();
        initPersistencePort();
        m_plotParametersChanged = true;
    }

    if (portPersistenceValue.isNew())
    {
        updateRelativePersistenceValue();
        m_plotParametersChanged = true;
    }

    if (portMinimalSegmentSize.isNew())
    {
        m_plotParametersChanged = true;
    }

    if (portThreshold.isNew())
    {
        m_plotParametersChanged = true;
        m_needFastSegmentationInitialization = true;
    }

    if (portSortNeighborsBy.isNew())
    {
        m_plotParametersChanged = true;
    }
}

void
HxContourTreeSegmentation::initThresholdPort()
{
    portThreshold.setMinMax(m_minValue, m_maxValue);
    portThreshold.setValue((m_minValue + m_maxValue) / 2.f);
}

void
HxContourTreeSegmentation::initPersistencePort()
{
    if (portPersistenceMode.getValue() == PERS_ABSOLUTE)
    {
        portPersistenceValue.setMinMax(0, m_scalarValueWidth);
        portPersistenceValue.setValue(m_scalarValueWidth * m_relativePersistenceValue);
    }
    else
    {
        portPersistenceValue.setMinMax(0.f, 1.f);
        portPersistenceValue.setValue(m_relativePersistenceValue);
    }
}

void
HxContourTreeSegmentation::updateRelativePersistenceValue()
{
    float persistenceValue = portPersistenceValue.getValue();
    if (portPersistenceMode.getValue() == PERS_ABSOLUTE)
        m_relativePersistenceValue = persistenceValue / m_scalarValueWidth;
    else
        m_relativePersistenceValue = persistenceValue;
}

float
HxContourTreeSegmentation::getPersistenceValue()
{
    if ((portPersistenceMode.getValue() == PERS_ABSOLUTE) ||
        (portPersistenceMode.getValue() == PERS_ADAPTIVE))
    {
        return portPersistenceValue.getValue();
    }
    else
    {
        const float persValue = portPersistenceValue.getValue();
        return (persValue * m_scalarValueWidth);
    }

    return 0.f;
}

void
HxContourTreeSegmentation::compute()
{
    if (portDoIt.wasHit())
    {
        if (!m_inputField)
        {
            return;
        }
        if (portFastSegmentation.getValue(0))
        {
            const float persistenceValue = getPersistenceValue();
            const mculong minNumberVoxels = (mculong)portMinimalSegmentSize.getValue();

            // Fast segmentation: iterate over contour tree merge points
            if (m_needFastSegmentationInitialization)
            {
                initOutputLabelField();
                initFastContourTreeSegmentation();
                fastContourTreeSegmentation(persistenceValue, minNumberVoxels);
            }
            else
            {
                initOutputLabelField();
                fastContourTreeSegmentation(persistenceValue, minNumberVoxels);
            }
        }
        else
        {
            // Normal segmentation: iterate over whole data set
            float persistenceValue = getPersistenceValue();
            const mculong minNumberVoxels = (mculong)portMinimalSegmentSize.getValue();
            initOutputLabelField();
            computeSegmentation(persistenceValue, minNumberVoxels);
        }
    }

    if (portCalculatePlot.wasHit())
    {
        if (m_plotParametersChanged)
        {
            if (!m_inputField)
            {
                return;
            }
            if (portFastSegmentation.getValue(0))
            {
                initOutputLabelField();
                initFastContourTreeSegmentation();
                calculatePlotForFastSegmentation();
            }
            else
            {
                calculatePlot();
            }
        }
        m_plot->show();
    }
}

void
HxContourTreeSegmentation::initOutputLabelField()
{
    theWorkArea->startWorking("Initializing data ...");

    if (!m_inputField)
        return;

    const McDim3l& dims = m_inputField->lattice().getDims();
    const McBox3f& bbox = m_inputField->getBoundingBox();

    const McPrimType outputFieldPrimType = McPrimType::MC_INT32;
    m_outputLabelField = dynamic_cast<HxUniformLabelField3*>(getResult());

    if (m_outputLabelField == 0)
    {
        m_outputLabelField = new HxUniformLabelField3(dims, outputFieldPrimType);
    }
    else
    {
        m_outputLabelField->lattice().setPrimType(outputFieldPrimType);
        m_outputLabelField->lattice().resize(dims);
    }

    m_outputLabelField->lattice().setBoundingBox(bbox);

    const mculong size = (mculong)dims[0] * (mculong)dims[1] * (mculong)dims[2];

    mcint32* data = (mcint32*)m_outputLabelField->lattice().dataPtr();

    const mcuint32 background = 0;
    for (mculong i = 0; i < size; ++i)
        data[i] = background;

    m_outputLabelField->composeLabel(m_inputField->getLabel(), "segmentation");

    theWorkArea->stopWorking();
}

// This plot computation is only working for the standard contour tree segmentation
// not the fast version
void
HxContourTreeSegmentation::calculatePlot()
{
    if (!m_inputField)
        return;

    McWatch watch;

    theWorkArea->startWorking("Initializing data ...");

    m_mesh->setThreshold(portThreshold.getValue());

    McDArray<mculong> sortedNodeIdx;
    const bool increasingOrder = false;
    m_mesh->getSortedListOfVertices(increasingOrder, sortedNodeIdx);

    SetUnionDataStructure setUnion;
    setUnion.setNumElements(sortedNodeIdx.size());

    McDArray<float> maxValuesOfComponent(sortedNodeIdx.size());
    McDArray<mculong> sizeOfComponent(sortedNodeIdx.size());
    McDArray<mculong> neighbors;
    McDArray<mculong> sortedNeighborIds;
    McDArray<float> values;

    McDArray<float> persistanceValues;
    McDArray<float> numberOfSegments;
    const float maxPersistanceValue = getPersistenceValue();
    const mclong numberOfPlotPoints = 10000;
    persistanceValues.resize(numberOfPlotPoints);
    numberOfSegments.resize(numberOfPlotPoints);
    float binSize = maxPersistanceValue / (numberOfPlotPoints - 1);
    for (int i = 0; i < numberOfPlotPoints; i++)
    {
        persistanceValues[i] = i * binSize;
        numberOfSegments[i] = 0;
    }

    mclong numberOfMaxima = 0;

    const mculong minNumberVoxels = (mculong)portMinimalSegmentSize.getValue();

    bool localRelative = (portPersistenceMode.getValue() == PERS_ADAPTIVE);
    float totalMin, totalMax;
    m_inputField->lattice().computeRange(totalMin, totalMax);

    theWorkArea->setProgressInfo("Calculating plot ...");

    for (mclong i = 0; i < sortedNodeIdx.size(); ++i)
    {
        if (i % 100 == 0)
        {
            theWorkArea->setProgressValue((float)i / (float)sortedNodeIdx.size());

            if (theWorkArea->wasInterrupted())
                break;
        }

        const mculong nodeIdx = sortedNodeIdx[i];
        const mculong component = nodeIdx;
        setUnion.setSetIdOfElement(nodeIdx, component);

        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        float value;
        m_inputField->lattice().eval(meshVertexId, &value);
        maxValuesOfComponent[nodeIdx] = value;
        sizeOfComponent[nodeIdx] = 1;

        neighbors.clear();
        values.clear();
        m_mesh->getNeighborsOfMeshVertex(increasingOrder, meshVertexId, neighbors, values);

        sortNeighbors(neighbors, values, sortedNeighborIds, setUnion);

        if (sortedNeighborIds.size() == 0)
        {
            numberOfMaxima++;
        }

        McDArray<mclong> ids;
        McDArray<float> values;
        McDArray<mculong> sizes;
        McDArray<bool> mergedBecauseOfSize;
        int numberOfMerges = 0;
        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);
            const mclong component = setUnion.findSetId(neighbor);

            if (ids.findSorted(component, &mcStandardCompare) == -1)
            {
                ids.append(component);
                values.append(maxValuesOfComponent[component]);
                if (j == 0)
                {
                    sizes.append(sizeOfComponent[component] + 1);
                }
                else
                {
                    sizes.append(sizeOfComponent[component]);
                }
                mergedBecauseOfSize.append(false);
                ids.sort(&mcStandardCompare);
            }
        }
        // Sort sizes (and values and components in the same way)
        for (int j = 0; j < sizes.sizeInt(); j++)
        {
            for (int k = j + 1; k < sizes.sizeInt(); k++)
            {
                if (sizes[k] < sizes[j])
                {
                    values.swap(j, k);
                    sizes.swap(j, k);
                }
            }
        }
        // Deal with minimal segment size
        // Largest component cannot merge into another
        for (mclong j = 0; j < sizes.size() - 1; ++j)
        {
            if (sizes[j] < minNumberVoxels)
            {
                mergedBecauseOfSize[j] = true;
                numberOfMerges++;
                for (mclong k = 0; k < numberOfPlotPoints; k++)
                {
                    numberOfSegments[k] = numberOfSegments[k] + 1;
                }
            }
        }
        // Sort values (and sizes and components in the same way)
        for (int j = 0; j < values.sizeInt(); j++)
        {
            for (int k = j + 1; k < values.sizeInt(); k++)
            {
                if (values[k] < values[j])
                {
                    values.swap(j, k);
                    sizes.swap(j, k);
                    mergedBecauseOfSize.swap(j, k);
                }
            }
        }
        // Deal with persistence value
        // Segments that are already merged using the segment size are not taken into account
        for (mclong j = 0; j < sizes.size() - 1; ++j)
        {
            if ((!mergedBecauseOfSize[j]) && (numberOfMerges < sizes.sizeInt() - 1))
            {
                float persistenceComponent = 0;
                if (localRelative)
                {
                    persistenceComponent = (values[j] - value) / (values[j] - totalMin);
                }
                else
                {
                    persistenceComponent = (values[j] - value);
                }
                for (mclong k = ((int)(persistenceComponent / binSize)) + 1; k < numberOfPlotPoints; k++)
                {
                    numberOfSegments[k] = numberOfSegments[k] + 1;
                }
                numberOfMerges++;
            }
        }

        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);
            const mclong component1 = setUnion.findSetId(nodeIdx);
            const mclong component2 = setUnion.findSetId(neighbor);

            if (component1 == component2)
                continue;

            const float maximumValue =
                (maxValuesOfComponent[component1] > maxValuesOfComponent[component2] ? maxValuesOfComponent[component1] : maxValuesOfComponent[component2]);

            setUnion.mergeSetsOfElements(nodeIdx, neighbor);
            maxValuesOfComponent[component1] = maximumValue;
            maxValuesOfComponent[component2] = maximumValue;
            const mculong sizeOfComponents = sizeOfComponent[component2] + sizeOfComponent[component1];
            sizeOfComponent[component1] = sizeOfComponents;
            sizeOfComponent[component2] = sizeOfComponents;
        }
    }

    for (mclong i = 0; i < numberOfPlotPoints; ++i)
    {
        numberOfSegments[i] = numberOfMaxima - numberOfSegments[i];
    }

    theWorkArea->stopWorking();

    theMsg->printf("plot calculation took %f seconds", watch.stop());

    updatePlot(numberOfPlotPoints, persistanceValues, numberOfSegments, maxPersistanceValue);

    m_plotParametersChanged = false;
}

void
HxContourTreeSegmentation::calculatePlotForFastSegmentation()
{
    McDArray<float> persistanceValues;
    McDArray<float> numberOfSegments;
    const float maxPersistanceValue = getPersistenceValue();
    const mculong minNumberVoxels = (mculong)portMinimalSegmentSize.getValue();

    const mclong numberOfPlotPoints = 200;
    persistanceValues.resize(numberOfPlotPoints);
    numberOfSegments.resize(numberOfPlotPoints);
    float binSize = maxPersistanceValue / (numberOfPlotPoints - 1);
    for (int i = 0; i < numberOfPlotPoints; i++)
    {
        persistanceValues[i] = i * binSize;
        numberOfSegments[i] = 0;
    }
    for (int i = 0; i < numberOfPlotPoints; i++)
    {
        const int numberSegments = fastContourTreeSegmentation(persistanceValues[i], minNumberVoxels);
        numberOfSegments[i] = numberSegments;
    }

    updatePlot(numberOfPlotPoints, persistanceValues, numberOfSegments, maxPersistanceValue);
}

// Returns number of segments
int
HxContourTreeSegmentation::computeSegmentation(const float persistenceValue, const mculong minNumberVoxels, AugmentedContourTree* augmentedJoinTree)
{
    if (!m_inputField)
        return -1;

    McWatch watch;

    theWorkArea->startWorking("Initializing data ...");

    m_mesh->setThreshold(portThreshold.getValue());

    McDArray<mculong> sortedNodeIdx;
    const bool increasingOrder = false;
    if (augmentedJoinTree)
    {
        augmentedJoinTree->getSortedListOfMeshVertices(sortedNodeIdx);
    }
    else
    {
        m_mesh->getSortedListOfVertices(increasingOrder, sortedNodeIdx);
    }

    SetUnionDataStructure setUnion;
    setUnion.setNumElements(sortedNodeIdx.size());

    McDArray<float> maxValuesOfComponent(sortedNodeIdx.size());
    McDArray<mculong> maxIdOfComponent(sortedNodeIdx.size());
    McDArray<mculong> sizeOfComponent(sortedNodeIdx.size());

    McDArray<mculong> neighbors;
    McDArray<mculong> sortedNeighborIds;
    McDArray<float> values;

    mcint32* data = (mcint32*)m_outputLabelField->lattice().dataPtr();

    theWorkArea->setProgressInfo("Propagating contours ...");

    bool localRelative = (portPersistenceMode.getValue() == PERS_ADAPTIVE);
    float totalMin, totalMax;
    m_inputField->lattice().computeRange(totalMin, totalMax);

    int numberSegments = 0;

    for (mclong i = 0; i < sortedNodeIdx.size(); ++i)
    {
        if (i % 100 == 0)
        {
            theWorkArea->setProgressValue((float)i / (float)sortedNodeIdx.size());

            if (theWorkArea->wasInterrupted())
                break;
        }

        const mculong nodeIdx = sortedNodeIdx[i];
        const mculong component = nodeIdx;
        setUnion.setSetIdOfElement(nodeIdx, component);

        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        float value;
        m_inputField->lattice().eval(meshVertexId, &value);
        maxValuesOfComponent[nodeIdx] = value;
        maxIdOfComponent[nodeIdx] = nodeIdx;
        sizeOfComponent[nodeIdx] = 1;

        neighbors.clear();
        values.clear();
        m_mesh->getNeighborsOfMeshVertex(increasingOrder, meshVertexId, neighbors, values);

        sortNeighbors(neighbors, values, sortedNeighborIds, setUnion);

        if (sortedNeighborIds.sizeInt() == 0)
        {
            numberSegments++;
        }

        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);
            const mclong component1 = setUnion.findSetId(nodeIdx);
            const mclong component2 = setUnion.findSetId(neighbor);

            if (component1 == component2)
                continue;

            float persistenceComponent1 = 0;
            float persistenceComponent2 = 0;

            if (localRelative)
            {
                persistenceComponent1 = (maxValuesOfComponent[component1] - value) / (maxValuesOfComponent[component1] - totalMin);
                persistenceComponent2 = (maxValuesOfComponent[component2] - value) / (maxValuesOfComponent[component2] - totalMin);
            }
            else
            {
                persistenceComponent1 = (maxValuesOfComponent[component1] - value);
                persistenceComponent2 = (maxValuesOfComponent[component2] - value);
            }

            if ((((maxValuesOfComponent[component1] - value) >= 0) && ((maxValuesOfComponent[component2] - value) >= 0)) && (((persistenceComponent1 <= persistenceValue) || (sizeOfComponent[component1] < minNumberVoxels)) || ((persistenceComponent2 <= persistenceValue) || (sizeOfComponent[component2] < minNumberVoxels))))
            {
                if (j != 0)
                {
                    numberSegments--;
                }

                setUnion.mergeSetsOfElements(nodeIdx, neighbor);
                if (maxValuesOfComponent[component1] >= maxValuesOfComponent[component2])
                {
                    maxValuesOfComponent[component2] = maxValuesOfComponent[component1];
                    maxIdOfComponent[component2] = maxIdOfComponent[component1];
                }
                else
                {
                    maxValuesOfComponent[component1] = maxValuesOfComponent[component2];
                    maxIdOfComponent[component1] = maxIdOfComponent[component2];
                }
                mculong sizeOfComponents = sizeOfComponent[component2] + sizeOfComponent[component1];
                sizeOfComponent[component2] = sizeOfComponents;
                sizeOfComponent[component1] = sizeOfComponents;
            }
        }

        const mcuint32 label = (mcuint32)(setUnion.findSetId(nodeIdx) + 1);
        data[meshVertexId] = label;
    }

    theWorkArea->stopWorking();

    relabelOutputLabelField(setUnion);

    m_setUnionFinestSegmentation = setUnion;
    m_maxValuesOfComponentFinestSegmentation = maxValuesOfComponent;
    m_sizeOfComponentFinestSegmentation = sizeOfComponent;

    theMsg->printf("segmentation took %f seconds", watch.stop());

    setResult(m_outputLabelField);

    return numberSegments;
}

void
HxContourTreeSegmentation::sortNeighbors(
    const McDArray<mculong>& neighbors,
    const McDArray<float>& values,
    McDArray<mculong>& sortedNeighborIds,
    SetUnionDataStructure& setUnion)
{
    if (portSortNeighborsBy.getValue() == 0)
        sortNeighborsAccordingToValues(values, sortedNeighborIds);
    else
        sortNeighborsAccordingToMajorityVote(neighbors, sortedNeighborIds, setUnion);
}

void
HxContourTreeSegmentation::sortNeighborsAccordingToValues(
    const McDArray<float>& values,
    McDArray<mculong>& sortedNeighborIds)
{
    sortedNeighborIds.resize(values.size());

    if (values.size() == 0)
        return;

    for (mclong i = 0; i < sortedNeighborIds.size(); ++i)
        sortedNeighborIds[i] = i;

    McIndexByValueDescendingComparator<float> comp(&values[0]);
    sort(&sortedNeighborIds[0],
         sortedNeighborIds.size(),
         comp);
}

void
HxContourTreeSegmentation::sortNeighborsAccordingToMajorityVote(
    const McDArray<mculong>& neighbors,
    McDArray<mculong>& sortedNeighborIds,
    SetUnionDataStructure& setUnion)
{
    sortedNeighborIds.resize(neighbors.size());

    if (neighbors.size() == 0)
        return;

    McDArray<mcuint32> mapToLabel;
    McDArray<int> labelCounter;
    McDArray<int> labelIds;

    for (int i = 0; i < neighbors.size(); ++i)
    {
        const mcuint32 label = (mcuint32)(setUnion.findSetId(m_mesh->getNodeIdx(neighbors[i])) + 1);

        if (label > 0)
        {
            int arrayId = -1;
            for (int j = 0; j < mapToLabel.size(); ++j)
            {
                if (mapToLabel[j] == label)
                {
                    arrayId = j;
                    ++labelCounter[j];
                    break;
                }
            }

            if (arrayId == -1)
            {
                arrayId = mapToLabel.size();
                mapToLabel.append(label);
                labelCounter.append(0);
            }

            labelIds.append(arrayId);
        }
        else
        {
            labelIds.append(-1);
        }
    }

    McDArray<mculong> numNeighborsWithEqualLabel(neighbors.size());
    for (mclong i = 0; i < numNeighborsWithEqualLabel.size(); ++i)
    {
        sortedNeighborIds[i] = i;

        if (labelIds[i] == -1)
            numNeighborsWithEqualLabel[i] = 0;
        else
            numNeighborsWithEqualLabel[i] = labelCounter[labelIds[i]];
    }

    McIndexByValueDescendingComparator<mculong> comp(&numNeighborsWithEqualLabel[0]);
    sort(&sortedNeighborIds[0],
         sortedNeighborIds.size(),
         comp);
}

void
HxContourTreeSegmentation::relabelOutputLabelField(
    SetUnionDataStructure& setUnion)
{
    const McDim3l& dims = m_outputLabelField->lattice().getDims();

    mculong size = dims[0] * dims[1] * dims[2];

    mcint32* data = (mcint32*)m_outputLabelField->lattice().dataPtr();

    McDArray<mclong> oldLabelToNewLabel(setUnion.getNumElements());
    oldLabelToNewLabel.fill(-1);

    // There is already the background label that will not be changed;
    // hence, numLabels equals 1.
    mculong numLabels = 1;
    for (mculong i = 0; i < size; ++i)
    {
        const mcuint32 label =
            evalOutputLabelField(m_outputLabelField->lattice(), i);

        if (label > 0)
        {
            const mculong componentLabel = setUnion.findSetId(label - 1);

            if (oldLabelToNewLabel[componentLabel] < 0)
            {
                oldLabelToNewLabel[componentLabel] = numLabels;
                ++numLabels;
            }

            const mculong newLabel = oldLabelToNewLabel[componentLabel];
            data[i] = newLabel;
        }
    }

    m_outputLabelField->touchMinMax();
}

mcuint32
HxContourTreeSegmentation::evalOutputLabelField(
    HxLattice3& lattice,
    mculong meshVertexId)
{
    const mcint32* data = (mcint32*)lattice.dataPtr();

    return data[meshVertexId];
}

void
HxContourTreeSegmentation::recomputeRange()
{
    if (!m_inputField)
        return;

    m_minValue = 0.f;
    m_maxValue = 0.f;
    m_inputField->lattice().computeRange(m_minValue, m_maxValue);
}

void
HxContourTreeSegmentation::setScalarValueWidth()
{
    if (portPersistenceOptions.getValue() == PERS_RANGE)
    {
        m_scalarValueWidth = m_maxValue - m_minValue;
    }
    else
    {
        if (m_recomputeHistogram)
            recomputeHistogram();
        m_scalarValueWidth = m_maxValueHistogram - m_minValueHistogram;
    }
}

void
HxContourTreeSegmentation::recomputeHistogram()
{
    if (!m_inputField)
        return;

    theWorkArea->startWorkingNoStop("Computing histogram ...");

    const McDim3l& dims = m_inputField->lattice().getDims();
    const mculong size = (mculong)dims[0] * (mculong)dims[1] * (mculong)dims[2];

    const int numBins = 256;

    McHistogram histogram(m_minValue, m_maxValue, numBins);

    McRawData::getHistogram(McTypedPointer(m_inputField->lattice().dataPtr(),
                                           m_inputField->lattice().primType()),
                            size,
                            histogram,
                            false,
                            0.0,
                            theWorkArea);

    getRange(histogram, m_minValueHistogram, m_maxValueHistogram, 0.05f, 0.95f);

    theMsg->printf("histogram range: [%.4f, %.4f]", m_minValueHistogram, m_maxValueHistogram);

    theWorkArea->stopWorking();

    m_recomputeHistogram = false;
}

void
HxContourTreeSegmentation::updatePlot(
    const mclong numberOfPlotPoints,
    McDArray<float>& persistanceValues,
    McDArray<float>& numberOfSegments,
    const float maxPersistanceValue)
{
    m_plot->autoUpdate(0);
    PzCurve* curve = m_plot->putData("Number of segments", numberOfPlotPoints, persistanceValues.dataPtr(), numberOfSegments.dataPtr());
    PzMarkerline* markerLine = m_plot->setMarkerline("line", maxPersistanceValue / 4);
    markerLine->setProbe(curve);
    markerLine->setAttr("useprobe", 1);
    markerLine->setAttr("probemarkertype", 2);
    markerLine->setAttr("showlabel", 1);
    markerLine->setAttr("labelalign2", 2);
    m_plot->update();
}

void
HxContourTreeSegmentation::initFastContourTreeSegmentation()
{
    McWatch watch;

    if (!m_inputField)
        return;

    // Initialize mesh
    m_mesh->setThreshold(portThreshold.getValue());

    // Initialize contour tree
    AugmentedContourTree augmentedContourTree(m_mesh);
    augmentedContourTree.computeJoinTree();
    m_contourTree = new ContourTree(augmentedContourTree);

    // Initialize finest segmentation
    m_numberOfSegmentsFinestSegmentation = computeSegmentation(0.0F, 0, &augmentedContourTree);
    m_needFastSegmentationInitialization = false;

    theMsg->printf("Fast contour tree segmentation initialization took %f seconds", watch.stop());
}

int
HxContourTreeSegmentation::fastContourTreeSegmentation(const float persistenceValue, const mculong minNumberVoxels)
{
    if (!m_inputField)
        return -1;

    McDArray<mculong> maxima;
    McDArray<mculong> saddles;
    m_contourTree->getMaximaAndSaddles(maxima, saddles);

    // Sort saddles
    McDArray<McVec2d> tuples;
    for (mclong i = 0; i < saddles.size(); i++)
    {
        const mculong contourTreeNodeId = saddles[i];
        const mclong meshVertexId = m_contourTree->getMeshVertexId(contourTreeNodeId);
        float value;
        m_inputField->lattice().eval(meshVertexId, &value);

        tuples.append(McVec2d(double(value), double(meshVertexId)));
    }
    tuples.sort(mcStandardCompare);

    McDArray<mculong> neighbors;
    McDArray<float> values;
    McDArray<mculong> sortedNeighborIds;

    SetUnionDataStructure setUnion = m_setUnionFinestSegmentation;
    McDArray<float> maxValuesOfComponent = m_maxValuesOfComponentFinestSegmentation;
    McDArray<mculong> sizeOfComponent = m_sizeOfComponentFinestSegmentation;

    bool localRelative = (portPersistenceMode.getValue() == PERS_ADAPTIVE);
    float totalMin, totalMax;
    m_inputField->lattice().computeRange(totalMin, totalMax);

    int numberOfSegments = m_numberOfSegmentsFinestSegmentation;

    // Iterate over all merge nodes in join tree
    // From high mesh values to low mesh values
    for (mclong i = saddles.size() - 1; i >= 0; i--)
    {
        const mclong meshVertexId = tuples[i][1];
        const mclong nodeIdx = m_mesh->getNodeIdx(meshVertexId);
        float value;
        m_inputField->lattice().eval(meshVertexId, &value);

        // Look at larger neighbors, sort them according to number of neighbors in current segmentation
        neighbors.clear();
        values.clear();
        m_mesh->getNeighborsOfMeshVertex(false, meshVertexId, neighbors, values);
        sortNeighbors(neighbors, values, sortedNeighborIds, setUnion);

        // Merge according to persistance value
        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);
            const mclong component1 = setUnion.findSetId(nodeIdx);
            const mclong component2 = setUnion.findSetId(neighbor);

            if (component1 == component2)
                continue;

            float persistenceComponent1 = 0;
            float persistenceComponent2 = 0;

            if (localRelative)
            {
                persistenceComponent1 = (maxValuesOfComponent[component1] - value) / (maxValuesOfComponent[component1] - totalMin);
                persistenceComponent2 = (maxValuesOfComponent[component2] - value) / (maxValuesOfComponent[component2] - totalMin);
            }
            else
            {
                persistenceComponent1 = (maxValuesOfComponent[component1] - value);
                persistenceComponent2 = (maxValuesOfComponent[component2] - value);
            }

            if ((((maxValuesOfComponent[component1] - value) >= 0) && ((maxValuesOfComponent[component2] - value) >= 0)) && (((persistenceComponent1 <= persistenceValue) || (sizeOfComponent[component1] < minNumberVoxels)) || ((persistenceComponent2 <= persistenceValue) || (sizeOfComponent[component2] < minNumberVoxels))))
            {
                setUnion.mergeSetsOfElements(nodeIdx, neighbor);
                if (maxValuesOfComponent[component1] >= maxValuesOfComponent[component2])
                {
                    maxValuesOfComponent[component2] = maxValuesOfComponent[component1];
                }
                else
                {
                    maxValuesOfComponent[component1] = maxValuesOfComponent[component2];
                }
                mculong sizeOfComponents = sizeOfComponent[component2] + sizeOfComponent[component1];
                sizeOfComponent[component2] = sizeOfComponents;
                sizeOfComponent[component1] = sizeOfComponents;

                numberOfSegments--;
            }
        }
    }

    mcint32* data = (mcint32*)m_outputLabelField->lattice().dataPtr();
    const McDim3l& dims = m_outputLabelField->lattice().getDims();
    mculong size = dims[0] * dims[1] * dims[2];
    for (mculong i = 0; i < size; ++i)
    {
        mclong nodeId = m_mesh->getNodeIdx(i);
        if (nodeId != -1)
        {
            const mcuint32 label = (mcuint32)(setUnion.findSetId(nodeId) + 1);
            data[i] = label;
        }
    }

    relabelOutputLabelField(setUnion);

    setResult(m_outputLabelField);

    return numberOfSegments;
}

void
HxContourTreeSegmentation::getRange(
    McHistogram& histogram,
    float& vMin,
    float& vMax,
    const float lowerBound,
    const float upperBound)
{
    const double min = histogram.getMin();
    const double max = histogram.getMax();

    if (lowerBound == 0.f && upperBound == 1.f)
    {
        vMin = static_cast<float>(min);
        vMax = static_cast<float>(max);
        return;
    }

    const mcuint64 overallCounts = histogram.getOverallCounts();

    mcuint64 minCount = McRoundEven(float(overallCounts) * lowerBound);
    mcuint64 maxCount = McRoundEven(float(overallCounts) * upperBound);

    const int numBins = histogram.getNumBins();

    int idxMax = numBins - 1;
    mcuint64 total = overallCounts;
    for (int i = numBins - 1; i >= 0; i--)
    {
        total -= histogram.getCount(i);
        if (total < maxCount)
            break;
        idxMax = i;
    }

    int idxMin = 0;
    total = 0;
    for (int i = 0; i < numBins; ++i)
    {
        total += histogram.getCount(i);
        if (total > minCount)
            break;
        idxMin = i;
    }

    const double binSize = fabs(max - min) / numBins;
    vMin = static_cast<float>(min + binSize * idxMin);
    vMax = static_cast<float>(min + binSize * idxMax);
}
