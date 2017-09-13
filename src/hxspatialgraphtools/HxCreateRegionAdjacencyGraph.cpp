#include "HxCreateRegionAdjacencyGraph.h"

#include <hxcore/HxSettingsMgr.h>

#include <mclib/internal/McWatch.h>
#include <mclib/McBox3i.h>

#ifdef _OPENMP
#include <omp.h>
#endif

HX_INIT_CLASS(HxCreateRegionAdjacencyGraph, HxCompModule)

HxCreateRegionAdjacencyGraph::HxCreateRegionAdjacencyGraph()
    : HxCompModule(HxUniformLabelField3::getClassTypeId())
    , portNeighborhood(this, "neighborhood", tr("Neighborhood"), 3)
    , portDoIt(this, "action", tr("Action"))
{
    portNeighborhood.setLabel(0, "6");
    portNeighborhood.setLabel(1, "18");
    portNeighborhood.setLabel(2, "26");
}

HxCreateRegionAdjacencyGraph::~HxCreateRegionAdjacencyGraph()
{
}

void
HxCreateRegionAdjacencyGraph::update()
{
}

void
HxCreateRegionAdjacencyGraph::compute()
{
    if (portDoIt.wasHit())
    {
        HxUniformLabelField3* labelField = hxconnection_cast<HxUniformLabelField3>(portData);

        if (!labelField)
        {
            theMsg->printf("Attach label field");
            return;
        }

        HxSpatialGraph* outputGraph = 0;
        if (portNeighborhood.getValue() == 0)
        {
            outputGraph = computeRAG(labelField, HxCreateRegionAdjacencyGraph::NEIGHBORHOOD_6);
        }
        if (portNeighborhood.getValue() == 1)
        {
            outputGraph = computeRAG(labelField, HxCreateRegionAdjacencyGraph::NEIGHBORHOOD_18);
        }
        if (portNeighborhood.getValue() == 2)
        {
            outputGraph = computeRAG(labelField, HxCreateRegionAdjacencyGraph::NEIGHBORHOOD_26);
        }
        setResult(outputGraph);
    }
}

HxSpatialGraph*
HxCreateRegionAdjacencyGraph::initOutput(const QString& inputLabelFieldName)
{
    // Prepare output
    HxSpatialGraph* outputGraph = (HxSpatialGraph*)getResult();
    if (!outputGraph)
    {
        outputGraph = HxSpatialGraph::createInstance();
        outputGraph->composeLabel(inputLabelFieldName, "RAG");
        outputGraph->addAttribute("Label", HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1);
    }
    else
    {
        outputGraph->clear();
        outputGraph->addAttribute("Label", HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1);
    }
    return outputGraph;
}

HxSpatialGraph*
HxCreateRegionAdjacencyGraph::computeRAG(HxUniformLabelField3* labelField, HxCreateRegionAdjacencyGraph::Neighborhood neighborhood)
{
    McWatch watch;

    float minLabel, maxLabel;
    labelField->getRange(minLabel, maxLabel);
    const int nLabels = (int)maxLabel + 1;
    const McDim3l& dims = labelField->lattice().getDims();

    // Init neighborhood
    McDArray<McVec3i> offsets;
    offsets.remax(13);
    offsets.append(McVec3i(1, 0, 0));
    offsets.append(McVec3i(0, 1, 0));
    offsets.append(McVec3i(0, 0, 1));
    if (neighborhood == HxCreateRegionAdjacencyGraph::NEIGHBORHOOD_18)
    {
        offsets.append(McVec3i(1, 1, 0));
        offsets.append(McVec3i(1, -1, 0));
        offsets.append(McVec3i(1, 0, 1));
        offsets.append(McVec3i(0, 1, 1));
        offsets.append(McVec3i(-1, 0, 1));
        offsets.append(McVec3i(0, -1, 1));
    }
    if (neighborhood == HxCreateRegionAdjacencyGraph::NEIGHBORHOOD_26)
    {
        offsets.append(McVec3i(1, 1, 1));
        offsets.append(McVec3i(1, -1, 1));
        offsets.append(McVec3i(-1, 1, 1));
        offsets.append(McVec3i(-1, -1, 1));
    }

    HxSpatialGraph* outputGraph = initOutput(labelField->getLabel());

    McDArray<McVec3f> centers;
    McDArray<int> numberVoxelsPerLabel;
    McBitfield isNeighboredBitfield; // with background label 0
    initCentersAndNeighborhood(labelField, dims, nLabels, offsets, centers, numberVoxelsPerLabel, isNeighboredBitfield);

    McDArray<int> labelToVertexId; // with background label 0
    createGraphVertices(outputGraph, nLabels, centers, numberVoxelsPerLabel, labelToVertexId);

    createGraphEdges(outputGraph, nLabels, labelToVertexId, isNeighboredBitfield);

    theMsg->printf("Time for RAG computation: %f seconds", watch.stop());

    return outputGraph;
}

void
HxCreateRegionAdjacencyGraph::initCentersAndNeighborhood(HxUniformLabelField3* labelField,
                                                         const McDim3l& dims,
                                                         const int nLabels,
                                                         McDArray<McVec3i>& offsets,
                                                         McDArray<McVec3f>& centers,
                                                         McDArray<int>& numberVoxelsPerLabel,
                                                         McBitfield& isNeighboredBitfield)
{
    isNeighboredBitfield.resize(nLabels * nLabels);
    isNeighboredBitfield.unsetAll();

    centers.resize(nLabels);
    numberVoxelsPerLabel.resize(nLabels);
    centers.fill(McVec3f(0, 0, 0));
    numberVoxelsPerLabel.fill(0);

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
                McVec3i index(i, j, k);
                float label;
                labelField->lattice().eval(i, j, k, &label);
                int labelInt = (int)label;
                if (labelInt == 0)
                {
                    continue;
                }

                McVec3f pos = labelField->coords()->pos(index);
#ifdef _OPENMP
#pragma omp critical
#endif
                centers[labelInt] += pos;
#ifdef _OPENMP
#pragma omp critical
#endif
                numberVoxelsPerLabel[labelInt]++;

                for (int l = 0; l < offsets.sizeInt(); l++)
                {
                    McVec3i neighborIndex = index + offsets[l];
                    if (isIndexInDataset(dims, neighborIndex))
                    {
                        float neighborLabel;
                        labelField->lattice().eval(neighborIndex[0], neighborIndex[1], neighborIndex[2], &neighborLabel);
                        int neighborLabelInt = (int)neighborLabel;
                        if ((labelInt != neighborLabelInt) && (neighborLabelInt > 0))
                        {
                            isNeighboredBitfield.set(labelInt * nLabels + neighborLabelInt);
                            isNeighboredBitfield.set(neighborLabelInt * nLabels + labelInt);
                        }
                    }
                }
            }
        }
    }

    for (int i = 1; i < nLabels; i++)
    {
        if (numberVoxelsPerLabel[i] > 0)
        {
            centers[i] = centers[i] / numberVoxelsPerLabel[i];
        }
    }
}

void
HxCreateRegionAdjacencyGraph::createGraphVertices(HxSpatialGraph* graph,
                                                  const int nLabels,
                                                  McDArray<McVec3f>& centers,
                                                  McDArray<int>& numberVoxelsPerLabel,
                                                  McDArray<int>& labelToVertexId)
{
    labelToVertexId.resize(nLabels); // with background label 0
    labelToVertexId.fill(-1);
    EdgeVertexAttribute* labelAttribute = graph->findVertexAttribute("Label");
    int numberGraphVertices = 0;
    for (int i = 1; i < nLabels; i++)
    {
        if (numberVoxelsPerLabel[i] > 0)
        {
            graph->addVertex(centers[i]);
            labelAttribute->setIntDataAtIdx(numberGraphVertices, i);
            labelToVertexId[i] = numberGraphVertices;
            numberGraphVertices++;
        }
    }
}

void
HxCreateRegionAdjacencyGraph::createGraphEdges(HxSpatialGraph* graph,
                                               const int nLabels,
                                               McDArray<int>& labelToVertexId,
                                               McBitfield& isNeighboredBitfield)
{
    for (int i = 1; i < nLabels; i++)
    {
        for (int j = 1; j < nLabels; j++)
        {
            if (isNeighboredBitfield[i * nLabels + j])
            {
                createEdge(graph, labelToVertexId[i], labelToVertexId[j]);
            }
        }
    }
}

void
HxCreateRegionAdjacencyGraph::createEdge(HxSpatialGraph* graph, const int vertexId1, const int vertexId2)
{
    if (!graph->hasEdge(vertexId1, vertexId2) && !graph->hasEdge(vertexId2, vertexId1))
    {
        graph->addEdge(vertexId1, vertexId2);
    }
}

bool
HxCreateRegionAdjacencyGraph::isIndexInDataset(const McDim3l& dims, McVec3i& index)
{
    return ((index[0] > 0) && (index[1] > 0) && (index[2] > 0) &&
            (index[0] < dims[0]) && (index[1] < dims[1]) && (index[2] < dims[2]));
}

void
HxCreateRegionAdjacencyGraph::getLatticeNeighbors(const McDim3l& dims, const mcint64 idx, const bool isSeedIncluded, McDArray<mcint64>& neighbors, Neighborhood neighborhood)
{
    const McVec3i gridPos = HxCreateRegionAdjacencyGraph::indexToGridPosition(dims, idx);
    const McBox3i bbox(0, (int)dims.nx - 1, 0, (int)dims.ny - 1, 0, (int)dims.nz - 1);

    // Init neighborhood offsets
    // clang-format off
    static const McVec3i s_neighborhoodOffsets[26] = {
        // 6 neighborhood offsets
        McVec3i( 1,  0,  0),
        McVec3i( 0,  1,  0),
        McVec3i( 0,  0,  1),
        McVec3i(-1,  0,  0),
        McVec3i( 0, -1,  0),
        McVec3i( 0,  0, -1),

        // 18 neighborhood offsets (that are not already in 6 neighborhood)
        McVec3i( 1,  1,  0),
        McVec3i(-1, -1,  0),
        McVec3i( 0,  1,  1),
        McVec3i( 0, -1, -1),
        McVec3i( 1,  0,  1),
        McVec3i(-1,  0, -1),
        McVec3i( 1, -1,  0),
        McVec3i(-1,  1,  0),
        McVec3i( 0,  1, -1),
        McVec3i( 0, -1,  1),
        McVec3i( 1,  0, -1),
        McVec3i(-1,  0,  1),

        // 26 neighborhood offsets (that are not already in 6 or 18 neighborhood)
        McVec3i( 1,  1,  1),
        McVec3i( 1, -1,  1),
        McVec3i(-1,  1,  1),
        McVec3i(-1, -1,  1),
        McVec3i( 1,  1, -1),
        McVec3i( 1, -1, -1),
        McVec3i(-1,  1, -1),
        McVec3i(-1, -1, -1)
    };
    // clang-format on

    int numberNeighbors = 6;
    if (neighborhood == NEIGHBORHOOD_18)
    {
        numberNeighbors = 18;
    }
    if (neighborhood == NEIGHBORHOOD_26)
    {
        numberNeighbors = 26;
    }

    neighbors.clear();
    neighbors.remax(numberNeighbors + 1);
    if (isSeedIncluded)
    {
        neighbors.append(idx);
    }
    for (int i = 0; i < numberNeighbors; i++)
    {
        const McVec3i neighborGridPos = gridPos + s_neighborhoodOffsets[i];

        if (bbox.contains(neighborGridPos))
        {
            // neighborGridPos is not outside of dataset
            const mcint64 neighborLatticeIdx = HxCreateRegionAdjacencyGraph::gridPositionToIndex(dims, neighborGridPos);
            neighbors.append(neighborLatticeIdx);
        }
    }
}

mcint64
HxCreateRegionAdjacencyGraph::gridPositionToIndex(const McDim3l& dims, const McVec3i& gridPos)
{
    return (gridPos[2] * dims.ny + gridPos[1]) * dims.nx + gridPos[0];
}

McVec3i
HxCreateRegionAdjacencyGraph::indexToGridPosition(const McDim3l& dims, const mcint64 idx)
{
    const int z = idx / (dims[1] * dims[0]);
    const int rest = idx - z * (dims[1] * dims[0]);
    const int y = rest / dims[0];
    const int x = rest - y * dims[0];
    return McVec3i(x, y, z);
}
