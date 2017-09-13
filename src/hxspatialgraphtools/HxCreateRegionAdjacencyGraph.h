#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortRadioBox.h>

#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <hxfield/HxUniformLabelField3.h>

/**
 * @brief The HxCreateRegionAdjacencyGraph class
 * Creates region adjacency graph (RAG) from given label field.
 * The RAG contains one vertex for each label in the label field
 * (placed in label centers; just labels with corresponding number of voxels > 0)
 * If two labels are neighbored, an edge is created.
 * Kind of neighborhood can be chosen between 6, 18 or 26 connectivity.
 */
class HxCreateRegionAdjacencyGraph : public HxCompModule
{
    HX_HEADER(HxCreateRegionAdjacencyGraph);

public:
    void compute();
    void update();

    HxPortRadioBox portNeighborhood;
    HxPortDoIt portDoIt;

private:
    /**
     * Type of neighborhood used by getLatticeNeighbors
     */
    enum Neighborhood
    {
        /// Voxels are neighbored if their faces are connected
        NEIGHBORHOOD_6,
        /// Voxels are neighbored if their edges are connected
        NEIGHBORHOOD_18,
        /// Voxels are neighbored if their corners are connected
        NEIGHBORHOOD_26,
    };

    /**
     * Returns directly connected neighbor lattice nodes of a given lattice node idx.
     * @param dims [in] dimensions of the lattice
     * @param idx [in] index of the node
     * @param isSeedIncluded [in] if true the input lattice node idx is included into the output neighbor nodes (at first position)
     * @param neighbors [out] nodes ids of lattice neighbors
     * @param neighborhood [in] type of lattice neighborhood
     */
    static void getLatticeNeighbors(const McDim3l& dims, const mcint64 idx, const bool isSeedIncluded, McDArray<mcint64>& neighbors, Neighborhood neighborhood = NEIGHBORHOOD_26);

    /**
     * Convert lattice position from (i,j,k) to idx.
     * There is no test whether input lattice position exists.
     * @param dims [in] dimensions of lattice
     * @param gridPos [in] triplet coordinates (i,j,k) of lattice node
     * @return lattice node index idx
     */
    static mcint64 gridPositionToIndex(const McDim3l& dims, const McVec3i& gridPos);

    /**
     * Convert lattice position from idx to (i,j,k).
     * There is no test whether input lattice position exists.
     * @param dims [in] dimensions of lattice
     * @param idx [in] lattice index idx of lattice node
     * @return lattice node triplet coordinates (i,j,k)
     */
    static McVec3i indexToGridPosition(const McDim3l& dims, const mcint64 idx);

    HxSpatialGraph* initOutput(const QString& inputLabelFieldName);

    HxSpatialGraph* computeRAG(HxUniformLabelField3* labelField, HxCreateRegionAdjacencyGraph::Neighborhood neighborhood);
    void initCentersAndNeighborhood(HxUniformLabelField3* labelField,
                                    const McDim3l& dims,
                                    const int nLabels,
                                    McDArray<McVec3i>& offsets,
                                    McDArray<McVec3f>& centers,
                                    McDArray<int>& numberVoxelsPerLabel,
                                    McBitfield& isNeighboredBitfield);
    void createGraphVertices(HxSpatialGraph* graph,
                             const int nLabels,
                             McDArray<McVec3f>& centers,
                             McDArray<int>& numberVoxelsPerLabel,
                             McDArray<int>& labelToVertexId);
    void createGraphEdges(HxSpatialGraph* graph,
                          const int nLabels,
                          McDArray<int>& labelToVertexId,
                          McBitfield& isNeighboredBitfield);
    void createEdge(HxSpatialGraph* graph, const int vertexId1, const int vertexId2);
    bool isIndexInDataset(const McDim3l& dims, McVec3i& index);
};
