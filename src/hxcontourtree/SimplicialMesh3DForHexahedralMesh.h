#ifndef SIMPLICIAL_MESH_3D_FOR_HEXAHEDRAL_MESH_H
#define SIMPLICIAL_MESH_3D_FOR_HEXAHEDRAL_MESH_H

#include <mclib/McHandle.h>
#include <mclib/McDArray.h>
#include <mclib/McPrimType.h>
#include <mclib/McVec2.h>
#include <mclib/McVec3i.h>

#include "SimplicialMesh.h"
#include "AugmentedContourTree.h"

#include "api.h"

class HxLattice3;

/**
 * @brief The SimplicialMesh3DForHexahedralMesh class
 * Class for a 3D (simplicial) mesh with possible threshold.
 * If a threshold is set (setThreshold() is called) the method getSortedListOfVertices() initializes a mapping
 * between values larger than this threshold (one nodeId for each) and all mesh nodes (called meshVertexIds).
 * getNeighborsOfNode() and getNeighborsOfMeshVertex() return all neighbors larger than this threshold as
 * meshVertexIds. This is not working if getSortedListOfVertices() was not called before.
 * If no threshold is used then meshVertexIds and nodeIds are the same.
 * Type of neighborhood is defined by m_neighborhood, which is initialized as simplicial mesh
 * but can also be a 6 or 26 neighborhood.
 */

class HXCONTOURTREE_API SimplicialMesh3DForHexahedralMesh : public SimplicialMesh
{
public:
    SimplicialMesh3DForHexahedralMesh(HxLattice3* lattice);
    ~SimplicialMesh3DForHexahedralMesh();
    /* Overloaded functions */
    void setThreshold(const double threshold);
    /// Initializes node IDs (mesh vertices with values larger than threshold),
    /// computes mapping from mesh IDs to node IDs
    /// and returns a sorted list of node IDs
    void getSortedListOfVertices(const bool computeSplitTree,
                                 McDArray<mculong>& sortedListOfVertices);
    /// Returns only neighbors with smaller or larger values than nodeIdx
    /// and which are above the threshold
    /// Uses m_neighborhood to determine the neighborhood
    /// Neighbors are mesh IDs
    /// Neighborhood determined by m_neighborhood
    void getNeighborsOfNode(const bool neighborsWithSmallerValue,
                            const mculong nodeIdx,
                            McDArray<mculong>& neighbors,
                            McDArray<float>& values) const;
    mclong getMeshVertexIdx(const mclong nodeIdx) const;
    mclong getNodeIdx(const mclong vertexIdx) const;
    void setNeighborhood(const int neighborhood);
    int getNeighborhood() const;

    /* Non-overloaded functions */
    void getGridPos(McDArray<McVec3i>& gridPos) const;
    void getGridPos(const McDArray<mclong>& vertexIds,
                    McDArray<McVec3i>& gridPos) const;
    /// Returns only neighbors with smaller or larger values than vertexIdx
    /// and which are above the threshold
    /// Neighborhood determined by m_neighborhood
    void getNeighborsOfMeshVertex(const bool neighborsWithSmallerValue,
                                  const mculong vertexIdx,
                                  McDArray<mculong>& neighbors,
                                  McDArray<float>& values) const;
    /// Returns all nodes, it does not matter whether the values are larger or smaller
    /// Neighbors are mesh IDs
    void get26NeighborHoodOfMeshVertexIgnoringThreshold(const mculong vertexIdx,
                                                        McDArray<mculong>& neighbors);

    void getGridPosFromIdx(const mculong idx,
                           McVec3i& gridPos) const;
    mculong getIdxFromGridPos(const mculong x,
                              const mculong y,
                              const mculong z) const;

public:
    enum
    {
        NEIGHBORHOOD_6 = 0,
        NEIGHBORHOOD_SIMPLICIALMESH = 2,
        NEIGHBORHOOD_26 = 3,
    };

private:
    void addNeighbor(const bool neighborsWithSmallerValue,
                     const McVec2d& vertex,
                     const McVec3i& gridPos,
                     const McVec3i& offSet,
                     McDArray<mculong>& neighbors,
                     McDArray<float>& values) const;
    void addNeighbor(const McVec3i& gridPos,
                     const McVec3i& offSet,
                     McDArray<mculong>& neighbors) const;
    bool isValidGridPos(const McVec3i& gridPos) const;

private:
    bool m_useThreshold;
    double m_threshold;
    HxLattice3* m_lattice;
    mculong m_dims[3];
    float m_bbox[6];
    mculong m_numNodes;
    McDArray<mclong> m_mapNodeToMeshVertex;
    McDArray<mclong> m_mapMeshVertexToNode;

    int m_neighborhood;
};

#endif
