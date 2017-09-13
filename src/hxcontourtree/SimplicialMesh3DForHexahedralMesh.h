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

class HXCONTOURTREE_API SimplicialMesh3DForHexahedralMesh : public SimplicialMesh
{
public:
    SimplicialMesh3DForHexahedralMesh(HxLattice3* lattice);
    ~SimplicialMesh3DForHexahedralMesh();
    /* Overloaded functions */
    void setThreshold(const double threshold);
    /// Initialized node IDs (mesh vertices with values larger than threshold),
    /// computes mapping from mesh IDs to node IDs
    /// and returns a sorted list of node IDs
    void getSortedListOfVertices(const bool computeSplitTree,
                                 McDArray<mculong>& sortedListOfVertices);
    /// Returns only neighbors with smaller or larger values than nodeIdx
    /// and which are above the threshold
    /// Uses m_neighborhood to determine the neighborhood
    /// Neighbors are mesh IDs
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
    void getNeighborsOfMeshVertex(const bool neighborsWithSmallerValue,
                                  const mculong vertexIdx,
                                  McDArray<mculong>& neighbors,
                                  McDArray<float>& values,
                                  const int typeOfNeighborhood = NEIGHBORHOOD_SIMPLICIALMESH) const;
    /// Returns all nodes, it does not matter whether the values are larger or smaller
    void get26NeighborHoodOfMeshVertex(const mculong vertexIdx,
                                       McDArray<mculong>& neighbors);
    /// Returns all nodes, it does not matter whether the values are larger or smaller
    void getAllNeighborsOfMeshVertex(const mculong vertexIdx,
                                     McDArray<mculong>& neighbors,
                                     const int typeOfNeighborHood);
    void getGridPosFromIdx(const mculong idx,
                           McVec3i& gridPos) const;
    mculong getIdxFromGridPos(const mculong x,
                              const mculong y,
                              const mculong z) const;

public:
    enum
    {
        NEIGHBORHOOD_6 = 0,
        NEIGHBORHOOD_18 = 1,
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
