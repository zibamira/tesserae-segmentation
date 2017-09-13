#ifndef AUGMENTED_CONTOUR_TREE_H
#define AUGMENTED_CONTOUR_TREE_H

#include <mclib/internal/McSmallArray.h>
#include <mclib/McDArray.h>
#include <mclib/McVec3.h>

#include "SimplicialMesh.h"
#include "SetUnionDataStructure.h"

#include "api.h"

class McBitfield;

class AugmentedContourTreeNode
{
public:
    friend class AugmentedContourTree;

    AugmentedContourTreeNode();

private:
    bool isUpperLeaf() const;

private:
    McDArray< mclong >               m_parentNodes;
    McDArray< mclong >               m_childNodes;
};

class HXCONTOURTREE_API AugmentedContourTree
{
public:
    AugmentedContourTree(SimplicialMesh * mesh);

    void computeJoinTree();
    void computeSplitTree();
    void computeContourTree();

    mclong getNumNodes() const;
    void getEdges(McDArray< McSmallArray< mclong, 2 > > & edges) const;
    void getComponents(McDArray< mclong > & components) const;
    void getCriticalPointNodesAndEdges(McDArray< mculong > & maxima,
                                       McDArray< mculong > & minima,
                                       McDArray< mculong > & saddles,
                                       McDArray< McSmallArray< mclong, 2 > > & edges) const;
    void getEdgesBetweenCriticalPoints(McDArray< McSmallArray< mclong, 2 > > & edges) const;
    void getContourTree(McDArray< mclong >  & vertexIds,
                        McDArray< McSmallArray< mclong, 2 > > & edges) const;
    void getSortedListOfMeshVertices(McDArray< mculong > & vertexIds);

private:
    McDArray< AugmentedContourTreeNode > & getNodes();
    void computeJoinSplitTree(const bool computeSplitTree);
    void mergeJoinAndSplitTree(AugmentedContourTree & joinTree,
                               AugmentedContourTree & splitTree);
    void unionFindMerge(const mculong nodeIdx, const mculong newComponent);
    mculong findNodeComponent(const mculong nodeIdx);
    void addEdgeToTree(const mculong nodeIdx1, const mculong nodeIdx2);
    void deleteNode(const mclong nodeId);
    void findComponent(const mclong         startNode,
                       McBitfield         & hasBeenProcessed,
                       McDArray< mclong > & components) const;
    void findComponentElements(mclong             & node,
                               McBitfield         & hasBeenProcessed,
                               McDArray< mclong > & components,
                               McDArray< mclong > & elements) const;
    void getCriticalPointNodesAndEdges(mculong               nodeIdx,
                                       mculong               nodeIdx2,
                                       McBitfield          & hasBeenProcessed,
                                       McDArray< mculong > & maxima,
                                       McDArray< mculong > & saddles,
                                       McDArray< McSmallArray< mclong, 2 > > & edges) const;
    mclong handleLeafNode(const AugmentedContourTreeNode & joinTreeNode,
                          const AugmentedContourTreeNode & splitTreeNode,
                          const mclong                     nodeId);

private:
    SimplicialMesh                          * m_mesh;
    SetUnionDataStructure                     m_setUnion;
    McDArray< AugmentedContourTreeNode >      m_nodes;
    McDArray< mculong >                       m_lowestNodeOfComponent;
    McDArray< mculong >                       m_sortedNodeIdx;
};

#endif
