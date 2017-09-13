#ifndef CONTOUR_TREE_H
#define CONTOUR_TREE_H

#include <mclib/McDArray.h>
#include <mclib/McVec2i.h>
#include <mclib/McHandable.h>

#include "api.h"

class AugmentedContourTree;

class ContourTreeNode
{
public:
    friend class ContourTree;

    ContourTreeNode();
    ContourTreeNode(const mclong meshVertexId);

private:
    mclong                           m_meshVertexId;
    McDArray< mclong >               m_parentNodes;
    McDArray< mclong >               m_childNodes;
};

class HXCONTOURTREE_API ContourTree : public McHandable
{
public:
    ContourTree(const AugmentedContourTree & augContourTree);

    mculong getNumNodes() const;

    void getNodesAndEdges(McDArray< mculong > & maxima,
                          McDArray< mculong > & minima,
                          McDArray< mculong > & saddles,
                          McDArray< McVec2i > & edges) const;

    void getMaximaAndSaddles(McDArray< mculong > & maxima,
                             McDArray< mculong > & saddles) const;

    void getMeshVertexIds(McDArray< mclong > & meshVertexIds) const;
    mclong getMeshVertexId(const mclong & nodeId) const;

private:
    void   addEdge(const mclong node1Id, const mclong node2Id);

private:
    McDArray< ContourTreeNode >      m_nodes;
};

#endif
