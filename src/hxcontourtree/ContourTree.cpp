#include <mclib/internal/McAssert.h>
#include <mclib/McBitfield.h>

#include "ContourTree.h"
#include "AugmentedContourTree.h"

ContourTreeNode::ContourTreeNode()
{
}

ContourTreeNode::ContourTreeNode(const mclong meshVertexId)
{
    m_meshVertexId = meshVertexId;
}

ContourTree::ContourTree(const AugmentedContourTree & augContourTree)
{
    McDArray< mclong > mapAugmentedTreeNodesToNodes(augContourTree.getNumNodes());
    mapAugmentedTreeNodesToNodes.fill(-1);

    McDArray< mclong >  meshVertexIds;
    McDArray< McSmallArray< mclong, 2 > > edges;
    augContourTree.getContourTree(meshVertexIds, edges);

    for ( mclong i=0; i<meshVertexIds.size(); ++i )
        m_nodes.append(ContourTreeNode(meshVertexIds[i]));

    for ( mclong i=0; i<edges.size(); ++i )
        addEdge(edges[i][0], edges[i][1]);
}

mculong
ContourTree::getNumNodes() const
{
    return m_nodes.size();
}

void
ContourTree::addEdge(
    const mclong node1Id,
    const mclong node2Id)
{   
    m_nodes[node1Id].m_parentNodes.append(node2Id);
    m_nodes[node2Id].m_childNodes.append(node1Id);
}

void
ContourTree::getNodesAndEdges(
    McDArray< mculong > & maxima,
    McDArray< mculong > & minima,
    McDArray< mculong > & saddles,
    McDArray< McVec2i > & edges) const
{
    maxima.clear();
    minima.clear();
    saddles.clear();
    edges.clear();

    for ( mclong i=0; i<m_nodes.size(); ++i )
    {
        if (    m_nodes[i].m_parentNodes.size() > 1
             || m_nodes[i].m_childNodes.size() > 1 )
            saddles.append(i);
        else if ( m_nodes[i].m_parentNodes.size() == 0 )
            maxima.append(i);
        else if ( m_nodes[i].m_childNodes.size()  == 0 )
            minima.append(i);

        for ( mclong j=0; j<m_nodes[i].m_parentNodes.size(); ++j )
        {
            const mclong parentNode = m_nodes[i].m_parentNodes[j];
            edges.append( McVec2i(i, parentNode) );
        }
    }
}

void
ContourTree::getMaximaAndSaddles(
    McDArray< mculong > & maxima,
    McDArray< mculong > & saddles) const
{
    maxima.clear();
    saddles.clear();

    for ( mclong i=0; i<m_nodes.size(); ++i )
    {
        if ( m_nodes[i].m_parentNodes.size() == 0 )
            maxima.append(i);

        if ( m_nodes[i].m_parentNodes.size() > 1 )
            saddles.append(i);
    }
}

void
ContourTree::getMeshVertexIds(
    McDArray< mclong > & meshVertexIds) const
{
    const mculong numNodes = m_nodes.size();
    meshVertexIds.resize(numNodes);
    
    for ( mculong i=0; i<numNodes; ++i )
        meshVertexIds[i] = m_nodes[i].m_meshVertexId;
}

mclong
ContourTree::getMeshVertexId(const mclong & nodeId) const
{
    return m_nodes[nodeId].m_meshVertexId;
}
