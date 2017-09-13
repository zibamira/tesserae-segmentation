#include <mclib/internal/McAssert.h>
#include <mclib/McBitfield.h>
#include <mclib/internal/McListQueue.h>

#include "SimplicialMesh.h"
#include "AugmentedContourTree.h"

AugmentedContourTreeNode::AugmentedContourTreeNode()
{
}

mclong
AugmentedContourTree::getNumNodes() const
{
    return m_nodes.size();
}

AugmentedContourTree::AugmentedContourTree(
    SimplicialMesh * mesh)
{
    m_mesh = mesh;
}

void
AugmentedContourTree::computeJoinTree()
{
    computeJoinSplitTree(false);
}

void
AugmentedContourTree::computeSplitTree()
{
    computeJoinSplitTree(true);
}

void
AugmentedContourTree::computeJoinSplitTree(
    const bool computeSplitTree)
{
    mcassert(m_mesh);

    m_mesh->getSortedListOfVertices(computeSplitTree, m_sortedNodeIdx);

    m_lowestNodeOfComponent.resize(m_sortedNodeIdx.size());
    m_nodes.resize(m_sortedNodeIdx.size());
    m_setUnion.setNumElements(m_sortedNodeIdx.size());

    McDArray< mculong > neighbors;
    McDArray< float >   values;
    for ( mclong i=0; i<m_sortedNodeIdx.size(); ++i )
    {
        const mculong nodeIdx   = m_sortedNodeIdx[i];
        const mculong component = nodeIdx;
        m_lowestNodeOfComponent[component] = nodeIdx;
        m_setUnion.setSetIdOfElement(nodeIdx, component);

        neighbors.clear();
        m_mesh->getNeighborsOfNode(computeSplitTree, nodeIdx, neighbors, values);
        
        for ( mclong j=0; j<neighbors.size(); ++j )
        {
            const mclong neighbor          = m_mesh->getNodeIdx(neighbors[j]);
            const mclong component         = m_setUnion.findSetId(nodeIdx);
            const mclong neighborComponent = m_setUnion.findSetId(neighbor);
            if ( component == neighborComponent )
                continue;

            m_setUnion.mergeSetsOfElements(nodeIdx, neighbor);

            addEdgeToTree(nodeIdx, m_lowestNodeOfComponent[neighborComponent]);

            m_lowestNodeOfComponent[neighborComponent] = nodeIdx;
        }
    }
}

void
AugmentedContourTree::addEdgeToTree(
    const mculong nodeIdx1,
    const mculong nodeIdx2)
{
    m_nodes[nodeIdx1].m_parentNodes.append(nodeIdx2);
    m_nodes[nodeIdx2].m_childNodes.append(nodeIdx1);
}

void
AugmentedContourTree::deleteNode(
    const mclong nodeId)
{
    mcassert(   m_nodes[nodeId].m_parentNodes.size()<2
             && m_nodes[nodeId].m_childNodes.size()<2);

    mclong parentNodeId = -1;
    mclong childNodeId = -1;

    if ( m_nodes[nodeId].m_parentNodes.size() > 0 )
    {
        parentNodeId = m_nodes[nodeId].m_parentNodes[0];
        m_nodes[nodeId].m_parentNodes.remove(0);
    }

    if ( m_nodes[nodeId].m_childNodes.size() > 0 )
    {
        childNodeId = m_nodes[nodeId].m_childNodes[0];
        m_nodes[nodeId].m_childNodes.remove(0);
    }

    if ( parentNodeId != -1 && childNodeId != -1 )
    {
        for ( mclong i=0; i<m_nodes[parentNodeId].m_childNodes.size(); ++i )
        {
            if ( m_nodes[parentNodeId].m_childNodes[i] == nodeId )
            {
                m_nodes[parentNodeId].m_childNodes[i] = childNodeId;
                break;
            }
        }        
        for ( mclong i=0; i<m_nodes[childNodeId].m_parentNodes.size(); ++i )
        {
            if ( m_nodes[childNodeId].m_parentNodes[i] == nodeId )
            {
                m_nodes[childNodeId].m_parentNodes[i] = parentNodeId;
                return;
            }
        }
    }
    else if ( parentNodeId != -1 )
    {
        for ( mclong i=0; i<m_nodes[parentNodeId].m_childNodes.size(); ++i )
        {
            if ( m_nodes[parentNodeId].m_childNodes[i] == nodeId )
            {
                m_nodes[parentNodeId].m_childNodes.remove(i);
                return;
            }
        }
    }
    else if ( childNodeId != -1 )
    {
        for ( mclong i=0; i<m_nodes[childNodeId].m_parentNodes.size(); ++i )
        {
            if ( m_nodes[childNodeId].m_parentNodes[i] == nodeId )
            {
                m_nodes[childNodeId].m_parentNodes.remove(i);
                return;
            }
        }
    }
}

void
AugmentedContourTree::computeContourTree()
{
    AugmentedContourTree splitTree(m_mesh);
    splitTree.computeSplitTree();

    AugmentedContourTree joinTree(m_mesh);
    joinTree.computeJoinTree();

    m_mesh->getSortedListOfVertices(false, m_sortedNodeIdx);

    mergeJoinAndSplitTree(joinTree, splitTree);
}

bool
AugmentedContourTreeNode::isUpperLeaf() const
{
    if (    m_parentNodes.size() == 0
         && m_childNodes.size()  == 1 )
        return true;

    return false;
}

mclong
AugmentedContourTree::handleLeafNode(
    const AugmentedContourTreeNode & joinTreeNode,
    const AugmentedContourTreeNode & splitTreeNode,
    const mclong                     nodeId)
{
    mclong adjacentNodeId = -1;
    if ( joinTreeNode.isUpperLeaf() )
    {
        adjacentNodeId = joinTreeNode.m_childNodes[0];
        addEdgeToTree(adjacentNodeId, nodeId);
    }
    else if ( splitTreeNode.isUpperLeaf() )
    {
        adjacentNodeId = splitTreeNode.m_childNodes[0];
        addEdgeToTree(nodeId, adjacentNodeId);
    }

    return adjacentNodeId;
}

void
AugmentedContourTree::mergeJoinAndSplitTree(
    AugmentedContourTree & joinTree,
    AugmentedContourTree & splitTree)
{
    McDArray< AugmentedContourTreeNode > & joinTreeNodes  = joinTree.getNodes();
    McDArray< AugmentedContourTreeNode > & splitTreeNodes = splitTree.getNodes();

    mcrequire(joinTreeNodes.size() == splitTreeNodes.size());

    McListQueue<mclong> queue;

    int sizeOfQueue = 0;
    m_nodes.resize(joinTreeNodes.size());
    for ( mclong i=0; i<m_nodes.size(); ++i )
    {
        if ( joinTreeNodes[i].m_parentNodes.size() + splitTreeNodes[i].m_parentNodes.size() == 1 )
        {
            queue.put(i);
            ++sizeOfQueue;
        }
    }

    while ( sizeOfQueue > 1 )
    {
        const mclong nodeId = queue.get();
        --sizeOfQueue;

        const mclong adjacentNodeId =
            handleLeafNode(joinTreeNodes[nodeId], splitTreeNodes[nodeId], nodeId);

        if ( adjacentNodeId == -1 )
            continue;

        joinTree.deleteNode(nodeId);
        splitTree.deleteNode(nodeId);

        if (    (   joinTreeNodes[adjacentNodeId].isUpperLeaf()
                 || splitTreeNodes[adjacentNodeId].isUpperLeaf() )
             && joinTreeNodes[adjacentNodeId].m_parentNodes.size() + splitTreeNodes[adjacentNodeId].m_parentNodes.size() == 1 )
        {
            queue.put(adjacentNodeId);
            ++sizeOfQueue;
        }
    }
}

McDArray< AugmentedContourTreeNode > &
AugmentedContourTree::getNodes()
{
    return m_nodes;
}

void
AugmentedContourTree::getEdges(
    McDArray< McSmallArray< mclong, 2 > > & edges) const
{
    for ( mclong i=0; i<m_nodes.size(); ++i )
        for ( mclong j=0; j<m_nodes[i].m_parentNodes.size(); ++j )
        {
            McSmallArray< mclong, 2 > edge(2);
            edge[0] = i;
            edge[1] = m_nodes[i].m_parentNodes[j];
            edges.append( edge );
        }
}

void
AugmentedContourTree::getComponents(
    McDArray< mclong > & components) const
{
    components.resize(m_nodes.size());
    components.fill(-1);

    McBitfield hasBeenProcessed(m_nodes.size());
    hasBeenProcessed.unsetAll();

    for ( mclong i=m_sortedNodeIdx.size()-1; i>-1; --i )
    {
        const mclong node = m_sortedNodeIdx[i];
        if ( ! hasBeenProcessed[node] )
        {
            findComponent(node, hasBeenProcessed, components);
        }
    }
}

void
AugmentedContourTree::findComponent(
    const mclong         startNode,
    McBitfield         & hasBeenProcessed,
    McDArray< mclong > & components) const
{
    McDArray< mclong > elements;

    mclong node = startNode;
    findComponentElements(node, hasBeenProcessed, components, elements);

    for ( mclong i=0; i<elements.size(); ++i )
        components[elements[i]] = node;
}

void
AugmentedContourTree::findComponentElements(
    mclong             & node,
    McBitfield         & hasBeenProcessed,
    McDArray< mclong > & components,
    McDArray< mclong > & elements) const
{
    while ( m_nodes[node].m_parentNodes.size() == 1 )
    {
        if ( hasBeenProcessed[node] )
        {
            elements.clear();
            return;
        }

        hasBeenProcessed.set(node);
        
        elements.append(node);
        node = m_nodes[node].m_parentNodes[0];
    }

    if ( m_nodes[node].m_parentNodes.size() == 0 )
    {
        elements.append(node);
        return;
    }
    else if ( m_nodes[node].m_parentNodes.size() > 1 )
    {
        elements.clear();

        for ( mclong j=0; j<m_nodes[node].m_parentNodes.size(); ++j )
            findComponent(m_nodes[node].m_parentNodes[j], hasBeenProcessed, components);
    }
}

void 
AugmentedContourTree::getContourTree(
    McDArray< mclong >  & meshVertexIds,
    McDArray< McSmallArray< mclong, 2 > > & edges) const
{
    meshVertexIds.clear();
    edges.clear();

    McDArray< mculong > maxima;
    McDArray< mculong > minima;
    McDArray< mculong > saddles;

    getCriticalPointNodesAndEdges(maxima, minima, saddles, edges);

    const mclong numNodes = getNumNodes();
    McDArray< mclong > mapNodesToCriticalPoints(numNodes);
    mapNodesToCriticalPoints.fill(-1);

    mclong numCriticalPoints = 0;
    for ( mclong i=0; i<maxima.size(); ++i )
    {
        meshVertexIds.append(m_mesh->getMeshVertexIdx(maxima[i]));
        mapNodesToCriticalPoints[maxima[i]] = numCriticalPoints;
        ++numCriticalPoints;
    }

    for ( mclong i=0; i<saddles.size(); ++i )
    {
        meshVertexIds.append(m_mesh->getMeshVertexIdx(saddles[i]));
        mapNodesToCriticalPoints[saddles[i]] = numCriticalPoints;
        ++numCriticalPoints;
    }

    for ( mclong i=0; i<minima.size(); ++i )
    {
        meshVertexIds.append(m_mesh->getMeshVertexIdx(minima[i]));
        mapNodesToCriticalPoints[minima[i]] = numCriticalPoints;
        ++numCriticalPoints;
    }

    for ( mclong i=0; i<edges.size(); ++i )
    {
        edges[i][0] = mapNodesToCriticalPoints[edges[i][0]];
        edges[i][1] = mapNodesToCriticalPoints[edges[i][1]];
    }
}

void
AugmentedContourTree::getEdgesBetweenCriticalPoints(
    McDArray< McSmallArray< mclong, 2 > > & edges) const
{
    McDArray< mculong > maxima;
    McDArray< mculong > minima;
    McDArray< mculong > saddles;

    getCriticalPointNodesAndEdges(maxima, minima, saddles, edges);
}

void
AugmentedContourTree::getCriticalPointNodesAndEdges(
    McDArray< mculong > & maxima,
    McDArray< mculong > & minima,
    McDArray< mculong > & saddles,
    McDArray< McSmallArray< mclong, 2 > > & edges) const
{
    maxima.clear();
    minima.clear();
    saddles.clear();
    edges.clear();

    const mclong numNodes = m_sortedNodeIdx.size();

    McBitfield hasBeenProcessed(numNodes);
    hasBeenProcessed.unsetAll();

    for ( mclong i=numNodes-1; i>=0; --i )
    {
        const mculong nodeIdx = m_sortedNodeIdx[i];
        if ( ! hasBeenProcessed[nodeIdx] )
        {
            hasBeenProcessed.set(nodeIdx);

            if (    m_nodes[nodeIdx].m_parentNodes.size() > 1
                 || m_nodes[nodeIdx].m_childNodes.size() > 1 )
                saddles.append(nodeIdx);
            else if ( m_nodes[nodeIdx].m_parentNodes.size() == 0 )
                maxima.append(nodeIdx);
            else if ( m_nodes[nodeIdx].m_childNodes.size() == 0 )
                minima.append(nodeIdx);
            else
                printf("ERROR in AugmentedContourTree::getCriticalPointNodesAndEdges.\n");

            for ( mclong i=0; i<m_nodes[nodeIdx].m_parentNodes.size(); ++i )
            {
                const mculong nodeIdx2 = m_nodes[nodeIdx].m_parentNodes[i];
                getCriticalPointNodesAndEdges(nodeIdx,
                                              nodeIdx2,
                                              hasBeenProcessed,
                                              maxima,
                                              saddles,
                                              edges);
            }
        }
    }
}

void
AugmentedContourTree::getCriticalPointNodesAndEdges(
    mculong               nodeIdx,
    mculong               nodeIdx2,
    McBitfield          & hasBeenProcessed,
    McDArray< mculong > & maxima,
    McDArray< mculong > & saddles,
    McDArray< McSmallArray< mclong, 2 > > & edges) const
{
    McDArray< mculong > nodeIdxPairs;
    nodeIdxPairs.append(nodeIdx);
    nodeIdxPairs.append(nodeIdx2);

    for ( mclong firstNodeIdx = 0; firstNodeIdx < nodeIdxPairs.size(); firstNodeIdx += 2 )
    {
        nodeIdx  = nodeIdxPairs[firstNodeIdx];
        nodeIdx2 = nodeIdxPairs[firstNodeIdx+1];

        while (    m_nodes[nodeIdx2].m_parentNodes.size() == 1
                && m_nodes[nodeIdx2].m_childNodes.size()  <= 1 )
        {
            hasBeenProcessed.set(nodeIdx2);
            nodeIdx2 = m_nodes[nodeIdx2].m_parentNodes[0];
        }

        McSmallArray< mclong, 2 > edge(2);
        edge[0] = nodeIdx;
        edge[1] = nodeIdx2;
        edges.append( edge );

        if ( hasBeenProcessed[nodeIdx2] )
            continue;

        hasBeenProcessed.set(nodeIdx2);

        if (    m_nodes[nodeIdx2].m_parentNodes.size() == 0
             && m_nodes[nodeIdx2].m_childNodes.size()  <= 1 )
            maxima.append(nodeIdx2);
        else
            saddles.append(nodeIdx2);

        for ( mclong i=0; i<m_nodes[nodeIdx2].m_parentNodes.size(); ++i )
        {
            nodeIdxPairs.append(nodeIdx2);
            nodeIdxPairs.append(m_nodes[nodeIdx2].m_parentNodes[i]);
        }
    }
}

void
AugmentedContourTree::getSortedListOfMeshVertices(
    McDArray< mculong > & vertexIds)
{
    vertexIds = m_sortedNodeIdx;
}
