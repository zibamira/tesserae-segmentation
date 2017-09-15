#include "HxCorrectAndAnalyzeLabelfield.h"

#include <hxcore/HxObjectPool.h>
#include <hxcore/HxMessage.h>
#include <hxcore/HxMatDatabase.h>

#include <hxspatialgraph/internal/SpatialGraphSelection.h>

#include <amiramesh/HxParamBundle.h>

#include <mclib/internal/McDMatrix.h>
#include <mclib/internal/McSorter.h>
#include <mclib/internal/McComparators.h>
#include <mclib/McBitfield.h>
#include <mclib/internal/McVoxelChecker.h>
#include <mclib/McData3D.h>

#include <hxfield/HxUniformVectorField3.h>
#include <hxfield/HxRegField3.h>

#include <hxlineviewer/HxLineRaycast.h>

#include <hxsurftools/internal/HxDisplaySurface.h>
#include <hxsurftools/internal/HxGMC.h>

#include <hxsurface/HxSurfaceScalarField.h>

#include <mcsegalgo/internal/McFloodFill3D.h>

#include <hxvolumeviz2/internal/HxVoxelizedRender.h>
#include <hxvolumeviz2/internal/HxVolumeRenderingSettings.h>

HX_INIT_CLASS(HxCorrectAndAnalyzeLabelfield, HxCompModule)

HxCorrectAndAnalyzeLabelfield::HxCorrectAndAnalyzeLabelfield()
    : HxCompModule(HxSpatialGraph::getClassTypeId())
    , portLabelfield(this, "labelfield", tr("Label Field"), HxUniformLabelField3::getClassTypeId())
    , portIntensityField(this, "intensityField", tr("Intensity Field"), HxUniformScalarField3::getClassTypeId())
    , portDistanceMap(this, "distancemap", tr("Distance Map"), HxUniformScalarField3::getClassTypeId())
    , portReferenceLabelfield(this, "referenceLabelField", tr("Reference Label Field"), HxUniformLabelField3::getClassTypeId())
    , portMeld(this, "meld", tr("Meld"), 2)
    , portSplit(this, "split", tr("Split"), 3)
    , portRemove(this, "remove", tr("Remove"), 2)
    , portVisualization(this, "visualization", tr("Visualization"), 3)
    , portSurface(this, "surface", tr("Surface"))
    , portMisc(this, "misc", tr("Misc"), 5)
    , portEditBuffer(this, "editbuffer", tr("Edit Buffer"), 5)
    , portCriticalRegionSelection(this, "criticalRegionSelection", tr("Critical Region Selection"), 6)
    , portSurfaceField(this, "surfaceField", tr("Surface Field"))
    , portIntParameters(this, "intParameters", tr("Int Parameters"), 2)
    , portFloatParameters(this, "floatParameters", tr("Float Parameters"), 2)
    , portDoIt(this, "doIt", tr("Action"))
    , m_eng(0)
    , m_surface(NULL)
    , m_surfaceNotUpdated(true)
    , m_autoUpdateSurface(false)
    , m_initDone(false)
    , m_initSplitDone(false)
{
    portSurface.setLabel("Surface");
    portSurface.insertPushButton(0, "Create / Update");
    portSurface.insertCheckBox(1, "Auto update");

    portVisualization.setLabel("Show");
    portVisualization.setLabel(0, "All");
    portVisualization.setLabel(1, "Critical regions");
    portVisualization.setLabel(2, "Selected regions");

    portEditBuffer.setLabel("Critical regions");
    portEditBuffer.setLabel(0, "Clear");
    portEditBuffer.setLabel(1, "Add");
    portEditBuffer.setLabel(2, "Remove");
    portEditBuffer.setLabel(3, "Add selected");
    portEditBuffer.setLabel(4, "Remove selected");

    portRemove.setLabel(0, "All critical regions");
    portRemove.setLabel(1, "Selected regions");

    portMeld.setLabel(0, "All critical regions");
    portMeld.setLabel(1, "Selected regions");

    portMisc.setLabel(0, "Improve colors");
    portMisc.setLabel(1, "Print vertices");
    portMisc.setLabel(2, "Thresholding");
    portMisc.setLabel(3, "Get clean graph");
    portMisc.setLabel(4, "Remove regions using reference field");

    portCriticalRegionSelection.setLabel(0, "Small regions");
    portCriticalRegionSelection.setLabel(1, "Near regions");
    portCriticalRegionSelection.setLabel(2, "Inner regions");
    portCriticalRegionSelection.setLabel(3, "Number of neighbours");
    portCriticalRegionSelection.setLabel(4, "Size of connected components");
    portCriticalRegionSelection.setLabel(5, "Just connected to exterior");

    portSplit.setLabel(0, "Two regions contour tree");
    portSplit.setLabel(1, "Spectral clustering");
    portSplit.setLabel(2, "Standard contour tree");
}

HxCorrectAndAnalyzeLabelfield::~HxCorrectAndAnalyzeLabelfield()
{
    // Close matlab engine
    if (m_eng)
    {
        theMatlabEng->engClose(m_eng);
    }
}

void
HxCorrectAndAnalyzeLabelfield::update()
{
    if (portData.isNew(HxData::NEW_SOURCE))
    {
        m_graph = hxconnection_cast<HxSpatialGraph>(portData);
        m_initDone = false;
        m_initSplitDone = false;
        portDoIt.setSensitivity(0, 1);

        if (m_graph)
        {
            QStringList surfaceFieldLabels;
            const int numberOfAttributes = m_graph->numAttributes(HxSpatialGraph::VERTEX);
            for (int i = 0; i < numberOfAttributes; i++)
            {
                GraphAttribute* attribute = m_graph->attribute(HxSpatialGraph::VERTEX, i);
                surfaceFieldLabels << attribute->getName();
            }
            const int numEntries = portSurfaceField.getNum();
            for (int i = numEntries - 1; i >= 0; i--)
            {
                portSurfaceField.deleteItem(i);
            }
            portSurfaceField.insertComboBox(0, surfaceFieldLabels);
            portSurfaceField.insertPushButton(1, "Create");
        }
        else
        {
            const int numEntries = portSurfaceField.getNum();
            for (int i = numEntries - 1; i >= 0; i--)
            {
                portSurfaceField.deleteItem(i);
            }
        }
    }
    if (portLabelfield.isNew(HxData::NEW_SOURCE))
    {
        m_labelField = hxconnection_cast<HxUniformLabelField3>(portLabelfield);
        m_initDone = false;
        m_initSplitDone = false;
        portDoIt.setSensitivity(0, 1);
    }
    if (portDistanceMap.isNew(HxData::NEW_SOURCE))
    {
        m_initSplitDone = false;
    }
    if (portSurface.isNew(HxData::NEW_SOURCE))
    {
        m_autoUpdateSurface = portSurface.getValue(1);
        if (m_initDone && m_autoUpdateSurface && m_surfaceNotUpdated)
        {
            computeSurface();
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::compute()
{
    if (portDoIt.wasHit())
    {
        if (!m_initDone)
        {
            init();
        }
    }

    if (!m_initDone)
    {
        return;
    }

    if (portSurface.getValue(0))
    {
        computeSurface();
    }
    if (portVisualization.wasHit(0))
    {
        showAllLabels(!m_surfaceNotUpdated);
    }
    if (portVisualization.wasHit(1))
    {
        showCriticalLabels(!m_surfaceNotUpdated);
    }
    if (portVisualization.wasHit(2))
    {
        showSelectedLabels(!m_surfaceNotUpdated);
    }
    if (portEditBuffer.wasHit(0))
    {
        m_criticalNodes.clear();
    }
    if (portEditBuffer.wasHit(1))
    {
        addToCriticalNodes();
    }
    if (portEditBuffer.wasHit(2))
    {
        removeFromCriticalNodes();
    }
    if (portEditBuffer.wasHit(3))
    {
        addSelectedNodesToCriticalNodes();
    }
    if (portEditBuffer.wasHit(4))
    {
        removeSelectedNodesFromCriticalNodes();
    }
    if (portRemove.wasHit(0))
    {
        removeAllCriticalRegions();
    }
    if (portRemove.wasHit(1))
    {
        removeSelectedRegions();
    }
    if (portMeld.wasHit(0))
    {
        if (m_surfaceNotUpdated)
        {
            theMsg->printf("Update the surface");
            return;
        }
        meldAllCriticalRegions();
    }
    if (portMeld.wasHit(1))
    {
        meldSelectedRegionsToOneRegion();
    }
    if (portMisc.wasHit(0))
    {
        adjustLabelColors();
    }
    if (portMisc.wasHit(1))
    {
        printSelectedVertices();
    }
    if (portMisc.wasHit(2))
    {
        m_surfaceNotUpdated = true;
        float threshold = portFloatParameters.getValue(0);
        thresholdingOnSelectedSegments(threshold);
    }
    if (portMisc.wasHit(3))
    {
        createGraphForCleanLabels();
    }
    if (portMisc.wasHit(4))
    {
        removeRegionsAccordingToReferenceLabelField();
    }
    if (portSurfaceField.getValue(1))
    {
        if (m_surfaceNotUpdated)
        {
            theMsg->printf("Update the surface");
            return;
        }
        QString surfaceFieldString;
        portSurfaceField.getStringValue(0, surfaceFieldString);
        createSurfaceField(surfaceFieldString);
    }
    if (portSplit.wasHit(0))
    {
        if (initSplit())
        {
            m_surfaceNotUpdated = true;
            twoRegionsContourTreeSplit();
        }
    }
    if (portSplit.wasHit(1))
    {
        if (initSplit())
        {
            m_surfaceNotUpdated = true;
            spectralClusteringSplit();
        }
    }
    if (portSplit.wasHit(2))
    {
        if (initSplit())
        {
            m_surfaceNotUpdated = true;
            standardContourTreeSplit();
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::init()
{
    if ((!m_labelField) || (!m_graph))
    {
        theMsg->printf("Attach graph and label field to initialize the module.");
        return;
    }

    McPrimType dataType = m_labelField->lattice().primType();
    if (dataType != McPrimType::MC_UINT16)
    {
        theMsg->printf("Cast label field to 16 bit");
        return;
    }

    // Non-statistic attributes
    EdgeVertexAttribute* labelAttribute = m_graph->findVertexAttribute("Label");
    EdgeVertexAttribute* criticalAttribute = dynamic_cast<EdgeVertexAttribute*>(m_graph->addAttribute(
        "Critical", HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));

    if (!labelAttribute)
    {
        theMsg->printf("Attached graph has to have an attribute called Label containing the corresponding labels in the attached label field.");
        return;
    }

    // Initialize label to graph node mapping
    m_labelToNode.clear();
    m_nodeToLabel.clear();

    float minLabel, maxLabel;
    m_labelField->getRange(minLabel, maxLabel);

    const int nVertices = m_graph->getNumVertices();
    const int nLabels = (int)maxLabel;
    m_labelToNode.resize(nLabels);
    m_nodeToLabel.resize(nVertices);
    m_labelToNode.fill(-1);

    for (int i = 0; i < nVertices; i++)
    {
        const int nodeLabel = labelAttribute->getIntDataAtIdx(i);
        m_labelToNode[nodeLabel - 1] = i;
        m_nodeToLabel[i] = nodeLabel - 1;

        criticalAttribute->setIntDataAtIdx(i, 0);
    }

    m_criticalNodes.resize(0);

    calculateSizeOfRegion();
    createCenters();

    // Create volume rendering
    m_labelField->labelLattice().makeColormap();
    McString result = m_labelField->getLabel();
    // Search for the last dot, if any.
    int dot = result.rindex('.');
    if (dot >= 0)
    {
        // Remove the last present extension.
        result.removeLast(result.length() - dot - 1);
    }
    else
    {
        result += ".";
    }
    result += "colors";
    m_colormap = (HxColormap256*)theObjectPool->findObject(result);
    m_colormap->setChannel(0, 3, 0.0F);
    m_colormap->setChannel(1, 3, 0.0F);
    theObjectPool->addObject((HxObject*)m_colormap);

    m_volumeRenderingSettings = HxVolumeRenderingSettings::createInstance();
    m_volumeRenderingSettings->setLabel("VolumeRenderingSettings");
    m_volumeRenderingSettings->portData.connect(m_labelField);
    m_volumeRenderingSettings->fire();
    theObjectPool->addObject(m_volumeRenderingSettings);

    HxVoxelizedRender* voxelizedRenderer = HxVoxelizedRender::createInstance();
    theObjectPool->addObject(voxelizedRenderer);
    voxelizedRenderer->setLabel("VoxelizedRenderer");
    voxelizedRenderer->portVolumeRenderingSettings.connect(m_volumeRenderingSettings);
    voxelizedRenderer->fire();
    voxelizedRenderer->portColormap.connect((HxObject*)m_colormap);
    voxelizedRenderer->fire();

    // Create lineraycast
    m_lineRayCast = HxLineRaycast::createInstance();
    m_lineRayCast->setLabel("LineRaycast");
    m_lineRayCast->portData.connect(m_graph);
    m_lineRayCast->fire();
    theObjectPool->addObject(m_lineRayCast);

    m_initDone = true;
    portDoIt.setSensitivity(0, 0);
}

void
HxCorrectAndAnalyzeLabelfield::addToCriticalNodes()
{
    McDArray<mculong> nodes;
    getSelectedCriticalNodes(nodes);
    m_criticalNodes.appendArray(nodes);
}

void
HxCorrectAndAnalyzeLabelfield::addSelectedNodesToCriticalNodes()
{
    McDArray<int> selectedVerticesInt = m_lineRayCast->getSelectedVertices();
    McDArray<mculong> selectedVertices(selectedVerticesInt.size());
    for (int i = 0; i < selectedVerticesInt.size(); i++)
    {
        selectedVertices[i] = (mculong)selectedVerticesInt[i];
    }
    m_criticalNodes.appendArray(selectedVertices);
}

void
HxCorrectAndAnalyzeLabelfield::removeSelectedNodesFromCriticalNodes()
{
    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    for (mclong i = 0; i < selectedVertices.size(); i++)
    {
        removeEntryFromCriticalNodes(selectedVertices[i]);
    }
}

void
HxCorrectAndAnalyzeLabelfield::removeFromCriticalNodes()
{
    McDArray<mculong> nodes;
    getSelectedCriticalNodes(nodes);
    for (mclong i = 0; i < nodes.size(); i++)
    {
        removeEntryFromCriticalNodes(nodes[i]);
    }
}

void
HxCorrectAndAnalyzeLabelfield::removeEntryFromCriticalNodes(mculong removeId)
{
    for (mclong i = 0; i < m_criticalNodes.size(); i++)
    {
        if (m_criticalNodes[i] == removeId)
        {
            m_criticalNodes.remove(i);
            return;
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getSelectedCriticalNodes(McDArray<mculong>& nodes)
{
    if (portCriticalRegionSelection.getValue(0) == SELECTION_SMALL)
    {
        getCriticalSmallLabels(nodes);
    }
    if (portCriticalRegionSelection.getValue(0) == SELECTION_NEAR)
    {
        getCriticalNearNodes(nodes);
    }
    if (portCriticalRegionSelection.getValue(0) == SELECTION_INSIDE)
    {
        getCriticalInsideLabels(nodes);
    }
    if (portCriticalRegionSelection.getValue(0) == SELECTION_NEIGHBOURS)
    {
        getCriticalNumberOfNeighbours(nodes);
    }
    if (portCriticalRegionSelection.getValue(0) == SELECTION_CONNECTED_COMPONENTS)
    {
        getCriticalConnectedComponentNodes(nodes);
    }
    if (portCriticalRegionSelection.getValue(0) == SELECTION_EXTERIOR)
    {
        getCriticalJustExteriorLabels(nodes);
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalSmallLabels(McDArray<mculong>& nodes)
{
    const mculong minSize = portIntParameters.getValue(0);

    for (mclong i = 1; i < m_sizeOfRegion.size(); i++)
    {
        if ((m_sizeOfRegion[i] < minSize) && (m_labelToNode[i - 1] != -1))
        {
            nodes.append(m_labelToNode[i - 1]);
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalInsideLabels(McDArray<mculong>& nodes)
{
    if (m_surfaceNotUpdated)
    {
        theMsg->printf("Update the surface");
        return;
    }

    const int nVertices = m_graph->getNumVertices();

    McDArray<bool> labelConnectedToExterior; // without ext
    labelConnectedToExterior.resize(nVertices);
    labelConnectedToExterior.fill(0);

    const McDArray<Surface::Patch*> patches = m_surface->patches;
    for (int i = 0; i < patches.size(); i++)
    {
        const int nonExteriorMaterial = getNonExteriorMaterial(patches[i]);
        if (nonExteriorMaterial != -1)
        {
            // nonExteriorMaterial has connection to exterior
            labelConnectedToExterior[nonExteriorMaterial - 1] = 1;
        }
    }
    for (int i = 0; i < labelConnectedToExterior.size(); i++)
    {
        if (!labelConnectedToExterior[i])
        {
            nodes.append(m_labelToNode[i]);
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalJustExteriorLabels(McDArray<mculong>& nodes)
{
    if (m_surfaceNotUpdated)
    {
        theMsg->printf("Update the surface");
        return;
    }

    const int nVertices = m_graph->getNumVertices();

    McDArray<bool> labelJustConnectedToExterior; // without ext
    labelJustConnectedToExterior.resize(nVertices);
    labelJustConnectedToExterior.fill(1);

    const McDArray<Surface::Patch*> patches = m_surface->patches;
    for (int i = 0; i < patches.size(); i++)
    {
        const int nonExteriorMaterial = getNonExteriorMaterial(patches[i]);
        if (nonExteriorMaterial == -1)
        {
            // Unset both materials
            labelJustConnectedToExterior[patches[i]->innerRegion - 1] = 0;
            labelJustConnectedToExterior[patches[i]->outerRegion - 1] = 0;
        }
    }
    for (int i = 0; i < labelJustConnectedToExterior.size(); i++)
    {
        if (labelJustConnectedToExterior[i])
        {
            nodes.append(m_labelToNode[i]);
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalConnectedComponentNodes(McDArray<mculong>& nodes)
{
    const int threshold = portIntParameters.getValue(0);

    const int nVertices = m_graph->getNumVertices();

    for (int i = 0; i < nVertices; i++)
    {
        SpatialGraphSelection connectedComponent = m_graph->getConnectedComponent(i);

        const int sizeOfConnectedComponent = connectedComponent.getNumSelectedVertices();
        if (sizeOfConnectedComponent < threshold)
        {
            nodes.append(i);
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalNearNodes(McDArray<mculong>& nodes)
{
    const float minDistance = portFloatParameters.getValue(0);

    const int nEdges = m_graph->getNumEdges();
    const int nVertices = m_graph->getNumVertices();

    McDArray<bool> nodeAdded; // without ext
    nodeAdded.resize(nVertices);
    nodeAdded.fill(0);

    for (int i = 0; i < nEdges; i++)
    {
        const int sourceNodeId = m_graph->getEdgeSource(i);
        const int targetNodeId = m_graph->getEdgeTarget(i);
        const McVec3f sourceNode = m_graph->getVertexCoords(sourceNodeId);
        const McVec3f targetNode = m_graph->getVertexCoords(targetNodeId);
        const McVec3f distanceVector = targetNode - sourceNode;

        if (distanceVector.length() < minDistance)
        {
            if (!nodeAdded[sourceNodeId])
            {
                nodes.append(sourceNodeId);
                nodeAdded[sourceNodeId] = 1;
            }
            if (!nodeAdded[targetNodeId])
            {
                nodes.append(targetNodeId);
                nodeAdded[targetNodeId] = 1;
            }
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::getCriticalNumberOfNeighbours(McDArray<mculong>& nodes)
{
    const int maxNeighbours = portIntParameters.getValue(1);
    const int minNeighbours = portIntParameters.getValue(0);

    const int nVertices = m_graph->getNumVertices();

    for (int i = 0; i < nVertices; i++)
    {
        const IncidenceList incidenceList = m_graph->getIncidentEdges(i);
        if ((incidenceList.size() > maxNeighbours) ||
            (incidenceList.size() < minNeighbours))
        {
            nodes.append(i);
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::showCriticalLabels(const bool useSurface)
{
    McDArray<mclong> criticalLabels;
    for (mclong i = 0; i < m_criticalNodes.size(); i++)
    {
        criticalLabels.insertSorted(m_nodeToLabel[m_criticalNodes[i]] + 1, &mcStandardCompare);
    }

    showLabels(criticalLabels, useSurface);
}

void
HxCorrectAndAnalyzeLabelfield::showSelectedLabels(const bool useSurface)
{
    McDArray<mclong> selectedLabels;
    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        selectedLabels.insertSorted(m_nodeToLabel[selectedVertices[i]] + 1, &mcStandardCompare);
    }

    showLabels(selectedLabels, useSurface);
}

void
HxCorrectAndAnalyzeLabelfield::showLabels(McDArray<mclong> sortedLabels, const bool useSurface)
{
    // Volume rendering of label field
    for (int i = 0; i < m_colormap->getSize(); i++)
    {
        m_colormap->setChannel(i, 3, 0.0F);
    }
    for (int i = 0; i < sortedLabels.sizeInt(); i++)
    {
        m_colormap->setChannel(sortedLabels[i] + 1, 3, 1.0F);
    }
    m_colormap->touch();
    m_volumeRenderingSettings->fire();

    // Surface visualization
    if (useSurface)
    {
        m_surfaceDisplay->triSurface().selectAllTriangles(0);

        const mclong nTriangles = m_surface->triangles().size();
        for (mclong i = 0; i < nTriangles; i++)
        {
            Surface::Patch* patch = m_surface->patches[m_surface->triangles()[i].patch];
            if ((sortedLabels.findSorted(patch->innerRegion, &mcStandardCompare) != -1) ||
                (sortedLabels.findSorted(patch->outerRegion, &mcStandardCompare) != -1))
            {
                // patch->innerRegion or patch->outerRegion aus kritschem label
                m_surfaceDisplay->triSurface().selectTriangle(i, 1);
            }
        }

        m_surfaceDisplay->portBuffer.setValue(3);
        m_surfaceDisplay->fire();
        m_surfaceDisplay->portBuffer.setValue(3);
        m_surfaceDisplay->fire();
    }
}

void
HxCorrectAndAnalyzeLabelfield::showAllLabels(const bool useSurface)
{
    // Volume rendering of label field
    for (int i = 2; i < m_colormap->getSize(); i++)
    {
        m_colormap->setChannel(i, 3, 1.0F);
    }
    m_colormap->touch();
    m_volumeRenderingSettings->fire();

    // Surface visualization
    if (useSurface)
    {
        m_surfaceDisplay->triSurface().selectAllTriangles(1);

        m_surfaceDisplay->portBuffer.setValue(3);
        m_surfaceDisplay->fire();
        m_surfaceDisplay->portBuffer.setValue(3);
        m_surfaceDisplay->fire();
    }
}

void
HxCorrectAndAnalyzeLabelfield::computeSurface()
{
    if (m_surface)
    {
        m_surfaceGen->portAction.setValue(0);
        m_surfaceGen->fire();

        m_surface->touch();
        m_surface->fire();
    }
    else
    {
        // Create surface
        m_surfaceGen = HxGMC::createInstance();
        m_surfaceGen->setLabel("SurfaceGen");
        theObjectPool->addObject(m_surfaceGen);
        m_surfaceGen->portData.connect(m_labelField);
        m_surfaceGen->fire();
        m_surfaceGen->portSmoothing.setValue(0);
        m_surfaceGen->portAction.setValue(0);
        m_surfaceGen->fire();
        m_surface = dynamic_cast<HxSurface*>(m_surfaceGen->getResult());
        theObjectPool->addObject(m_surface);

        m_surfaceDisplay = HxDisplaySurface::createInstance();
        m_surfaceDisplay->setLabel("SurfaceView");
        m_surfaceDisplay->portData.connect(m_surface);
        m_surfaceDisplay->fire();
        theObjectPool->addObject(m_surfaceDisplay);
    }
    m_surfaceNotUpdated = false;
}

void
HxCorrectAndAnalyzeLabelfield::touchInputLabelField()
{
    m_labelField->touchMinMax();
    m_labelField->touch();
    m_labelField->fire();
}

void
HxCorrectAndAnalyzeLabelfield::recomputeColormap()
{
    touchInputLabelField();

    m_labelField->fire();
    m_labelField->labelLattice().createMissingMaterials(m_labelField->labelLattice().materials(), HxLabelLattice3::RANDOM_COLORS);
    m_labelField->fire();
    const int nMaterials = m_labelField->labelLattice().materials()->getNumberOfBundles();

    // Update colormap
    m_colormap->resize(nMaterials + 2);
    m_colormap->setLabelField(false);
    m_colormap->setInterpolate(true);
    m_colormap->setMinMax(0, nMaterials - 1);
    m_colormap->setLabelField(true);
    m_colormap->setInterpolate(false);

    SbColor c;
    for (int i = 1; i <= nMaterials; i++)
    {
        HxParamBundle* m = m_labelField->labelLattice().materials()->getBundle(i - 1);
        if (!m->findColor(&c[0]))
            c = theDatabase->getColor(m->getName().toLatin1().constData());
        m_colormap->setRGBA(i, c[0], c[1], c[2], 1);
    }
    m_colormap->setChannel(0, 3, 0.0F);
    m_colormap->setChannel(1, 3, 0.0F);

    touchInputLabelField();
}

void
HxCorrectAndAnalyzeLabelfield::calculateSizeOfRegion()
{
    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();
    const McDim3l& dims = m_labelField->lattice().getDims();
    const mculong size = dims[0] * dims[1] * dims[2];

    m_sizeOfRegion.resize(m_labelToNode.sizeInt() + 1);
    m_sizeOfRegion.fill(0);

    for (mculong i = 0; i < size; i++)
    {
        mcuint16 label = labelData[i];
        m_sizeOfRegion[label] += 1;
    }
}

void
HxCorrectAndAnalyzeLabelfield::createCenters()
{
    m_centers.clear();
    m_centers.resize(m_sizeOfRegion.size());
    m_centers.fill(McVec3f(0, 0, 0));

    const McDim3l& dims = m_labelField->lattice().getDims();

    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float label;
                m_labelField->lattice().eval(i, j, k, &label);
                int labelInt = (int)label;
                if (labelInt != 0)
                {
                    McVec3f pos = m_labelField->coords()->pos(McVec3i(i, j, k));
                    m_centers[labelInt] += pos;
                }
            }
        }
    }
    for (mclong i = 1; i < m_centers.size(); i++)
    {
        if (m_sizeOfRegion[i] > 0)
        {
            m_centers[i] = m_centers[i] / m_sizeOfRegion[i];
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::createCentersAfterSplit(const int numberOfNewLabels, const int firstNewLabel, const int labelToSplit)
{
    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    m_centers[labelToSplit] = McVec3f(0, 0, 0);
    for (mclong i = 0; i < numberOfNewLabels; i++)
    {
        m_centers.append(McVec3f(0, 0, 0));
    }

    for (mclong i = 0; i < m_nodeIdsOfSplitLabel.size(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        const int label = labelData[meshVertexId];

        McVec3i gridPos(0, 0, 0);
        m_mesh->getGridPosFromIdx(meshVertexId, gridPos);
        McVec3f pos = m_labelField->coords()->pos(gridPos);
        m_centers[label] += pos;
    }

    if (m_sizeOfRegion[labelToSplit] > 0)
    {
        m_centers[labelToSplit] = m_centers[labelToSplit] / m_sizeOfRegion[labelToSplit];
    }
    for (mclong i = firstNewLabel; i < firstNewLabel + numberOfNewLabels; i++)
    {
        if (m_sizeOfRegion[i] > 0)
        {
            m_centers[i] = m_centers[i] / m_sizeOfRegion[i];
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::removeAllCriticalRegions()
{
    while (m_criticalNodes.size() > 0)
    {
        mclong sourceRegion = m_nodeToLabel[m_criticalNodes[0]] + 1;
        mclong targetRegion = 0;
        meldRegions(sourceRegion, targetRegion);
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
    recomputeColormap();
}

void
HxCorrectAndAnalyzeLabelfield::removeSelectedRegions()
{
    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    if (selectedVertices.sizeInt() == 0)
    {
        theMsg->printf("Select at least one graph vertex");
        return;
    }

    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        mclong sourceRegion = m_nodeToLabel[selectedVertices[i]] + 1;
        mclong targetRegion = 0;
        meldRegions(sourceRegion, targetRegion);

        for (int j = i + 1; j < selectedVertices.sizeInt(); j++)
        {
            if (selectedVertices[j] > selectedVertices[i])
            {
                selectedVertices[j] = selectedVertices[j] - 1;
            }
        }
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    recomputeColormap();
}

void
HxCorrectAndAnalyzeLabelfield::meldAllCriticalRegions()
{
    while (m_criticalNodes.size() > 0)
    {
        mclong sourceRegion = m_nodeToLabel[m_criticalNodes[0]] + 1;
        mclong targetRegion = getLargestNeighborRegion(sourceRegion);
        meldRegions(sourceRegion, targetRegion);
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
    recomputeColormap();
}

void
HxCorrectAndAnalyzeLabelfield::meldSelectedRegionsToOneRegion()
{
    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    if (selectedVertices.sizeInt() < 2)
    {
        theMsg->printf("Select at least two graph vertices");
        return;
    }

    for (int i = 0; i < selectedVertices.sizeInt() - 1; i++)
    {
        mclong sourceRegion = m_nodeToLabel[selectedVertices[i]] + 1;
        mclong targetRegion = m_nodeToLabel[selectedVertices[i + 1]] + 1;
        meldRegions(sourceRegion, targetRegion);

        for (int j = i + 1; j < selectedVertices.sizeInt(); j++)
        {
            if (selectedVertices[j] > selectedVertices[i])
            {
                selectedVertices[j] = selectedVertices[j] - 1;
            }
        }
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    recomputeColormap();
}

void
HxCorrectAndAnalyzeLabelfield::meldRegions(mclong sourceRegion, mclong targetRegion)
{
    const mclong sourceNode = m_labelToNode[sourceRegion - 1];
    mclong targetNode = -1;
    if (targetRegion != 0)
    {
        targetNode = m_labelToNode[targetRegion - 1];
    }

    // Update labelfield
    bool floodSuccessfull = floodRegionAndSetValue(sourceRegion, targetRegion);
    if (!floodSuccessfull)
    {
        // Area is not connected, start iteration over whole dataset
        mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();
        const McDim3l& dims = m_labelField->lattice().getDims();
        const mculong size = dims[0] * dims[1] * dims[2];
        for (mculong i = 0; i < size; i++)
        {
            if (labelData[i] == sourceRegion)
            {
                labelData[i] = (mcuint16)targetRegion;
            }
        }
    }
    m_labelField->touch();
    m_surfaceNotUpdated = true;

    // Update surface graph
    // Add graph edges from sourceNode to targetNode
    if (targetRegion != 0)
    {
        const IncidenceList incidenceListSourceNode = m_graph->getIncidentEdges(sourceNode);
        for (int i = 0; i < incidenceListSourceNode.size(); i++)
        {
            const int edge = incidenceListSourceNode[i];
            const int edgeSourceNode = m_graph->getEdgeSource(edge);
            const int edgeTargetNode = m_graph->getEdgeTarget(edge);
            if (((edgeSourceNode != sourceNode) || (edgeTargetNode != targetNode)) && ((edgeSourceNode != targetNode) || (edgeTargetNode != sourceNode)))
            {
                if (edgeSourceNode != sourceNode)
                {
                    if (!areVerticesConnected(edgeSourceNode, targetNode))
                    {
                        m_graph->setEdgeTarget(edge, targetNode);
                    }
                }
                else
                {
                    if (!areVerticesConnected(edgeTargetNode, targetNode))
                    {
                        m_graph->setEdgeSource(edge, targetNode);
                    }
                }
            }
        }
    }

    // Remove remaining edges from sourceNode
    const IncidenceList incidenceListSourceNodeAfterEdgeRearranging = m_graph->getIncidentEdges(sourceNode);
    McDArray<int> deleteEdges;
    for (int i = 0; i < incidenceListSourceNodeAfterEdgeRearranging.size(); i++)
    {
        deleteEdges.insertSorted(incidenceListSourceNodeAfterEdgeRearranging[i], &mcStandardCompare);
    }
    m_graph->removeMultipleEdges(deleteEdges);

    // Update m_sizeOfRegion and m_centers
    // Needed for center calculations
    McVec3f newCenter;
    if (targetRegion != 0)
    {
        newCenter = m_centers[sourceRegion] * m_sizeOfRegion[sourceRegion] + m_centers[targetRegion] * m_sizeOfRegion[targetRegion];
        m_sizeOfRegion[targetRegion] += m_sizeOfRegion[sourceRegion];
        newCenter = newCenter / m_sizeOfRegion[targetRegion];
        m_centers[targetRegion] = newCenter;
    }
    m_sizeOfRegion[sourceRegion] = 0;
    m_centers[sourceRegion] = McVec3f(0.0, 0.0, 0.0);
    for (mclong i = m_sizeOfRegion.size() - 1; i >= 0; i--)
    {
        if (m_sizeOfRegion[i] != 0)
        {
            break;
        }
        m_sizeOfRegion.removeLast();
        m_centers.removeLast();
    }

    // Replace targetNode
    if (targetRegion != 0)
    {
        m_graph->setVertexCoords(targetNode, m_centers[targetRegion]);
    }

    // Update internal point positions of edges of targetNode
    if (targetRegion != 0)
    {
        const IncidenceList incidenceListTargetNode = m_graph->getIncidentEdges(targetNode);
        for (int i = 0; i < incidenceListTargetNode.size(); i++)
        {
            const int edge = incidenceListTargetNode[i];
            if (m_graph->getEdgeSource(edge) == targetNode)
            {
                m_graph->setEdgeSource(edge, targetNode);
            }
            else
            {
                m_graph->setEdgeTarget(edge, targetNode);
            }
        }
    }

    // Remove sourceNode
    m_graph->removeVertex(sourceNode);

    // Update lineRayCast
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    // Update data structures
    m_labelToNode[sourceRegion - 1] = -1;
    for (mclong i = 0; i < m_labelToNode.size(); i++)
    {
        if (m_labelToNode[i] > sourceNode)
        {
            m_labelToNode[i] = m_labelToNode[i] - 1;
        }
    }
    for (mclong i = m_labelToNode.size() - 1; i >= 0; i--)
    {
        if (m_labelToNode[i] != -1)
        {
            break;
        }
        m_labelToNode.removeLast();
    }

    m_nodeToLabel.remove(sourceNode);

    for (mclong i = 0; i < m_criticalNodes.size(); i++)
    {
        if (m_criticalNodes[i] > (mculong)sourceNode)
        {
            m_criticalNodes[i] = m_criticalNodes[i] - 1;
        }
        else if (m_criticalNodes[i] == (mculong)sourceNode)
        {
            m_criticalNodes.remove(i);
            i--;
        }
    }
}

mclong
HxCorrectAndAnalyzeLabelfield::getLargestNeighborRegion(const int region)
{
    mclong largestNeighborRegion = 0;
    mculong maxSize = 0;

    mclong node = m_labelToNode[region - 1];

    // Check all connected vertices and return largest region
    const IncidenceList incidenceList = m_graph->getIncidentEdges(node);
    for (int i = 0; i < incidenceList.size(); i++)
    {
        int edge = incidenceList[i];
        int edgeSourceNode = m_graph->getEdgeSource(edge);
        int edgeTargetNode = m_graph->getEdgeTarget(edge);
        if (edgeSourceNode != node)
        {
            if (m_sizeOfRegion[m_nodeToLabel[edgeSourceNode] + 1] > maxSize)
            {
                maxSize = m_sizeOfRegion[m_nodeToLabel[edgeSourceNode] + 1];
                largestNeighborRegion = m_nodeToLabel[edgeSourceNode] + 1;
            }
        }
        else
        {
            if (m_sizeOfRegion[m_nodeToLabel[edgeTargetNode] + 1] > maxSize)
            {
                maxSize = m_sizeOfRegion[m_nodeToLabel[edgeTargetNode] + 1];
                largestNeighborRegion = m_nodeToLabel[edgeTargetNode] + 1;
            }
        }
    }

    return largestNeighborRegion;
}

void
HxCorrectAndAnalyzeLabelfield::printSelectedVertices()
{
    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        theMsg->printf("%i", selectedVertices[i]);
    }
}

void
HxCorrectAndAnalyzeLabelfield::info()
{
    infoTag("Number of critical regions: ", QString::number(m_criticalNodes.sizeInt()));
}

int
HxCorrectAndAnalyzeLabelfield::getNonExteriorMaterial(HxSurface::Patch* patch)
{
    int innerMaterial = patch->innerRegion;
    int outerMaterial = patch->outerRegion;
    if ((innerMaterial == 0) || (outerMaterial == 0))
    {
        if (innerMaterial == 0)
        {
            return outerMaterial;
        }
        else
        {
            return innerMaterial;
        }
    }
    else
    {
        return -1;
    }
}

bool
HxCorrectAndAnalyzeLabelfield::areVerticesConnected(int v1, int v2)
{
    const IncidenceList incidenceList = m_graph->getIncidentEdges(v1);
    for (int i = 0; i < incidenceList.size(); i++)
    {
        int edge = incidenceList[i];
        int edgeSourceNode = m_graph->getEdgeSource(edge);
        int edgeTargetNode = m_graph->getEdgeTarget(edge);
        if (((edgeSourceNode == v1) && (edgeTargetNode == v2)) || ((edgeSourceNode == v2) && (edgeTargetNode == v1)))
        {
            return true;
        }
    }
    return false;
}

int
HxCorrectAndAnalyzeLabelfield::getGraphEdgeNeighbor(const int edge, const int vertex)
{
    const int edgeSource = m_graph->getEdgeSource(edge);
    const int edgeTarget = m_graph->getEdgeTarget(edge);
    if (edgeSource == vertex)
    {
        return edgeTarget;
    }
    else
    {
        return edgeSource;
    }
}

int
HxCorrectAndAnalyzeLabelfield::getConnectingEdge(const int vertex1, const int vertex2)
{
    const IncidenceList incidentEdges = m_graph->getIncidentEdges(vertex1);
    for (int i = 0; i < incidentEdges.size(); i++)
    {
        if ((m_graph->getEdgeSource(incidentEdges[i]) == vertex2) ||
            (m_graph->getEdgeTarget(incidentEdges[i]) == vertex2))
        {
            return incidentEdges[i];
        }
    }
    return -1;
}

bool
HxCorrectAndAnalyzeLabelfield::floodRegionAndSetValue(const int sourceRegion, const int targetRegion)
{
    const McDim3l& dims = m_labelField->lattice().getDims();
    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    McData3D<mcuint16> labelVolume(dims[0], dims[1], dims[2], labelData);
    McEqualityChecker<mcuint16> checker;
    checker.init(sourceRegion);
    McBitfield output;
    output.resize((mculong)dims[0] * dims[1] * dims[2]);
    output.unsetAll();
    int box[6];

    McVec3f centerGlobalCoords = m_centers[sourceRegion];
    HxLoc3Regular* loc = m_labelField->lattice().coords()->createLocation();
    loc->set(centerGlobalCoords);

    float labelValue;
    int ix = loc->getIx();
    int iy = loc->getIy();
    int iz = loc->getIz();
    m_labelField->lattice().eval(ix, iy, iz, &labelValue);
    if (labelValue != sourceRegion)
    {
        theMsg->printf("Attention: center is not inside the label. Calculation will last longer");
        findLatticeCoordsOfLabel(sourceRegion, ix, iy, iz);
    }

    mculong numberOfFloodedElements = mcFloodFill3D26(labelVolume, checker, ix, iy, iz, output, box);
    if (numberOfFloodedElements != m_sizeOfRegion[sourceRegion])
    {
        theMsg->printf("Attention: area not connected. No flooding possible. Calculation will be slower.");
        return false;
    }

    mculong outputIdx = 0;
    for (int kk = box[4]; kk <= box[5]; kk++)
    {
        for (int jj = box[2]; jj <= box[3]; jj++)
        {
            for (int ii = box[0]; ii <= box[1]; ii++)
            {
                outputIdx = (mculong)ii + (mculong)dims[0] * ((mculong)jj + (mculong)dims[1] * (mculong)kk);
                if (output[outputIdx])
                {
                    labelData[outputIdx] = targetRegion;
                }
            }
        }
    }

    return true;
}

void
HxCorrectAndAnalyzeLabelfield::findLatticeCoordsOfLabel(const int label, int& i, int& j, int& k)
{
    const McDim3l& dims = m_labelField->lattice().getDims();

    for (i = 0; i < dims[0]; i++)
    {
        for (j = 0; j < dims[1]; j++)
        {
            for (k = 0; k < dims[2]; k++)
            {
                float currentLabel;
                m_labelField->lattice().eval(i, j, k, &currentLabel);
                int currentLabelInt = (int)currentLabel;
                if (currentLabelInt == label)
                {
                    return;
                }
            }
        }
    }
}

void
HxCorrectAndAnalyzeLabelfield::adjustLabelColors()
{
    const int nEdges = m_graph->getNumEdges();

    const float eps = 0.5;

    HxParamBundle* materials = m_labelField->parameters.getMaterials();

    for (int i = 0; i < nEdges; i++)
    {
        const int sourceNodeId = m_graph->getEdgeSource(i);
        const int targetNodeId = m_graph->getEdgeTarget(i);
        mclong sourceRegion = m_nodeToLabel[sourceNodeId] + 1;
        mclong targetRegion = m_nodeToLabel[targetNodeId] + 1;

        // Get colors
        float colorSourceRegion[3];
        float colorTargetRegion[3];
        HxParamBundle* materialSourceRegion = materials->getBundle(sourceRegion);
        HxParamBundle* materialTargetRegion = materials->getBundle(targetRegion);
        materialSourceRegion->findColor(colorSourceRegion);
        materialTargetRegion->findColor(colorTargetRegion);

        McVec3f colorSourceRegionVector(colorSourceRegion[0], colorSourceRegion[1], colorSourceRegion[2]);
        McVec3f colorTargetRegionVector(colorTargetRegion[0], colorTargetRegion[1], colorTargetRegion[2]);
        McVec3f colorDistanceVector = colorTargetRegionVector - colorSourceRegionVector;

        // Check if colors are too similar
        if (colorDistanceVector.length() < eps)
        {
            // Change color of targetRegion
            // New color has to be different to all neighbors
            float newColor[3];
            while (true)
            {
                newColor[0] = (float)rand() / RAND_MAX;
                newColor[1] = (float)rand() / RAND_MAX;
                newColor[2] = (float)rand() / RAND_MAX;
                McVec3f newColorVector(newColor[0], newColor[1], newColor[2]);

                if (isColorDifferentToNeighbors(targetNodeId, newColorVector, eps))
                {
                    break;
                }
            }
            materialTargetRegion->setColor(newColor);
        }
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
    recomputeColormap();
}

bool
HxCorrectAndAnalyzeLabelfield::isColorDifferentToNeighbors(const int nodeId, const McVec3f color, const float eps)
{
    HxParamBundle* materials = m_labelField->parameters.getMaterials();

    const IncidenceList incidenceList = m_graph->getIncidentEdges(nodeId);
    for (int i = 0; i < incidenceList.size(); i++)
    {
        int edgeSource = m_graph->getEdgeSource(incidenceList[i]);
        int edgeTarget = m_graph->getEdgeTarget(incidenceList[i]);

        mclong neighborRegion;
        if (edgeSource == nodeId)
        {
            neighborRegion = m_nodeToLabel[edgeTarget] + 1;
        }
        else
        {
            neighborRegion = m_nodeToLabel[edgeSource] + 1;
        }

        float colorNeighbor[3];
        HxParamBundle* materialNeighbor = materials->getBundle(neighborRegion);
        materialNeighbor->findColor(colorNeighbor);

        McVec3f colorNeighborVector(colorNeighbor[0], colorNeighbor[1], colorNeighbor[2]);
        McVec3f colorDistanceVector = colorNeighborVector - color;

        if (colorDistanceVector.length() < eps)
        {
            return false;
        }
    }
    return true;
}

bool
HxCorrectAndAnalyzeLabelfield::initSplit(const int nodeToSplit)
{
    HxUniformScalarField3* inputField = hxconnection_cast<HxUniformScalarField3>(portDistanceMap);
    if (!inputField)
    {
        theMsg->printf("Please attach a distance map to perform a split");
        return false;
    }

    int labelToSplit;
    if (nodeToSplit == -1)
    {
        McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
        if (selectedVertices.sizeInt() == 0)
        {
            theMsg->printf("Select one graph vertex");
            return false;
        }
        labelToSplit = m_nodeToLabel[selectedVertices[0]] + 1;
    }
    else
    {
        labelToSplit = m_nodeToLabel[nodeToSplit] + 1;
    }

    if (!m_initSplitDone)
    {
        m_mesh = new SimplicialMesh3DForHexahedralMesh(&(inputField->lattice()));
        m_mesh->setNeighborhood(SimplicialMesh3DForHexahedralMesh::NEIGHBORHOOD_26);
        m_mesh->setThreshold(0.000001);

        const bool increasingOrder = false;
        m_mesh->getSortedListOfVertices(increasingOrder, m_sortedNodeIdx);
    }

    m_nodeIdsOfSplitLabel.clear();
    m_nodeIdToPositionInLabel.clear();
    m_nodeIdToPositionInLabel.resize(m_sortedNodeIdx.size());
    m_nodeIdToPositionInLabel.fill(-1);

    for (mclong i = 0; i < m_sortedNodeIdx.size(); ++i)
    {
        const mculong nodeIdx = m_sortedNodeIdx[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        float labelNumber;
        m_labelField->lattice().eval(meshVertexId, &labelNumber);

        // Do just consider nodes with the given label
        if (labelToSplit == labelNumber)
        {
            m_nodeIdToPositionInLabel[nodeIdx] = m_nodeIdsOfSplitLabel.size();
            m_nodeIdsOfSplitLabel.append(nodeIdx);
        }
    }

    m_initSplitDone = true;

    return true;
}

void
HxCorrectAndAnalyzeLabelfield::updateInternalStructuresAfterSplit(McDArray<int> numberOfVoxels, const int numberOfNewLabels, const int firstNewLabel, const int labelToSplit)
{
    // Update nodes of graph and connections and set label attribute
    // Needs centers -> needs m_sizeOfRegion
    const mclong selectedNode = m_labelToNode[labelToSplit - 1];

    // Update m_sizeOfRegion
    m_sizeOfRegion[labelToSplit] = numberOfVoxels[0];
    for (int i = 1; i <= numberOfNewLabels; i++)
    {
        m_sizeOfRegion.append(numberOfVoxels[i]);
    }

    // Create new vertices on right positions
    EdgeVertexAttribute* labelAttribute = m_graph->findVertexAttribute("Label");
    createCentersAfterSplit(numberOfNewLabels, firstNewLabel, labelToSplit);
    for (int i = 1; i <= numberOfNewLabels; i++)
    {
        m_labelToNode.append(m_graph->getNumVertices());
        m_graph->addVertex(m_centers[firstNewLabel + i - 1]);
        labelAttribute->setIntDataAtIdx(m_graph->getNumVertices() - 1, firstNewLabel + i - 1);
        m_nodeToLabel.append(firstNewLabel + i - 2);
    }

    // Remove edges
    McDArray<int> deleteEdges;
    const IncidenceList incidenceListSelectedNode = m_graph->getIncidentEdges(selectedNode);
    for (int i = 0; i < incidenceListSelectedNode.size(); i++)
    {
        deleteEdges.insertSorted(incidenceListSelectedNode[i], &mcStandardCompare);
    }
    m_graph->removeMultipleEdges(deleteEdges);

    // Move selected node
    m_graph->setVertexCoords(selectedNode, m_centers[labelToSplit]);

    // Update edges
    McDArray<mculong> neighbors;
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        float currentLabel;
        m_labelField->lattice().eval(meshVertexId, &currentLabel);

        neighbors.clear();
        m_mesh->get26NeighborHoodOfMeshVertexIgnoringThreshold(meshVertexId, neighbors);

        for (int j = 0; j < neighbors.sizeInt(); ++j)
        {
            const mclong neighborNodeIdx = m_mesh->getNodeIdx(neighbors[j]);
            if (neighborNodeIdx < 0)
            {
                continue;
            }
            float neighborLabel;
            m_labelField->lattice().eval(neighbors[j], &neighborLabel);
            if (neighborLabel == 0.0)
            {
                continue;
            }
            if ((neighborLabel != currentLabel) &&
                (!m_graph->hasEdge(m_labelToNode[(int)neighborLabel - 1], m_labelToNode[(int)currentLabel - 1])))
            {
                // Connect vertex of neighborLabel with vertex of label currentLabel
                m_graph->addEdge(m_labelToNode[(int)neighborLabel - 1], m_labelToNode[(int)currentLabel - 1]);
            }
        }
    }

    recomputeColormap();
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
}

int
HxCorrectAndAnalyzeLabelfield::splitAtGivenContourTreeNode(const int labelToSplit, const int mergeAfterNumberOfNodes, SetUnionDataStructure& setUnion)
{
    HxUniformScalarField3* inputField = hxconnection_cast<HxUniformScalarField3>(portDistanceMap);
    if ((!inputField) || (!m_labelField))
    {
        return -1;
    }

    setUnion.setNumElements(m_nodeIdsOfSplitLabel.size());

    McDArray<mculong> neighbors;
    McDArray<mculong> sortedNeighborIds;
    McDArray<float> values;

    McDArray<float> maxValuesOfComponent(m_nodeIdsOfSplitLabel.size());
    McDArray<mculong> numberOfComponentNodes(m_nodeIdsOfSplitLabel.size());

    bool alreadySplit = false;
    mculong doNotMergeNodeId1 = -1;
    mculong doNotMergeNodeId2 = -1;
    int numberOfVisitedMergeNodes = 0;

    const bool increasingOrder = false;

    int totalNumberOfComponents = 0;

    // Iterate over all mesh nodes with given label
    for (mclong i = 0; i < m_nodeIdsOfSplitLabel.size(); ++i)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        const mculong component = i;
        setUnion.setSetIdOfElement(i, component);

        float value;
        inputField->lattice().eval(meshVertexId, &value);
        maxValuesOfComponent[component] = value;
        numberOfComponentNodes[component] = 1;
        totalNumberOfComponents++;

        neighbors.clear();
        values.clear();
        m_mesh->getNeighborsOfMeshVertex(increasingOrder, meshVertexId, neighbors, values);

        // Remove nodes of other labels
        for (mclong j = 0; j < neighbors.size(); ++j)
        {
            float neighborLabelNumber;
            m_labelField->lattice().eval(neighbors[j], &neighborLabelNumber);

            // Do just consider neighbor-nodes with the given label
            if (labelToSplit != neighborLabelNumber)
            {
                neighbors.remove(j);
                values.remove(j);
                j--;
            }
        }

        sortNeighbors(neighbors, values, sortedNeighborIds, setUnion);

        // Iterate over all neighbors
        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);

            const mclong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
            const mclong neighborComponent = setUnion.findSetId(m_nodeIdToPositionInLabel[neighbor]);
            if (component == neighborComponent)
                continue;

            if (j > 0)
            {
                numberOfVisitedMergeNodes++;
            }

            const float maximumValue =
                (maxValuesOfComponent[component] > maxValuesOfComponent[neighborComponent] ? maxValuesOfComponent[component] : maxValuesOfComponent[neighborComponent]);
            const mculong totalNumberOfComponentNodes =
                numberOfComponentNodes[neighborComponent] + numberOfComponentNodes[component];

            // If there was a split, it is not aloud to merge the 2 splitted regions
            if (alreadySplit)
            {
                const mclong doNotMergeComponent1 = setUnion.findSetId(m_nodeIdToPositionInLabel[doNotMergeNodeId1]);
                const mclong doNotMergeComponent2 = setUnion.findSetId(m_nodeIdToPositionInLabel[doNotMergeNodeId2]);
                if (((component == doNotMergeComponent1) || (component == doNotMergeComponent2)) &&
                    ((neighborComponent == doNotMergeComponent1) || (neighborComponent == doNotMergeComponent2)))
                {
                    continue;
                }
            }
            // If there was no split, try to figure out whether a split is useful
            if ((!alreadySplit) &&
                (j > 0) && (numberOfVisitedMergeNodes == mergeAfterNumberOfNodes))
            {
                alreadySplit = true;
                doNotMergeNodeId1 = nodeIdx;
                doNotMergeNodeId2 = neighbor;
                //theMsg->printf("split: %i : %i : %i", i, numberOfComponentNodes[neighborComponent], numberOfComponentNodes[component]);
                continue;
            }

            // Merge the components
            setUnion.mergeSetsOfElements(m_nodeIdToPositionInLabel[nodeIdx], m_nodeIdToPositionInLabel[neighbor]);
            maxValuesOfComponent[neighborComponent] = maximumValue;
            maxValuesOfComponent[component] = maximumValue;
            numberOfComponentNodes[neighborComponent] = totalNumberOfComponentNodes;
            numberOfComponentNodes[component] = totalNumberOfComponentNodes;
            totalNumberOfComponents--;
        }
    }

    if (alreadySplit)
    {
        // There was a split, return split quality
        int numberOfNodesInFirstComponent = numberOfComponentNodes[setUnion.findSetId(m_nodeIdToPositionInLabel[doNotMergeNodeId1])];
        int numberOfNodesInSecondComponent = numberOfComponentNodes[setUnion.findSetId(m_nodeIdToPositionInLabel[doNotMergeNodeId2])];
        return abs(numberOfNodesInFirstComponent - numberOfNodesInSecondComponent);
    }
    else
    {
        // There was no split
        if (totalNumberOfComponents > 1)
        {
            theMsg->printf("Be careful: it seems that the segment consists of multiple not connected parts");
        }
        return -1;
    }
}

void
HxCorrectAndAnalyzeLabelfield::twoRegionsContourTreeSplit(const int nodeToSplit)
{
    int labelToSplit;
    if (nodeToSplit == -1)
    {
        McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
        if (selectedVertices.sizeInt() == 0)
        {
            theMsg->printf("Select one graph vertex");
            return;
        }
        labelToSplit = m_nodeToLabel[selectedVertices[0]] + 1;
    }
    else
    {
        labelToSplit = m_nodeToLabel[nodeToSplit] + 1;
    }

    int mergeAfterNumberOfNodes = 1;

    int currentDiff = 1;
    int bestDiff = 10000000;
    int bestMergeNode = -1;

    SetUnionDataStructure setUnion;

    // Find the best split node
    while (currentDiff >= 0)
    {
        currentDiff = splitAtGivenContourTreeNode(labelToSplit, mergeAfterNumberOfNodes, setUnion);

        if ((currentDiff < bestDiff) && (currentDiff >= 0))
        {
            bestDiff = currentDiff;
            bestMergeNode = mergeAfterNumberOfNodes;
        }

        mergeAfterNumberOfNodes++;
    }

    // Recompute best result
    currentDiff = splitAtGivenContourTreeNode(labelToSplit, bestMergeNode, setUnion);

    setTwoRegionsContourTreeSplitResult(setUnion, labelToSplit);
}

void
HxCorrectAndAnalyzeLabelfield::setTwoRegionsContourTreeSplitResult(SetUnionDataStructure& setUnion, const int labelToSplit)
{
    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    // Determine new label number (which was not in use)
    float lRange, rRange;
    m_labelField->getRange(lRange, rRange);
    const mcuint16 newLabel = (mcuint16)(((int)rRange) + 1);

    // Set new label number, use union-find information
    mclong firstComponent = -1;
    McDArray<int> numberOfVoxels;
    numberOfVoxels.resize(2);
    numberOfVoxels.fill(0);
    for (mclong i = 0; i < m_nodeIdsOfSplitLabel.size(); ++i)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        const mclong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
        if (firstComponent == -1)
        {
            firstComponent = component;
        }
        if (component != firstComponent)
        {
            numberOfVoxels[1] += 1;
            labelData[meshVertexId] = newLabel;
        }
        else
        {
            numberOfVoxels[0] += 1;
            labelData[meshVertexId] = (mcuint16)labelToSplit;
        }
    }

    const int numberOfNewLabels = 1;
    updateInternalStructuresAfterSplit(numberOfVoxels, numberOfNewLabels, newLabel, labelToSplit);

    theMsg->printf("First tessera has %i voxels, second tessera has %i voxels", numberOfVoxels[0], numberOfVoxels[1]);

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    // Temporary labelfield containing only split labels for fast surface generation
    /*HxUniformLabelField3 *temporaryLabelField = dynamic_cast<HxUniformLabelField3*> ( theObjectPool->findObject("SplitResult") );
    if ( !temporaryLabelField )
    {
        const int *dims = m_labelField->lattice.dims();
        const float * bbox = m_labelField->bbox();
        temporaryLabelField = new HxUniformLabelField3(dims, McPrimType::mc_uint16);
        temporaryLabelField->setLabel("SplitResult");
        temporaryLabelField->lattice.setBoundingBox(bbox);
        theObjectPool->addObject(temporaryLabelField);
    }
    mcuint16 *temporaryLabelData = (mcuint16 *)temporaryLabelField->lattice.dataPtr();
    for ( int i = 0; i<m_nodeIdsOfSplitLabel.sizeInt(); i++ )
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);
        const mculong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
        if ( component == firstComponent )
        {
            temporaryLabelData[meshVertexId] = 1;
        }
        else
        {
            temporaryLabelData[meshVertexId] = 2;
        }
    }*/
}

void
HxCorrectAndAnalyzeLabelfield::sortNeighbors(
    const McDArray<mculong>& neighbors,
    const McDArray<float>& values,
    McDArray<mculong>& sortedNeighborIds,
    SetUnionDataStructure& setUnion)
{
    sortNeighborsAccordingToMajorityVote(neighbors, sortedNeighborIds, setUnion);
}

void
HxCorrectAndAnalyzeLabelfield::sortNeighborsAccordingToValues(
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
HxCorrectAndAnalyzeLabelfield::sortNeighborsAccordingToMajorityVote(
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
        const mcuint32 label = (mcuint32)(setUnion.findSetId(m_nodeIdToPositionInLabel[m_mesh->getNodeIdx(neighbors[i])]) + 1);

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
HxCorrectAndAnalyzeLabelfield::spectralClusteringSplit(const int nodeToSplit)
{
    int labelToSplit;
    if (nodeToSplit == -1)
    {
        McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
        if (selectedVertices.sizeInt() == 0)
        {
            theMsg->printf("Select one graph vertex");
            return;
        }
        labelToSplit = m_nodeToLabel[selectedVertices[0]] + 1;
    }
    else
    {
        labelToSplit = m_nodeToLabel[nodeToSplit] + 1;
    }
    const int numberOfLabelNodes = m_sizeOfRegion[labelToSplit];

    const int wishedNumberOfSegments = portIntParameters.getValue(0);
    if (wishedNumberOfSegments < 2)
    {
        theMsg->printf("Write the number of labels you want the selected label to split into (at least 2) in the first <Int Parameters> port.");
        return;
    }

    if (!m_eng)
    {
        m_eng = theMatlabEng->getMatlabEngine(NULL);
#ifdef HX_OS_WIN
        if (!m_eng)
        {
            theMsg->printf("There was an error opening the Matlab engine.");
            theMsg->printf("Please check your Matlab installation. On Windows, you may need to register the Matlab engine by typing");
            theMsg->printf("    matlab /regserver");
            theMsg->printf("on the command line.\n");
            return;
        }

        theMatlabEng->engSetVisible(m_eng, 0);
#else
        if (!m_eng)
        {
            theMsg->printf("There was an error opening the Matlab engine.");
            theMsg->printf("Please check your Matlab installation. Consider the important notes in the module help.");
            return;
        }
#endif
    }
    theMatlabEng->engEvalString(m_eng, "clear all;");
    char str[256];

    mxArray* indicesMxArray1 = theMatlabEng->mxCreateNumericMatrix(27 * numberOfLabelNodes, 1, mxDOUBLE_CLASS, mxREAL);
    mxArray* indicesMxArray2 = theMatlabEng->mxCreateNumericMatrix(27 * numberOfLabelNodes, 1, mxDOUBLE_CLASS, mxREAL);
    mxArray* valuesMxArray = theMatlabEng->mxCreateNumericMatrix(27 * numberOfLabelNodes, 1, mxDOUBLE_CLASS, mxREAL);

    double* indices1 = (double*)theMatlabEng->mxGetData(indicesMxArray1);
    double* indices2 = (double*)theMatlabEng->mxGetData(indicesMxArray2);
    double* values = theMatlabEng->mxGetPr(valuesMxArray);

    int numberOfElements = 0;
    McDArray<mculong> neighbors;

    // Create laplace matrix A
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        neighbors.clear();
        m_mesh->get26NeighborHoodOfMeshVertexIgnoringThreshold(meshVertexId, neighbors);

        int numberRelevantNeighbors = 0;

        for (int j = 0; j < neighbors.sizeInt(); ++j)
        {
            const mclong neighborNodeIdx = m_mesh->getNodeIdx(neighbors[j]);
            if (neighborNodeIdx < 0)
            {
                continue;
            }
            const mclong neighborPosInLabel = m_nodeIdToPositionInLabel[neighborNodeIdx];
            if (neighborPosInLabel != -1)
            {
                numberRelevantNeighbors++;

                indices1[numberOfElements] = i + 1;
                indices2[numberOfElements] = (int)neighborPosInLabel + 1;
                values[numberOfElements] = -1;
                numberOfElements++;
            }
        }
        indices1[numberOfElements] = i + 1;
        indices2[numberOfElements] = i + 1;
        values[numberOfElements] = numberRelevantNeighbors;
        numberOfElements++;
    }
    for (int i = numberOfElements; i < 27 * numberOfLabelNodes; i++)
    {
        indices1[i] = 1;
        indices2[i] = 1;
        values[i] = 0;
    }

    // Do computations in matlab
    /*QString matlabCode = QString("A = double(A);"
                                 "A = sparse(A);"
                                 "[V,D] = eigs(A,%1,0);"
                                 "resultClustering = kmeans(V,%1);"
                                 "clear V;"
                                 "clear D;"
                                 "clear A;").arg(wishedNumberOfSegments);*/

    // Results of eigs are the smallest eigenvalues
    // but they are not sorted (not needed for kmeans)
    theMatlabEng->engPutVariable(m_eng, "v", valuesMxArray);
    theMatlabEng->engPutVariable(m_eng, "i", indicesMxArray1);
    theMatlabEng->engPutVariable(m_eng, "j", indicesMxArray2);
    sprintf(str, "A = sparse(i,j,v,%i,%i);", numberOfLabelNodes, numberOfLabelNodes);
    theMatlabEng->engEvalString(m_eng, McString(str).getString());
    sprintf(str, "[V,D] = eigs(A,%i,0);", wishedNumberOfSegments);
    theMatlabEng->engEvalString(m_eng, McString(str).getString());
    sprintf(str, "resultClustering = kmeans(V,%i);", wishedNumberOfSegments);
    theMatlabEng->engEvalString(m_eng, McString(str).getString());
    theMatlabEng->engEvalString(m_eng, "clear V;");
    theMatlabEng->engEvalString(m_eng, "clear D;");
    theMatlabEng->engEvalString(m_eng, "clear A;");
    theMatlabEng->engEvalString(m_eng, "clear v;");
    theMatlabEng->engEvalString(m_eng, "clear i;");
    theMatlabEng->engEvalString(m_eng, "clear j;");

    // Get result
    setSpectralClusteringSplitResult(labelToSplit);
    theMatlabEng->engEvalString(m_eng, "clear resultClustering;");

    // Clean up
    theMatlabEng->mxDestroyArray(indicesMxArray1);
    theMatlabEng->mxDestroyArray(indicesMxArray2);
    theMatlabEng->mxDestroyArray(valuesMxArray);
}

void
HxCorrectAndAnalyzeLabelfield::spectralClusteringSplitCalculateFullMatrix(const int nodeToSplit)
{
    int labelToSplit;
    if (nodeToSplit == -1)
    {
        McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
        if (selectedVertices.sizeInt() == 0)
        {
            theMsg->printf("Select one graph vertex");
            return;
        }
        labelToSplit = m_nodeToLabel[selectedVertices[0]] + 1;
    }
    else
    {
        labelToSplit = m_nodeToLabel[nodeToSplit] + 1;
    }
    const int numberOfLabelNodes = m_sizeOfRegion[labelToSplit];

    const int wishedNumberOfSegments = portIntParameters.getValue(0);
    if (wishedNumberOfSegments < 2)
    {
        theMsg->printf("Write the number of labels you want the selected label to split into (at least 2) in the first <Int Parameters> port.");
        return;
    }

    if (!m_eng)
    {
        m_eng = theMatlabEng->getMatlabEngine(NULL);
#ifdef HX_OS_WIN
        if (!m_eng)
        {
            theMsg->printf("There was an error opening the Matlab engine.");
            theMsg->printf("Please check your Matlab installation. On Windows, you may need to register the Matlab engine by typing");
            theMsg->printf("    matlab /regserver");
            theMsg->printf("on the command line.\n");
            return;
        }

        theMatlabEng->engSetVisible(m_eng, 0);
#else
        if (!m_eng)
        {
            theMsg->printf("There was an error opening the Matlab engine.");
            theMsg->printf("Please check your Matlab installation. Consider the important notes in the module help.");
            return;
        }
#endif
    }
    theMatlabEng->engEvalString(m_eng, "clear all;");
    char str[256];

    // Create sparse matrix later
    mxArray* laplaceMatrix = theMatlabEng->mxCreateNumericMatrix(numberOfLabelNodes, numberOfLabelNodes, mxDOUBLE_CLASS, mxREAL);
    double* laplaceMatrixPtr = theMatlabEng->mxGetPr(laplaceMatrix);

    McDArray<mculong> neighbors;

    // Create laplace matrix A
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        neighbors.clear();
        m_mesh->get26NeighborHoodOfMeshVertexIgnoringThreshold(meshVertexId, neighbors);

        int numberRelevantNeighbors = 0;

        for (int j = 0; j < neighbors.sizeInt(); ++j)
        {
            const mclong neighborNodeIdx = m_mesh->getNodeIdx(neighbors[j]);
            if (neighborNodeIdx < 0)
            {
                continue;
            }
            const mclong neighborPosInLabel = m_nodeIdToPositionInLabel[neighborNodeIdx];
            if (neighborPosInLabel != -1)
            {
                numberRelevantNeighbors++;

                // Create sparse matrix later
                laplaceMatrixPtr[(int)neighborPosInLabel * numberOfLabelNodes + i] = -1; // A[i][(int) neighborPosInLabel] = -1;
                laplaceMatrixPtr[i * numberOfLabelNodes + (int)neighborPosInLabel] = -1; //A[(int) neighborPosInLabel][i] = -1;
            }
        }
        // Create sparse matrix later
        laplaceMatrixPtr[i * numberOfLabelNodes + i] = numberRelevantNeighbors; //A[i][i] = numberRelevantNeighbors;
    }

    // Do computations in matlab: Create sparse matrix later
    theMatlabEng->engPutVariable(m_eng, "A", laplaceMatrix);
    theMatlabEng->engEvalString(m_eng, "A = sparse(A);");
    sprintf(str, "[V,D] = eigs(A,%i,0);", wishedNumberOfSegments);
    theMatlabEng->engEvalString(m_eng, McString(str).getString());
    sprintf(str, "resultClustering = kmeans(V,%i);", wishedNumberOfSegments);
    theMatlabEng->engEvalString(m_eng, McString(str).getString());
    theMatlabEng->engEvalString(m_eng, "clear V;");
    theMatlabEng->engEvalString(m_eng, "clear D;");
    theMatlabEng->engEvalString(m_eng, "clear A;");

    // Get result
    setSpectralClusteringSplitResult(labelToSplit);
    theMatlabEng->engEvalString(m_eng, "clear resultClustering;");
}

void
HxCorrectAndAnalyzeLabelfield::setSpectralClusteringSplitResult(const int labelToSplit)
{
    if (!m_eng)
    {
        return;
    }

    mxArray* resultMxArray = theMatlabEng->engGetVariable(m_eng, "resultClustering");
    double* resultMxArrayPtr = theMatlabEng->mxGetPr(resultMxArray);

    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    // Determine new label number (which was not in use)
    float lRange, rRange;
    m_labelField->getRange(lRange, rRange);
    const mcuint16 firstNewLabel = (mcuint16)(((int)rRange) + 1);

    theMatlabEng->engEvalString(m_eng, "maxValue = max(resultClustering);");
    const mxArray* maxValueMxArray = theMatlabEng->engGetVariable(m_eng, "maxValue");
    int maxValue = theMatlabEng->mxGetScalar(maxValueMxArray);
    const int numberOfNewLabels = maxValue - 1;

    McDArray<mcuint16> clusteringResultToLabel;
    clusteringResultToLabel.resize(numberOfNewLabels + 1);
    clusteringResultToLabel[0] = (mcuint16)labelToSplit;
    for (int i = 1; i <= numberOfNewLabels; i++)
    {
        clusteringResultToLabel[i] = (mcuint16)(firstNewLabel + i - 1);
    }

    McDArray<int> numberOfVoxels;
    numberOfVoxels.resize(numberOfNewLabels + 1);
    numberOfVoxels.fill(0);

    // Update labelfield
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        double currentBin = resultMxArrayPtr[i];
        int index = (int)currentBin - 1;

        // Voxel with nodeId m_nodeIdsOfSplitLabel[i] belongs to resultClusteringField[i]
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);
        labelData[meshVertexId] = clusteringResultToLabel[index];
        numberOfVoxels[index] += 1;
    }

    updateInternalStructuresAfterSplit(numberOfVoxels, numberOfNewLabels, firstNewLabel, labelToSplit);

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    theMatlabEng->engEvalString(m_eng, "clear maxValue;");
}

void
HxCorrectAndAnalyzeLabelfield::standardContourTreeSplit(const int nodeToSplit)
{
    HxUniformScalarField3* inputField = hxconnection_cast<HxUniformScalarField3>(portDistanceMap);
    if ((!inputField) || (!m_labelField))
    {
        return;
    }

    int labelToSplit;
    if (nodeToSplit == -1)
    {
        McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
        if (selectedVertices.sizeInt() == 0)
        {
            theMsg->printf("Select one graph vertex");
            return;
        }
        labelToSplit = m_nodeToLabel[selectedVertices[0]] + 1;
    }
    else
    {
        labelToSplit = m_nodeToLabel[nodeToSplit] + 1;
    }

    SetUnionDataStructure setUnion;
    const bool increasingOrder = false;

    setUnion.setNumElements(m_nodeIdsOfSplitLabel.size());

    McDArray<mculong> neighbors;
    McDArray<mculong> sortedNeighborIds;
    McDArray<float> values;

    McDArray<float> maxValuesOfComponent(m_nodeIdsOfSplitLabel.size());
    McDArray<mculong> numberOfComponentNodes(m_nodeIdsOfSplitLabel.size());

    int numberOfVisitedMergeNodes = 0;

    int totalNumberOfComponents = 0;

    const float persistenceValue = portFloatParameters.getValue(0);

    // Iterate over all mesh nodes with given label
    for (mclong i = 0; i < m_nodeIdsOfSplitLabel.size(); ++i)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);

        const mculong component = i;
        setUnion.setSetIdOfElement(i, component);

        float value;
        inputField->lattice().eval(meshVertexId, &value);
        maxValuesOfComponent[component] = value;
        numberOfComponentNodes[component] = 1;
        totalNumberOfComponents++;

        neighbors.clear();
        values.clear();
        m_mesh->getNeighborsOfMeshVertex(increasingOrder, meshVertexId, neighbors, values);

        // Remove nodes of other labels
        for (mclong j = 0; j < neighbors.size(); ++j)
        {
            float neighborLabelNumber;
            m_labelField->lattice().eval(neighbors[j], &neighborLabelNumber);

            // Do just consider neighbor-nodes with the given label
            if (labelToSplit != neighborLabelNumber)
            {
                neighbors.remove(j);
                values.remove(j);
                j--;
            }
        }

        sortNeighbors(neighbors, values, sortedNeighborIds, setUnion);

        // Iterate over all neighbors
        for (mclong j = 0; j < sortedNeighborIds.size(); ++j)
        {
            const mclong neighbor = m_mesh->getNodeIdx(neighbors[sortedNeighborIds[j]]);

            const mclong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
            const mclong neighborComponent = setUnion.findSetId(m_nodeIdToPositionInLabel[neighbor]);
            if (component == neighborComponent)
                continue;

            if (j > 0)
            {
                numberOfVisitedMergeNodes++;
            }

            const float maximumValue =
                (maxValuesOfComponent[component] > maxValuesOfComponent[neighborComponent] ? maxValuesOfComponent[component] : maxValuesOfComponent[neighborComponent]);
            const mculong totalNumberOfComponentNodes =
                numberOfComponentNodes[neighborComponent] + numberOfComponentNodes[component];

            float persistenceComponent1 = (maxValuesOfComponent[component] - value);
            float persistenceComponent2 = (maxValuesOfComponent[neighborComponent] - value);

            if ((((maxValuesOfComponent[component] - value) >= 0) && ((maxValuesOfComponent[neighborComponent] - value) >= 0)) && ((persistenceComponent1 <= persistenceValue) || (persistenceComponent2 <= persistenceValue)))
            {
                setUnion.mergeSetsOfElements(m_nodeIdToPositionInLabel[nodeIdx], m_nodeIdToPositionInLabel[neighbor]);
                maxValuesOfComponent[neighborComponent] = maximumValue;
                maxValuesOfComponent[component] = maximumValue;
                numberOfComponentNodes[neighborComponent] = totalNumberOfComponentNodes;
                numberOfComponentNodes[component] = totalNumberOfComponentNodes;
                totalNumberOfComponents--;
            }
        }
    }

    setStandardContourtreeSplitResult(labelToSplit, setUnion);

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();
}

void
HxCorrectAndAnalyzeLabelfield::setStandardContourtreeSplitResult(const int labelToSplit, SetUnionDataStructure& setUnion)
{
    McDArray<int> componentToIndex;
    componentToIndex.resize(m_nodeIdsOfSplitLabel.sizeInt());
    componentToIndex.fill(-1);
    int numberOfComponents = 0;
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
        if (componentToIndex[component] == -1)
        {
            componentToIndex[component] = numberOfComponents;
            numberOfComponents++;
        }
    }

    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    // Determine new label number (which was not in use)
    float lRange, rRange;
    m_labelField->getRange(lRange, rRange);
    const mcuint16 firstNewLabel = (mcuint16)(((int)rRange) + 1);

    const int numberOfNewLabels = numberOfComponents - 1;

    McDArray<mcuint16> clusteringResultToLabel;
    clusteringResultToLabel.resize(numberOfNewLabels + 1);
    clusteringResultToLabel[0] = (mcuint16)labelToSplit;
    for (int i = 1; i <= numberOfNewLabels; i++)
    {
        clusteringResultToLabel[i] = (mcuint16)(firstNewLabel + i - 1);
    }

    McDArray<int> numberOfVoxels;
    numberOfVoxels.resize(numberOfNewLabels + 1);
    numberOfVoxels.fill(0);

    // Update labelfield
    for (int i = 0; i < m_nodeIdsOfSplitLabel.sizeInt(); i++)
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);
        const mculong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
        const int index = componentToIndex[component];

        labelData[meshVertexId] = clusteringResultToLabel[index];
        numberOfVoxels[index] += 1;
    }

    updateInternalStructuresAfterSplit(numberOfVoxels, numberOfNewLabels, firstNewLabel, labelToSplit);

    // Clear selection of graph nodes
    m_lineRayCast->clearSelection();
    m_lineRayCast->portData.touch();
    m_lineRayCast->fire();

    // Temporary labelfield containing only split labels for fast surface generation
    /*HxUniformLabelField3 *temporaryLabelField = dynamic_cast<HxUniformLabelField3*> ( theObjectPool->findObject("SplitResult") );
    if ( !temporaryLabelField )
    {
        const int *dims = m_labelField->lattice.dims();
        const float * bbox = m_labelField->bbox();
        temporaryLabelField = new HxUniformLabelField3(dims, McPrimType::mc_uint16);
        temporaryLabelField->setLabel("SplitResult");
        temporaryLabelField->lattice.setBoundingBox(bbox);
        theObjectPool->addObject(temporaryLabelField);
    }
    mcuint16 *temporaryLabelData = (mcuint16 *)temporaryLabelField->lattice.dataPtr();
    for ( int i = 0; i<m_nodeIdsOfSplitLabel.sizeInt(); i++ )
    {
        const mculong nodeIdx = m_nodeIdsOfSplitLabel[i];
        const mculong meshVertexId = m_mesh->getMeshVertexIdx(nodeIdx);
        const mculong component = setUnion.findSetId(m_nodeIdToPositionInLabel[nodeIdx]);
        const int index = componentToIndex[component];

        temporaryLabelData[meshVertexId] = index+1;
    }*/
}

void
HxCorrectAndAnalyzeLabelfield::thresholdingOnSelectedSegments(float threshold)
{
    HxUniformScalarField3* inputField = hxconnection_cast<HxUniformScalarField3>(portIntensityField);
    if (!inputField)
    {
        return;
    }

    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    if (selectedVertices.sizeInt() == 0)
    {
        theMsg->printf("Select at least one graph vertex");
        return;
    }

    McDArray<int> selectedLabels;
    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        selectedLabels.insertSorted(m_nodeToLabel[selectedVertices[i]] + 1, &mcStandardCompare);
    }

    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();
    const McDim3l& dims = m_labelField->lattice().getDims();
    const mculong size = dims[0] * dims[1] * dims[2];

    float value;
    for (mculong i = 0; i < size; i++)
    {
        mcuint16 label = labelData[i];
        inputField->lattice().eval(i, &value);
        if ((value < threshold) && (selectedLabels.findSorted(label, &mcStandardCompare) != -1))
        {
            labelData[i] = 0;
            m_sizeOfRegion[(mculong)label]--;
        }
    }

    createCenters();
    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        m_graph->setVertexCoords(selectedVertices[i], m_centers[selectedLabels[i]]);
    }
    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
}

void
HxCorrectAndAnalyzeLabelfield::createSurfaceField(QString attributeName)
{
    HxSurfaceScalarField* surfaceField = new HxSurfaceScalarField(m_surface, HxSurfaceScalarField::ON_TRIANGLES);
    float* surfaceFieldData = surfaceField->dataPtr();

    EdgeVertexAttribute* attribute = m_graph->findVertexAttribute(attributeName.toLatin1().constData());
    if (!attribute)
    {
        theMsg->printf("Attribute not existing");
        return;
    }

    const mclong nTriangles = m_surface->triangles().size();
    for (mclong i = 0; i < nTriangles; i++)
    {
        Surface::Patch* patch = m_surface->patches[m_surface->triangles()[i].patch];
        const int nonExteriorMaterial = getNonExteriorMaterial(patch);
        if (nonExteriorMaterial != -1)
        {
            float attributeValue = 0.0F;
            if (attribute->primType() == McPrimType::MC_FLOAT)
            {
                attributeValue = attribute->getFloatDataAtIdx(m_labelToNode[nonExteriorMaterial - 1]);
            }
            else
            {
                attributeValue = (float)attribute->getIntDataAtIdx(m_labelToNode[nonExteriorMaterial - 1]);
            }
            surfaceFieldData[i] = attributeValue;
        }
        else
        {
            surfaceFieldData[i] = 0.0F;
        }
    }

    surfaceField->setLabel("Surface Field");
    setResult(0, surfaceField);
}

void
HxCorrectAndAnalyzeLabelfield::createGraphForCleanLabels()
{
    HxSpatialGraph* outputGraph = m_graph->duplicate();
    outputGraph->composeLabel(m_graph->getLabel(), "cleanLabels");
    setResult(1, outputGraph);

    McDArray<int> emptyRegions;

    const int n = m_sizeOfRegion.sizeInt();
    for (int i = 1; i < n; i++)
    {
        if (m_sizeOfRegion[i] == 0)
        {
            emptyRegions.append(i);
        }
    }

    EdgeVertexAttribute* labelAttribute = outputGraph->findVertexAttribute("Label");
    const int nVertices = outputGraph->getNumVertices();

    for (int i = 0; i < nVertices; i++)
    {
        int label = labelAttribute->getIntDataAtIdx(i);
        int offset = 0;
        for (int j = 0; j < emptyRegions.sizeInt(); j++)
        {
            if (label > emptyRegions[j])
            {
                offset++;
            }
        }
        labelAttribute->setIntDataAtIdx(i, label - offset);
    }
}

void
HxCorrectAndAnalyzeLabelfield::removeRegionsAccordingToReferenceLabelField()
{
    HxUniformLabelField3* referenceLabelField = hxconnection_cast<HxUniformLabelField3>(portReferenceLabelfield);
    const McDim3l& dims = referenceLabelField->lattice().getDims();

    McDArray<mculong> newNumberOfVoxels(m_sizeOfRegion);

    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float labelReferenceField;
                referenceLabelField->lattice().eval(i, j, k, &labelReferenceField);
                if (labelReferenceField == 0)
                {
                    float label;
                    m_labelField->lattice().eval(i, j, k, &label);
                    if (label > 0)
                    {
                        newNumberOfVoxels[(mculong)label]--;
                    }
                }
            }
        }
    }

    for (int i = 1; i < newNumberOfVoxels.sizeInt(); i++)
    {
        if ((newNumberOfVoxels[i] == 0) && (m_labelToNode[i - 1] != -1) && (m_sizeOfRegion[i] > 0))
        {
            meldRegions(i, 0);
        }
    }

    if (m_autoUpdateSurface)
    {
        computeSurface();
    }
    recomputeColormap();
}

// From here: ---------------------  Tests ------------------------------

// Not used: test
// Floods splitRegion and sets m_nodeIdToPositionInLabel and m_nodeIdsOfSplitLabel
bool
HxCorrectAndAnalyzeLabelfield::floodRegionAndInitForSplit(const int splitRegion)
{
    HxUniformScalarField3* inputField = hxconnection_cast<HxUniformScalarField3>(portDistanceMap);

    const McDim3l& dims = m_labelField->lattice().getDims();
    mcuint16* labelData = (mcuint16*)m_labelField->lattice().dataPtr();

    McData3D<mcuint16> labelVolume(dims[0], dims[1], dims[2], labelData);
    McEqualityChecker<mcuint16> checker;
    checker.init(splitRegion);
    McBitfield output;
    output.resize((mculong)dims[0] * dims[1] * dims[2]);
    output.unsetAll();
    int box[6];

    McVec3f centerGlobalCoords = m_centers[splitRegion];
    HxLoc3Regular* loc = m_labelField->lattice().coords()->createLocation();
    loc->set(centerGlobalCoords);

    float labelValue;
    int ix = loc->getIx();
    int iy = loc->getIy();
    int iz = loc->getIz();
    m_labelField->lattice().eval(ix, iy, iz, &labelValue);
    if (labelValue != splitRegion)
    {
        theMsg->printf("Attention: center is not inside the label. Calculation will last longer");
        findLatticeCoordsOfLabel(splitRegion, ix, iy, iz);
    }

    mculong numberOfFloodedElements = mcFloodFill3D26(labelVolume, checker, ix, iy, iz, output, box);
    if (numberOfFloodedElements != m_sizeOfRegion[splitRegion])
    {
        theMsg->printf("Attention: area not connected. No flooding possible. Calculation will be slower.");
        return false;
    }

    float value;
    McDArray<float> values;
    McDArray<mclong> sortedNodeIdsOfLabel;
    mclong outputIdx = 0;
    for (int kk = box[4]; kk <= box[5]; kk++)
    {
        for (int jj = box[2]; jj <= box[3]; jj++)
        {
            for (int ii = box[0]; ii <= box[1]; ii++)
            {
                outputIdx = ii + dims[0] * (jj + dims[1] * kk);
                if (output[outputIdx])
                {
                    const mclong nodeIdx = m_mesh->getNodeIdx(outputIdx);
                    inputField->lattice().eval(outputIdx, &value);
                    const int insertPosition = values.insertSorted(value, &mcStandardCompare);
                    sortedNodeIdsOfLabel.insert(insertPosition, 1, &nodeIdx);
                }
            }
        }
    }
    for (int i = sortedNodeIdsOfLabel.sizeInt() - 1; i >= 0; i--)
    {
        const mclong nodeIdx = sortedNodeIdsOfLabel[i];
        m_nodeIdToPositionInLabel[nodeIdx] = m_nodeIdsOfSplitLabel.size();
        m_nodeIdsOfSplitLabel.append((mculong)nodeIdx);
    }

    return true;
}

// Not used: test
void
HxCorrectAndAnalyzeLabelfield::splitAll()
{
    const int nVertices = m_graph->getNumVertices();

    int numberVisitedNodeAtEnd = 0;
    for (int i = 0; i < nVertices - numberVisitedNodeAtEnd; i++)
    {
        const IncidenceList incidenceListNode = m_graph->getIncidentEdges(i);
        for (int j = 0; j < incidenceListNode.size(); j++)
        {
            const int edge = incidenceListNode[j];
            const int neighborNode = getGraphEdgeNeighbor(edge, i);

            mclong sourceRegion = m_nodeToLabel[neighborNode] + 1;
            mclong targetRegion = m_nodeToLabel[i] + 1;
            meldRegions(sourceRegion, targetRegion);
            if (neighborNode < i)
            {
                i--;
                numberVisitedNodeAtEnd++;
            }
            if (initSplit(i))
            {
                twoRegionsContourTreeSplit(i);
                //spectralClusteringSplit(i);
            }
        }
    }
}

// Not used: test
void
HxCorrectAndAnalyzeLabelfield::returnLocalLabelfield()
{
    const McDim3l& dims = m_labelField->lattice().getDims();
    const McBox3f& bbox = m_labelField->getBoundingBox();

    HxUniformLabelField3* localLabelField = new HxUniformLabelField3(dims, McPrimType::MC_UINT8);
    localLabelField->lattice().setBoundingBox(bbox);

    McDArray<int> selectedVertices = m_lineRayCast->getSelectedVertices();
    McDArray<int> selectedLabels;

    for (int i = 0; i < selectedVertices.sizeInt(); i++)
    {
        selectedLabels.insertSorted(m_nodeToLabel[selectedVertices[i]] + 1, &mcStandardCompare);
    }

    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float label;
                m_labelField->lattice().eval(i, j, k, &label);
                const int labelInt = (int)label;
                const float index = (float)(selectedLabels.findSorted(labelInt, &mcStandardCompare) + 1);
                if (index != -1)
                {
                    localLabelField->lattice().set(i, j, k, &index);
                }
            }
        }
    }

    localLabelField->setLabel("Local Label Field");
    theObjectPool->addObject(localLabelField);
}

// Not used: test
// Find nodes with user-defined attribute range
void
HxCorrectAndAnalyzeLabelfield::filter(McDArray<mculong>& nodes, QString attributeName)
{
    EdgeVertexAttribute* attribute = m_graph->findVertexAttribute(attributeName.toLatin1().constData());
    if (attribute->primType() == McPrimType::MC_FLOAT)
    {
        filterFloat(nodes, attribute);
    }
    else
    {
        filterInt(nodes, attribute);
    }
}

// Not used: test
void
HxCorrectAndAnalyzeLabelfield::filterFloat(McDArray<mculong>& nodes, EdgeVertexAttribute* attribute)
{
    const int nVertices = m_graph->getNumVertices();

    const float minValue = portFloatParameters.getValue(0);
    const float maxValue = portFloatParameters.getValue(1);

    for (int i = 0; i < nVertices; i++)
    {
        const float attributeValue = attribute->getFloatDataAtIdx(i);
        if ((attributeValue >= minValue) && (attributeValue <= maxValue))
        {
            nodes.append(i);
        }
    }
}

// Not used: test
void
HxCorrectAndAnalyzeLabelfield::filterInt(McDArray<mculong>& nodes, EdgeVertexAttribute* attribute)
{
    const int nVertices = m_graph->getNumVertices();

    const int minValue = portIntParameters.getValue(0);
    const int maxValue = portIntParameters.getValue(1);

    for (int i = 0; i < nVertices; i++)
    {
        const int attributeValue = attribute->getIntDataAtIdx(i);
        if ((attributeValue >= minValue) && (attributeValue <= maxValue))
        {
            nodes.append(i);
        }
    }
}
