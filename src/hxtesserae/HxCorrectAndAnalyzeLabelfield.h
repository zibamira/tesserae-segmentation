#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortButtonList.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortGeneric.h>
#include <hxcore/HxPortMultiMenu.h>

#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformLabelField3.h>

#include <hxsurface/HxSurface.h>

#include <mclib/McMat3d.h>
#include <mclib/McHandle.h>
#include <mclib/McPrimType.h>

#include <hxcontourtree/SimplicialMesh3DForHexahedralMesh.h>
#include <hxcontourtree/SetUnionDataStructure.h>
#include <hxcontourtree/ContourTree.h>
#include <hxcontourtree/AugmentedContourTree.h>

#include <hxmatlab/internal/HxMatlabEng.h>

#include <hxspatialgraph/internal/IncidenceList.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>

class HxVolumeRenderingSettings;

#include "api.h"

/**
 * This module was developed for the postprocessing / error correction of segmentations - especially for tesserae
 * (mineralized tiles on the endoskeleton of sharks and rays) segmentations - allowing
 * operations like label merging or splitting. Additionally, there are operations to highlight
 * unusual labels (for example regarding their size, these labels are then called "critical regions").
 * Main idea of the module: use the region adjacency graph (RAG or graph) of the label field for user interaction
 * and visualization. The RAG can be created by the module HxCreateRegionAdjacencyGraph,
 * must have a vertex attribute called "Label" containing the corresponding label
 * in the label field and is attached to the data port. The label attribute is only used
 * during the initialization step of this module. The label field can contain gaps, that means
 * labels with 0 voxels.
 *
 * This module can also be used for arbitrary, non-tesserae label fields (and their RAGs).
 *
 * The module has an initialization step (by pressing apply). For interaction with this module you have
 * to use the visualization modules directly created in this initialization. Visualization
 * is done with Voxelized Renderer or optional with a surface.
 *
 * It is not supported to save a project file and load it again. You have to use this module manually.
 *
 * Internal information:
 * - graph is updated during operations
 * - following datastructures are automatically updated:
 *     labelfield, surface graph (node positions, edges), m_sizeOfRegion,
 *     m_centers, m_criticalNodes, m_nodeToLabel, m_labelToNode, label attribute
 * - other attributes are not automatically updated after split/merge operations
 * - surface is only automatically updated after split/merge operations if m_autoUpdateSurface is true
 * - for split operations new labels are attached at the end of existing labels
 *   (example: highest label = 100, 2-split creates new label 101)
 */

class HxDisplaySurface;
class HxGMC;
class HxLineRaycast;
class HxSpatialGraph;

class HXTESSERAE_API HxCorrectAndAnalyzeLabelfield : public HxCompModule
{
    HX_HEADER(HxCorrectAndAnalyzeLabelfield);

public:
    void compute();
    void update();

    virtual void info();

    /// Mandatory
    HxConnection portLabelfield;
    /**
     * Optional
     * Used for thresholding of selected labels
     * The field is usually the input CT
     */
    HxConnection portIntensityField;
    /**
     * Optional
     * Used for all splits:
     * - contour tree based: distance map = 0 for background; other distance map values used for segmentation;
     *                       use 2d plane distance map
     * - spectral clustering: distance map = 0 for background;
     *                        use foreground segmentation or 2d distance map
     */
    HxConnection portDistanceMap;
    /**
     * Optional
     * Only used in removeRegionsAccordingToReferenceLabelField()
     */
    HxConnection portReferenceLabelfield;

    HxPortButtonList portMeld;
    HxPortButtonList portSplit;
    HxPortButtonList portRemove;
    HxPortButtonList portVisualization;
    HxPortGeneric portSurface;
    HxPortButtonList portMisc;
    HxPortButtonList portEditBuffer;
    HxPortMultiMenu portCriticalRegionSelection;
    HxPortGeneric portSurfaceField;
    HxPortIntTextN portIntParameters;
    HxPortFloatTextN portFloatParameters;
    HxPortDoIt portDoIt;

private:
    void init();

    // Computation / handling of critical regions
    void addToCriticalNodes();
    void addSelectedNodesToCriticalNodes();
    void removeSelectedNodesFromCriticalNodes();
    void removeFromCriticalNodes();
    void removeEntryFromCriticalNodes(mculong removeId);
    void getSelectedCriticalNodes(McDArray<mculong>& nodes);
    void getCriticalSmallLabels(McDArray<mculong>& nodes);
    /**
      * Look for labels with no connection to exterior
      * Needs surface
      */
    void getCriticalInsideLabels(McDArray<mculong>& nodes); // Needs surface
    void getCriticalNearNodes(McDArray<mculong>& nodes);
    void getCriticalNumberOfNeighbours(McDArray<mculong>& nodes);
    void getCriticalConnectedComponentNodes(McDArray<mculong>& nodes);
    /**
      * Look for labels that have only a connection to exterior
      * Needs surface
      */
    void getCriticalJustExteriorLabels(McDArray<mculong>& nodes); // Needs surface

    // Visualization with surface or voxelized renderer
    void showCriticalLabels(const bool useSurface);
    void showSelectedLabels(const bool useSurface);
    void showAllLabels(const bool useSurface);
    void showLabels(McDArray<mclong> sortedLabels, const bool useSurface);
    void computeSurface();
    void touchInputLabelField(); // Forces update of voxelized renderer after colormap has been changed
    /**
      * Called after each split / merge / remove operation
      * Uses material colors as new colormap entries
      */
    void recomputeColormap(); // Used after changes of material colors

    // Update of internal variables
    void calculateSizeOfRegion();
    /**
      * Needs updated m_sizeOfRegions and updated labelfield
      * Updates m_centers
      */
    void createCenters(); // with ext-label
    /**
      * Faster than createCenters
      * Only iterates over m_nodeIdsOfSplitLabel
      */
    void createCentersAfterSplit(const int numberOfNewLabels, const int firstNewLabel, const int labelToSplit);

    // Meld / remove methods
    void removeAllCriticalRegions();
    void removeSelectedRegions();
    /// Meld with largest neighbor
    void meldAllCriticalRegions();
    void meldSelectedRegionsToOneRegion();
    /**
      * Merge two regions
      * Updates labelfield, surface graph, m_sizeOfRegion, m_centers, m_criticalNodes, m_nodeToLabel, m_labelToNode
      * Input regions inclusive exterior region
      * Does not update graph attributes
      */
    void meldRegions(mclong sourceRegion, mclong targetRegion);
    /**
      * Input region inclusive exterior region
      * Output region inclusive exterior region
      */
    mclong getLargestNeighborRegion(int region);

    // Aux stuff
    void printSelectedVertices();
    int getNonExteriorMaterial(HxSurface::Patch* patch);
    bool areVerticesConnected(int v1, int v2);
    int getGraphEdgeNeighbor(const int edge, const int vertex);
    int getConnectingEdge(const int vertex1, const int vertex2);
    /**
      * Start flooding from center of sourceRegion and assign voxels to targetRegion
      * Does only work if sourceRegion is connected (function returns true)
      * Otherwise function doesn't change anything and returns false
      */
    bool floodRegionAndSetValue(const int sourseRegion, const int targetRegion);
    void findLatticeCoordsOfLabel(const int label, int& i, int& j, int& k);

    // Color adjustment
    /**
      * Changes material colors
      * Then the colormap is updated and if necessary the surface is updated
      */
    void adjustLabelColors();
    bool isColorDifferentToNeighbors(const int nodeId, const McVec3f color, const float eps);

    // For all splits
    /**
      * Uses portDistanceMap to initialize mesh, everything larger 0 in given label is used
      * All voxels in given label must have distance map value larger 0: Problems can occur with prior thresholding
      * Intializes m_mesh, m_nodeIdsOfSplitLabel, m_nodeIdToPositionInLabel
      * m_initSplitDone set to true
      */
    bool initSplit(const int nodeToSplit = -1);
    /**
      * Update surface graph, m_sizeOfRegion, m_centers, m_nodeToLabel, m_labelToNode, label attribute
      * Does not update graph attributes
      * Needs correct label field
      */
    void updateInternalStructuresAfterSplit(McDArray<int> numberOfVoxels, const int numberOfNewLabels, const int firstNewLabel, const int labelToSplit);

    // Functions used for contour tree splits:
    // 1. Region is exactly split in two (on node where size of resulting regions is as similar as possible)
    // 2. Standard cotour tree split with user-given persistence
    /**
      * Split one label in exactly 2 parts according to the contour tree of the
      * gray values (just on positions with the given label)
      */
    void twoRegionsContourTreeSplit(const int nodeToSplit = -1);
    int splitAtGivenContourTreeNode(const int labelToSplit, const int mergeAfterNumberOfNodes, SetUnionDataStructure& setUnion);
    /// Updates label field and calls updateInternalStructuresAfterSplit to update remaining structures
    void setTwoRegionsContourTreeSplitResult(SetUnionDataStructure& setUnion, const int labelToSplit);
    void standardContourTreeSplit(const int nodeToSplit = -1);
    /// Updates label field and calls updateInternalStructuresAfterSplit to update remaining structures
    void setStandardContourtreeSplitResult(const int labelToSplit, SetUnionDataStructure& setUnion);
    void sortNeighbors(const McDArray<mculong>& neighbors,
                       const McDArray<float>& values,
                       McDArray<mculong>& sortedNeighborIds,
                       SetUnionDataStructure& setUnion);
    void sortNeighborsAccordingToValues(const McDArray<float>& values,
                                        McDArray<mculong>& sortedNeighborIds);
    void sortNeighborsAccordingToMajorityVote(const McDArray<mculong>& neighbors,
                                              McDArray<mculong>& sortedNeighborIds,
                                              SetUnionDataStructure& setUnion);

    // Spectral clustering split
    /// Creates sparse matrix
    void spectralClusteringSplit(const int nodeToSplit = -1);
    /// Deprecated, creates full matrix
    void spectralClusteringSplitCalculateFullMatrix(const int nodeToSplit = -1);
    /// Updates label field and calls updateInternalStructuresAfterSplit to update remaining structures
    void setSpectralClusteringSplitResult(const int labelToSplit);

    // Thresholding
    /// If splits are used after a thresholding: recompute the distance map
    void thresholdingOnSelectedSegments(float threshold);

    // Surface fields
    void createSurfaceField(QString attributeName);

    /// Create graph for labelfield without empty labels
    void createGraphForCleanLabels();

    /// Remove all labels where all voxels have value 0 in reference label field
    void removeRegionsAccordingToReferenceLabelField();

    // Starting from here: all methods are not used and testcases
    bool floodRegionAndInitForSplit(const int splitRegion); // Not used: test

    void splitAll(); // Not used: test

    void returnLocalLabelfield(); // Not used: test

    void filter(McDArray<mculong>& nodes, QString attributeName);               // Not used: test
    void filterFloat(McDArray<mculong>& nodes, EdgeVertexAttribute* attribute); // Not used: test
    void filterInt(McDArray<mculong>& nodes, EdgeVertexAttribute* attribute);   // Not used: test

    HxSpatialGraph* m_graph;
    HxUniformLabelField3* m_labelField;

    // Mappings between labels in labelfield and nodes/vertices in graph
    McDArray<mclong> m_labelToNode; // without ext-label
    McDArray<mclong> m_nodeToLabel; // without ext-label

    McDArray<mculong> m_criticalNodes;

    McDArray<mculong> m_sizeOfRegion; // with ext-label, contains number of voxels
    McDArray<McVec3f> m_centers;      // with ext-label

    // Data for split methods
    McDArray<mculong> m_nodeIdsOfSplitLabel;
    McDArray<mclong> m_nodeIdToPositionInLabel;
    McHandle<SimplicialMesh3DForHexahedralMesh> m_mesh;
    McDArray<mculong> m_sortedNodeIdx;

    // Matlab engine for spectral clustering split
    Engine* m_eng;

    // Visualization of labelfield is done with voxelized renderer
    // or a surface
    HxColormap256* m_colormap;
    HxVolumeRenderingSettings* m_volumeRenderingSettings;
    HxSurface* m_surface;
    HxGMC* m_surfaceGen;
    HxDisplaySurface* m_surfaceDisplay;

    // Visualization of graph
    HxLineRaycast* m_lineRayCast;

    // m_surfaceNotUpdated and m_autoUpdateSurface control the interaction with the surface
    // surface computation may take a long time so this visualization is optional
    // it is possible to use only VoxelizedRender for visualization
    bool m_surfaceNotUpdated; // Set false only in computeSurface method
                              // Set true whenever labelfield is changed (meld, remove, split, thresholding)
                              // If true everything working on the surface is deactivated (show, some statistics, some critical region computations)
    bool m_autoUpdateSurface; // Set by user
                              // If true the surface is recomputed after all labelfield changes

    bool m_initDone;      // Set true after apply is pressed and module is initialized
                          // Set false if labelfield-port or graph-port is changed
    bool m_initSplitDone; // Set true after m_mesh is initialized for split methods
                          // Set false if distance map-port, labelfield-port or data-port is changed

    enum SelectionType
    {
        SELECTION_SMALL = 0,
        SELECTION_NEAR = 1,
        SELECTION_INSIDE = 2,
        SELECTION_NEIGHBOURS = 3,
        SELECTION_CONNECTED_COMPONENTS = 4,
        SELECTION_EXTERIOR = 5
    };
};
