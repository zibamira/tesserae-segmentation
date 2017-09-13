#include "HxSegmentationEvaluation.h"

#include <hxcore/HxMessage.h>

#include <mclib/internal/McDMatrix.h>

#include <hxcontourtree/SimplicialMesh3DForHexahedralMesh.h>

HX_INIT_CLASS(HxSegmentationEvaluation, HxCompModule)

HxSegmentationEvaluation::HxSegmentationEvaluation()
    : HxCompModule(HxUniformScalarField3::getClassTypeId())
    , portGroundTruthSegmentation(this, "groundTruthSegmentation", tr("Ground Truth Segmentation"), HxUniformScalarField3::getClassTypeId())
    , portMask(this, "mask", tr("Mask"), HxUniformScalarField3::getClassTypeId())
    , portUserSelection(this, "userSelection", tr("User Selection"), 3)
    , portAction(this, "action", tr("Action"))
{
    portUserSelection.setLabel(0, tr("VI"));
    portUserSelection.setLabel(1, tr("RAND"));
    portUserSelection.setLabel(2, tr("Dice"));
    portUserSelection.setValue(0, 1);
    portUserSelection.setValue(1, 0);
    portUserSelection.setValue(2, 0);
}

HxSegmentationEvaluation::~HxSegmentationEvaluation()
{
}

void
HxSegmentationEvaluation::update()
{
}

void
HxSegmentationEvaluation::compute()
{
    if (portAction.wasHit())
    {
        HxUniformScalarField3* segmentationField = hxconnection_cast<HxUniformScalarField3>(portData);
        HxUniformScalarField3* groundTruthField = hxconnection_cast<HxUniformScalarField3>(portGroundTruthSegmentation);
        HxUniformScalarField3* maskField = hxconnection_cast<HxUniformScalarField3>(portMask);
        if (!checkInputFields(segmentationField, groundTruthField, maskField))
        {
            return;
        }
        if (portUserSelection.getValue(0))
        {
            double vi = variationOfInformation(segmentationField, groundTruthField, maskField);
            theMsg->printf("Variation of information: %f", vi);
        }
        if (portUserSelection.getValue(1))
        {
            double rand = randIndex(segmentationField, groundTruthField, maskField);
            theMsg->printf("Rand index: %f", rand);
        }
        if (portUserSelection.getValue(2))
        {
            double dice = diceMeasureOnForeground(segmentationField, groundTruthField, maskField);
            theMsg->printf("Dice : %f", dice);
        }
    }
}

bool
HxSegmentationEvaluation::checkInputFields(HxUniformScalarField3* segmentationField,
                                           HxUniformScalarField3* groundTruthField,
                                           HxUniformScalarField3* maskField)
{
    if (!segmentationField || !groundTruthField || !maskField)
    {
        theMsg->printf("Attach all datasets");
        return false;
    }
    McPrimType dataTypeSegmentation = segmentationField->lattice().primType();
    McPrimType dataTypeSGroundTruth = groundTruthField->lattice().primType();
    if ((dataTypeSegmentation != McPrimType::MC_UINT16) || (dataTypeSGroundTruth != McPrimType::MC_UINT16))
    {
        theMsg->printf("Cast fields to unsigned 16 bit");
        return false;
    }
    return true;
}

double
HxSegmentationEvaluation::diceMeasureOnForeground(HxUniformScalarField3* segmentationField,
                                                  HxUniformScalarField3* groundTruthField,
                                                  HxUniformScalarField3* maskField)
{
    const McDim3l& dims = maskField->lattice().getDims();
    int numberVoxelsInSegmentation = 0;
    int numberVoxelsInGroundTruth = 0;
    int numberVoxelsInSegmentationAndGroundTruth = 0;

    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float inMask;
                maskField->lattice().eval(i, j, k, &inMask);
                if (inMask != 0)
                {
                    float segmentationLabel;
                    segmentationField->lattice().eval(i, j, k, &segmentationLabel);
                    float groundTruthLabel;
                    groundTruthField->lattice().eval(i, j, k, &groundTruthLabel);
                    if (segmentationLabel > 0)
                    {
                        numberVoxelsInSegmentation++;
                    }
                    if (groundTruthLabel > 0)
                    {
                        numberVoxelsInGroundTruth++;
                    }
                    if ((segmentationLabel > 0) && (groundTruthLabel > 0))
                    {
                        numberVoxelsInSegmentationAndGroundTruth++;
                    }
                }
            }
        }
    }

    return (2 * (double)numberVoxelsInSegmentationAndGroundTruth) /
           (numberVoxelsInSegmentation + numberVoxelsInGroundTruth);
}

double
HxSegmentationEvaluation::variationOfInformation(HxUniformScalarField3* segmentationField,
                                                 HxUniformScalarField3* groundTruthField,
                                                 HxUniformScalarField3* maskField)
{
    float minLabel, maxLabel;
    segmentationField->getRange(minLabel, maxLabel);
    const int numberSegmentsInSegmentation = (int)maxLabel + 1;
    groundTruthField->getRange(minLabel, maxLabel);
    const int numberSegmentsInGroundTruth = (int)maxLabel + 1;

    McDArray<int> X;
    X.resize(numberSegmentsInSegmentation);
    X.fill(0);
    McDArray<int> Y;
    Y.resize(numberSegmentsInGroundTruth);
    Y.fill(0);
    McDMatrix<int> A;
    A.resize(numberSegmentsInSegmentation, numberSegmentsInGroundTruth);
    A.fill(0);

    const McDim3l& dims = maskField->lattice().getDims();
    int numberVoxelsInStrip = 0;

    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int k = 0; k < dims[2]; k++)
            {
                float inStrip;
                maskField->lattice().eval(i, j, k, &inStrip);
                if (inStrip != 0)
                {
                    numberVoxelsInStrip++;
                    float segmentationLabel;
                    segmentationField->lattice().eval(i, j, k, &segmentationLabel);
                    float groundTruthLabel;
                    groundTruthField->lattice().eval(i, j, k, &groundTruthLabel);
                    X[(int)segmentationLabel]++;
                    Y[(int)groundTruthLabel]++;
                    A[(int)segmentationLabel][(int)groundTruthLabel]++;
                }
            }
        }
    }

    double h1 = 0, h2 = 0;
    double vi = 0, r = 0, p = 0, q = 0, v = 0;
    for (int i = 0; i < numberSegmentsInSegmentation; i++)
    {
        for (int j = 0; j < numberSegmentsInGroundTruth; j++)
        {
            r = (double)A[i][j] / (double)numberVoxelsInStrip;
            p = (double)X[i] / (double)numberVoxelsInStrip;
            q = (double)Y[j] / (double)numberVoxelsInStrip;

            if ((i == 0) && (q != 0))
            {
                h2 -= q * log(q);
            }

            if (r == 0)
            {
                continue;
            }
            v = r * (log(r / p) + log(r / q));
            // Used to find the worst VI pairs:
            //if (v < -0.004)
            //{
            //    theMsg->printf("v = %f, segment 1 = %i, segment 2 = %i", v, i, j);
            //}
            vi -= v;
        }
        if (p != 0)
        {
            h1 -= p * log(p);
        }
    }
    // Used to print entropies
    //theMsg->printf("H1: %f, H2: %f", h1, h2);

    return vi;
}

double
HxSegmentationEvaluation::randIndex(HxUniformScalarField3* segmentationField,
                                    HxUniformScalarField3* groundTruthField,
                                    HxUniformScalarField3* maskField)
{
    SimplicialMesh3DForHexahedralMesh maskMesh(&(maskField->lattice()));
    maskMesh.setNeighborhood(SimplicialMesh3DForHexahedralMesh::NEIGHBORHOOD_26);
    maskMesh.setThreshold(0.5);

    McDArray<mculong> sortedNodeIdx;
    const bool increasingOrder = false;
    maskMesh.getSortedListOfVertices(increasingOrder, sortedNodeIdx);

    mcuint16* segmentationData = (mcuint16*)segmentationField->lattice().dataPtr();
    mcuint16* groundTruthData = (mcuint16*)groundTruthField->lattice().dataPtr();

    mcuint64 a = 0, b = 0, c = 0, d = 0;

    for (mclong i = 0; i < sortedNodeIdx.size(); ++i)
    {
        const mculong nodeIdx1 = sortedNodeIdx[i];
        const mculong meshVertexId1 = maskMesh.getMeshVertexIdx(nodeIdx1);

        for (mclong j = 0; j < sortedNodeIdx.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }
            const mculong nodeIdx2 = sortedNodeIdx[j];
            const mculong meshVertexId2 = maskMesh.getMeshVertexIdx(nodeIdx2);

            const int segmentationLabel1 = segmentationData[meshVertexId1];
            const int segmentationLabel2 = segmentationData[meshVertexId2];
            const int groundTruthLabel1 = groundTruthData[meshVertexId1];
            const int groundTruthLabel2 = groundTruthData[meshVertexId2];

            if (segmentationLabel1 == segmentationLabel2)
            {
                if (groundTruthLabel1 == groundTruthLabel2)
                {
                    a++;
                }
                else
                {
                    c++;
                }
            }
            else
            {
                if (groundTruthLabel1 == groundTruthLabel2)
                {
                    d++;
                }
                else
                {
                    b++;
                }
            }
        }
    }

    //theMsg->printf("a: %llu, b: %llu, c: %llu, d: %llu", a, b, c, d);
    return ((double)(a + b)) / ((double)(a + b + c + d));
}

int
HxSegmentationEvaluation::parse(Tcl_Interp* t, int argc, char** argv)
{
    if (argc < 2)
        return TCL_OK;
    char* cmd = argv[1];

    if (CMD("getRandIndex"))
    {
        HxUniformScalarField3* segmentationField = hxconnection_cast<HxUniformScalarField3>(portData);
        HxUniformScalarField3* groundTruthField = hxconnection_cast<HxUniformScalarField3>(portGroundTruthSegmentation);
        HxUniformScalarField3* maskField = hxconnection_cast<HxUniformScalarField3>(portMask);
        if (!checkInputFields(segmentationField, groundTruthField, maskField))
        {
            return TCL_ERROR;
        }
        double rand = randIndex(segmentationField, groundTruthField, maskField);
        Tcl_VaSetResult(t, "%f", rand);
    }
    else if (CMD("getVI"))
    {
        HxUniformScalarField3* segmentationField = hxconnection_cast<HxUniformScalarField3>(portData);
        HxUniformScalarField3* groundTruthField = hxconnection_cast<HxUniformScalarField3>(portGroundTruthSegmentation);
        HxUniformScalarField3* maskField = hxconnection_cast<HxUniformScalarField3>(portMask);
        if (!checkInputFields(segmentationField, groundTruthField, maskField))
        {
            return TCL_ERROR;
        }
        double vi = variationOfInformation(segmentationField, groundTruthField, maskField);
        Tcl_VaSetResult(t, "%f", vi);
    }

    else
        return (HxCompModule::parse(t, argc, argv));

    return TCL_OK;
}
