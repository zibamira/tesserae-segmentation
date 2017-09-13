#include <gtest/internal/gtest.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/HxSpatialGraphViewTest.h>
#include <hxlines/internal/HxLineSet.h>
#include "HxLineRaycast.h"
#include <hxcore/HxObjectPool.h>
#include <hxcolor/HxColormap256.h>

TEST(HxLineRaycastTest, testAdjustRangeSpatialGraph)
{
    McHandle<HxSpatialGraph> graph = HxSpatialGraphViewTest::makeSpatialGraph();

    McHandle<HxLineRaycast> view = HxLineRaycast::createInstance();
    view->portData.connect(graph);
    view->fire();

    McHandle<HxColormap> color(HxColormap256::createInstance());
    EXPECT_TRUE(color);
    color->setMinMax(0.0f, 0.0f);

    view->mPortEndingColorStyle.setValue(1, 2);
    view->fire();
    view->mPortEndingColorMap.connect(color);
    view->fire();
    // Note that 'view->mPortEndingColorMap.adjustRange();' does not adjust the
    // range as one might expect.  Therefore, we provide the data object.
    view->mPortEndingColorMap.adjustRange(graph);
    EXPECT_FLOAT_EQ(1.f, view->mPortEndingColorMap.getMinValue());
    EXPECT_FLOAT_EQ(2.f, view->mPortEndingColorMap.getMaxValue());
    EXPECT_FLOAT_EQ(1.f, color->minCoord());
    EXPECT_FLOAT_EQ(2.f, color->maxCoord());

    view->mPortLineColorStyle.setValue(1, 2);
    view->fire();
    view->mPortLineColorMap.connect(color);
    view->fire();
    view->mPortLineColorMap.adjustRange(graph);
    EXPECT_FLOAT_EQ(4.f, color->minCoord());
    EXPECT_FLOAT_EQ(5.f, color->maxCoord());

    view->mPortLineColorStyle.setValue(1, 3);
    view->fire();
    view->mPortLineColorMap.connect(color);
    view->fire();
    view->mPortLineColorMap.adjustRange(graph);
    EXPECT_FLOAT_EQ(8.f, color->minCoord());
    EXPECT_FLOAT_EQ(9.f, color->maxCoord());

    view->mPortNodeColorStyle.setValue(1, 2);
    view->fire();
    view->mPortNodeColorMap.connect(color);
    view->fire();
    view->mPortNodeColorMap.adjustRange(graph);
    EXPECT_FLOAT_EQ(4.f, color->minCoord());
    EXPECT_FLOAT_EQ(5.f, color->maxCoord());

    view->mPortNodeColorStyle.setValue(1, 3);
    view->fire();
    view->mPortNodeColorMap.connect(color);
    view->fire();
    view->mPortNodeColorMap.adjustRange(graph);
    EXPECT_FLOAT_EQ(8.f, color->minCoord());
    EXPECT_FLOAT_EQ(9.f, color->maxCoord());
}

TEST(HxLineRaycastTest, testAdjustRangeLineSet)
{
    McHandle<HxLineSet> lines = HxLineSet::createInstance();
    const int nPoints = 2;
    McVec3f p[nPoints] = { McVec3f(0, 0, 1.f), McVec3f(0, 0, 2.f) };
    lines->addPoints(p, nPoints);
    lines->setNumDataValues(2);
    lines->data[0][0] = 1.0f;
    lines->data[0][1] = 2.0f;
    lines->data[1][0] = 3.0f;
    lines->data[1][1] = 4.0f;

    McHandle<HxLineRaycast> view = HxLineRaycast::createInstance();
    view->portData.connect(lines);
    view->fire();

    McHandle<HxColormap> color(HxColormap256::createInstance());
    EXPECT_TRUE(color);

    view->mPortEndingColorStyle.setValue(1, 2);
    view->fire();
    view->mPortEndingColorMap.connect(color);
    view->fire();
    color->setMinMax(0.0f, 0.0f);
    view->mPortEndingColorMap.adjustRange(lines);
    EXPECT_FLOAT_EQ(1.0f, color->minCoord());
    EXPECT_FLOAT_EQ(2.0f ,color->maxCoord());

    view->mPortEndingColorStyle.setValue(1, 3);
    view->fire();
    view->mPortEndingColorMap.connect(color);
    view->fire();
    color->setMinMax(0.0f, 0.0f);
    view->mPortEndingColorMap.adjustRange(lines);
    EXPECT_FLOAT_EQ(3.0f, color->minCoord());
    EXPECT_FLOAT_EQ(4.0f, color->maxCoord());
}
