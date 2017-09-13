#include "SetUnionDataStructure.h"

#include <mclib/internal/McAssert.h>

void
SetUnionDataStructure::setNumElements(
    const mclong numElements)
{
    m_setIds.resize(numElements);
    m_setIds.fill(-1);
}

mculong
SetUnionDataStructure::getNumElements() const
{
    return m_setIds.size();
}

void
SetUnionDataStructure::setSetIdOfElement(
    const mclong elemId,
    const mclong setId)
{
    mcassert(m_setIds.size() > elemId);
    mcassert(m_setIds.size() > setId);

    m_setIds[elemId] = setId;
}

void
SetUnionDataStructure::mergeSetsOfElements(
    const mculong elem1Id,
    const mculong elem2Id)
{
    const mculong oldSetId = findSetId(elem1Id);
    const mculong newSetId = findSetId(elem2Id);

    m_setIds[elem1Id] = newSetId;

    mergeSets(oldSetId, newSetId);
}

void
SetUnionDataStructure::mergeSets(
    const mculong setId,
    const mculong newSetId)
{
    const mculong oldSetId = m_setIds[setId];

    m_setIds[setId] = newSetId;

    if ( setId == oldSetId )
        return;

    mergeSets(oldSetId, newSetId);
}

mculong
SetUnionDataStructure::findSetId(
    const mculong elemId)
{
    mclong eId   = elemId;
    mclong setId = m_setIds[eId];

    while ( m_setIds[eId] != eId )
    {
        m_setIds[eId] = m_setIds[m_setIds[setId]];
        eId           = m_setIds[eId];
        setId         = m_setIds[eId];
    }

    m_setIds[elemId] = setId;

    return setId;
}
