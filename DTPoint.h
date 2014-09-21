#pragma once

#include "tools/Tools.h"

class DTPoint : public Tools::IObject, public virtual IShape
{
public:
    DTPoint();
    DTPoint(const Variable* pCoords, uint32_t dimension);
    DTPoint(const DTPoint& p);
    virtual ~DTPoint();

    virtual DTPoint& operator=(const DTPoint& p);
    virtual bool operator==(const DTPoint& p) const;

    //
    // IObject interface
    //
    virtual DTPoint* clone();

    //
    // ISerializable interface
    //
    virtual uint32_t getByteArraySize();
    virtual void loadFromByteArray(const byte* data);
    virtual void storeToByteArray(byte** data, uint32_t& length);

    //
    // IShape interface
    //
    virtual bool intersectsShape(const IShape& in) const;
    virtual bool containsShape(const IShape& in) const;
    virtual bool touchesShape(const IShape& in) const;
    virtual void getCenter(DTPoint& out) const;
    virtual uint32_t getDimension() const;
    virtual void getMBR(Region& out) const;
    virtual double getArea() const;
    virtual double getMinimumDistance(const IShape& in) const;

    virtual double getMinimumDistance(const DTPoint& p) const;

    virtual Variable getCoordinate(uint32_t index) const;

    virtual void makeInfinite(uint32_t dimension);
    virtual void makeDimension(uint32_t dimension);

public:
    uint32_t m_dimension;
    Variable* m_pCoords;

    friend class DTRegion;
    friend std::ostream& operator<<(std::ostream& os, const DTPoint& pt);
}; // DTPoint

typedef Tools::PoolPointer<DTPoint> DTPointPtr;

std::ostream& operator<<(std::ostream& os, const DTPoint& pt);
