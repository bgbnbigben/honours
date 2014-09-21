#pragma once

//class SIDX_DLL Region : public Tools::IObject, public virtual IShape
/* A dual-type'd region */
class DTRegion : public SpatialIndex::Region
{
public:
    DTRegion();
    DTRegion(const Variable* pLow, const Variable* pHigh, uint32_t dimension);
    DTRegion(const Point& low, const Point& high);
    DTRegion(const DTRegion& in);
    virtual ~DTRegion();

    virtual DTRegion& operator=(const DTRegion& r);
    virtual bool operator==(const DTRegion&) const;

    //
    // IObject interface
    //
    virtual DTRegion* clone();

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
    virtual void getCenter(Point& out) const;
    virtual uint32_t getDimension() const;
    virtual void getMBR(DTRegion& out) const;
    virtual double getArea() const;
    virtual double getMinimumDistance(const IShape& in) const;

    virtual bool intersectsRegion(const DTRegion& in) const;
    virtual bool containsRegion(const DTRegion& in) const;
    virtual bool touchesRegion(const DTRegion& in) const;
    virtual double getMinimumDistance(const DTRegion& in) const;

    virtual bool intersectsLineSegment(const LineSegment& in) const;

    virtual bool containsPoint(const Point& in) const;
    virtual bool touchesPoint(const Point& in) const;
    virtual double getMinimumDistance(const Point& in) const;

    virtual DTRegion getIntersectingRegion(const DTRegion& r) const;
    virtual double getIntersectingArea(const DTRegion& in) const;
    virtual double getMargin() const;

    virtual void combineRegion(const DTRegion& in);
    virtual void combinePoint(const Point& in);
    virtual void getCombinedRegion(DTRegion& out, const DTRegion& in) const;

    virtual Variable getLow(uint32_t index) const;
    virtual Variable getHigh(uint32_t index) const;

    virtual void makeInfinite(uint32_t dimension);
    virtual void makeDimension(uint32_t dimension);

private:
    void initialize(const Variable* pLow, const Variable* pHigh, uint32_t dimension);

public:
    //(TODO): Time for some blind hope
    uint32_t m_dimension;
    Variable* m_pLow;
    Variable* m_pHigh;

    friend std::ostream& operator<<(std::ostream& os, const DTRegion& r);
}; // DTRegion

typedef Tools::PoolPointer<DTRegion> DTRegionPtr;
std::ostream& operator<<(std::ostream& os, const DTRegion& r);
