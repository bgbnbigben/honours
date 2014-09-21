#include <spatialindex/SpatialIndex.h>

#include <cstring>
#include <cmath>
#include <limits>

DTRegion::DTRegion()
	: m_dimension(0), m_pLow(0), m_pHigh(0)
{
}

DTRegion::DTRegion(const Variable* pLow, const Variable* pHigh, uint32_t dimension)
{
	initialize(pLow, pHigh, dimension);
}

DTRegion::Region(const Point& low, const Point& high)
{
	if (low.m_dimension != high.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::Region: arguments have different number of dimensions."
		);

	initialize(low.m_pCoords, high.m_pCoords, low.m_dimension);
}

DTRegion::Region(const Region& r)
{
	initialize(r.m_pLow, r.m_pHigh, r.m_dimension);
}

void DTRegion::initialize(const Variable* pLow, const Variable* pHigh, uint32_t dimension)
{
	m_pLow = 0;
	m_dimension = dimension;

#ifndef NDEBUG
    for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
    {
     if ((pLow[cdim].type() == VariableType::Continuous && ((double)pLow[cDim] > (double)pHigh[cDim])) || (long long)pLow[cDim] > (long long)pHigh[cDim])
     {
         // check for infinitive region
         if ((pLow[cDim].type() == VariableType::Continuous && !(pLow[cDim] == std::numeric_limits<double>::max() ||
             pHigh[cDim] == -std::numeric_limits<double>::max() )) ||
                 !(pLow[cDim] == std::numeric_limits<long long>::max() || pHigh[cDim] == -std::numeric_limit<long long>::max()))
             throw Tools::IllegalArgumentException(
                 "DTRegion::initialize: Low point has larger coordinates than High point."
                 " Neither point is infinity."
             );
     }
    }
#endif

	try
	{
		m_pLow = new Variable[m_dimension];
		m_pHigh = new Variable[m_dimension];
	}
	catch (...)
	{
		delete[] m_pLow;
		throw;
	}

    // TODO deal with the memcpy
	memcpy(m_pLow, pLow, m_dimension * sizeof(Variable));
	memcpy(m_pHigh, pHigh, m_dimension * sizeof(Variable));
}

DTRegion::~Region()
{
	delete[] m_pLow;
	delete[] m_pHigh;
}

Region& DTRegion::operator=(const Region& r)
{
	if(this != &r)
	{
		makeDimension(r.m_dimension);
		memcpy(m_pLow, r.m_pLow, m_dimension * sizeof(Variable));
		memcpy(m_pHigh, r.m_pHigh, m_dimension * sizeof(Variable));
	}

	return *this;
}

bool DTRegion::operator==(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::operator==: Regions have different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (
			m_pLow[i] < r.m_pLow[i] - std::numeric_limits<double>::epsilon() ||
			m_pLow[i] > r.m_pLow[i] + std::numeric_limits<double>::epsilon() ||
			m_pHigh[i] < r.m_pHigh[i] - std::numeric_limits<double>::epsilon() ||
			m_pHigh[i] > r.m_pHigh[i] + std::numeric_limits<double>::epsilon())
			return false;
	}
	return true;
}

//
// IObject interface
//
Region* DTRegion::clone()
{
	return new Region(*this);
}

//
// ISerializable interface
//
uint32_t DTRegion::getByteArraySize()
{
	return (sizeof(uint32_t) + 2 * m_dimension * sizeof(Variable));
}

void DTRegion::loadFromByteArray(const byte* ptr)
{
	uint32_t dimension;
	memcpy(&dimension, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	makeDimension(dimension);
	memcpy(m_pLow, ptr, m_dimension * sizeof(Variable));
	ptr += m_dimension * sizeof(Variable);
	memcpy(m_pHigh, ptr, m_dimension * sizeof(Variable));
	//ptr += m_dimension * sizeof(Variable);
}

void DTRegion::storeToByteArray(byte** data, uint32_t& len)
{
	len = getByteArraySize();
	*data = new byte[len];
	byte* ptr = *data;

	memcpy(ptr, &m_dimension, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, m_pLow, m_dimension * sizeof(Variable));
	ptr += m_dimension * sizeof(Variable);
	memcpy(ptr, m_pHigh, m_dimension * sizeof(Variable));
	//ptr += m_dimension * sizeof(Variable);
}

//
// IShape interface
//
bool DTRegion::intersectsShape(const IShape& s) const
{
	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0) return intersectsRegion(*pr);

	const LineSegment* pls = dynamic_cast<const LineSegment*>(&s);
	if (pls != 0) return intersectsLineSegment(*pls);

	const Point* ppt = dynamic_cast<const Point*>(&s);
	if (ppt != 0) return containsPoint(*ppt);

	throw Tools::IllegalStateException(
		"DTRegion::intersectsShape: Not implemented yet!"
	);
}

bool DTRegion::containsShape(const IShape& s) const
{
	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0) return containsRegion(*pr);

	const Point* ppt = dynamic_cast<const Point*>(&s);
	if (ppt != 0) return containsPoint(*ppt);

	throw Tools::IllegalStateException(
		"DTRegion::containsShape: Not implemented yet!"
	);
}

bool DTRegion::touchesShape(const IShape& s) const
{
	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0) return touchesRegion(*pr);

	const Point* ppt = dynamic_cast<const Point*>(&s);
	if (ppt != 0) return touchesPoint(*ppt);

	throw Tools::IllegalStateException(
		"DTRegion::touchesShape: Not implemented yet!"
	);
}

void DTRegion::getCenter(Point& out) const
{
	out.makeDimension(m_dimension);
	for (uint32_t i = 0; i < m_dimension; ++i)
	{
        if (m_pLow[i].type() == VariableType::Continuous)
            out.m_pCoords[i] = (double)m_pLow[i] + ((double)m_pHigh[i] - (double)m_pLow[i]) / 2.0;
        else
            out.m_pCoords[i] = (long long)m_pLow[i] + ((long long)m_pHigh[i] - (long long)m_pLow[i]) / 2.0;
	}
}

uint32_t DTRegion::getDimension() const
{
	return m_dimension;
}

void DTRegion::getMBR(Region& out) const
{
	out = *this;
}

double DTRegion::getArea() const
{
	double area = 1.0;

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
        if (m_pLow[i].type() == VariableType::Continuous)
            area *= (double)m_pHigh[i] - (double)m_pLow[i];
        else
            // TODO this could be garbage.
            area *= std::stod(std::to_string((long long)m_pHigh[i] - (long long)m_pLow[i]));

	}

	return area;
}

double DTRegion::getMinimumDistance(const IShape& s) const
{
	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0) return getMinimumDistance(*pr);

	const Point* ppt = dynamic_cast<const Point*>(&s);
	if (ppt != 0) return getMinimumDistance(*ppt);

	throw Tools::IllegalStateException(
		"DTRegion::getMinimumDistance: Not implemented yet!"
	);
}

bool DTRegion::intersectsRegion(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::intersectsRegion: Regions have different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if ((m_pLow[i].type() == VariableType::Continuous && ((double)m_pLow[i] > (double)r.m_pHigh[i] || (double)m_pHigh[i] < (double)r.m_pLow[i])) || ((long long)m_pLow[i] > (long long)r.m_pHigh[i] || (long long)m_pHigh[i] < (long long)r.m_pLow[i])) return false;
	}
	return true;
}

bool DTRegion::containsRegion(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::containsRegion: Regions have different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (m_pLow[i].type() == VariableType::Continuous && ((double)m_pLow[i] > (double)r.m_pLow[i] || (double)m_pHigh[i] < (double)r.m_pHigh[i]) || ((long long)m_pLow[i] > (long long)r.m_pLow[i] || (long long)m_pHigh[i] < (long long)r.m_pHigh[i])) return false;
	}
	return true;
}

bool DTRegion::touchesRegion(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::touchesRegion: Regions have different number of dimensions."
		);
	
	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (
			m_pLow[i].type() == VariableType::Continuous && (((double)m_pLow[i] >= (double)r.m_pLow[i] + std::numeric_limits<double>::epsilon() &&
			(double)m_pLow[i] <= (double)r.m_pLow[i] - std::numeric_limits<double>::epsilon()) ||
			((double)m_pHigh[i] >= (double)r.m_pHigh[i] + std::numeric_limits<double>::epsilon() &&
			(double)m_pHigh[i] <= (double)r.m_pHigh[i] - std::numeric_limits<double>::epsilon())) ||
            (((long long)m_pLow[i] != (long long)r.m_pLow[i]) ||
			((long long)m_pHigh[i] != (long long)r.m_pHigh[i]))
        )
			return false;
	}
	return true;
}


double DTRegion::getMinimumDistance(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::getMinimumDistance: Regions have different number of dimensions."
		);

	double ret = 0.0;

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		double x = 0.0;

        if (m_pLow[i].type() == VariableType::Continuous) {
            if ((double)r.m_pHigh[i] < (double)m_pLow[i])
            {
                x = std::abs((double)r.m_pHigh[i] - (double)m_pLow[i]);
            }
            else if ((double)m_pHigh[i] < (double)r.m_pLow[i])
            {
                x = std::abs((double)r.m_pLow[i] - (double)m_pHigh[i]);
            }
        } else {
            if ((long long)r.m_pHigh[i] < (long long)m_pLow[i])
            {
                x = std::stod(std::tostring(std::abs((long long)r.m_pHigh[i] - (long long)m_pLow[i])));
            }
            else if ((long long)m_pHigh[i] < (long long)r.m_pLow[i])
            {
                x = std::stod(std::to_string(std::abs((long long)r.m_pLow[i] - (long long)m_pHigh[i])));
            }
        }

		ret += x * x;
	}

	return std::sqrt(ret);
}

bool DTRegion::intersectsLineSegment(const LineSegment& in) const
{
    throw Tools::NotSupportedException("This function is not supported for a dual-type region");
	if (m_dimension != 2)
		throw Tools::NotSupportedException(
			"DTRegion::intersectsLineSegment: only supported for 2 dimensions"
		);

	if (m_dimension != in.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::intersectsRegion: Region and LineSegment have different number of dimensions."
		);

    // there may be a more efficient method, but this suffices for now
    /*
    Point ll = Point(m_pLow, 2);
    Point ur = Point(m_pHigh, 2);
    // fabricate ul and lr coordinates and points
    double c_ul[2] = {m_pLow[0], m_pHigh[1]};
    double c_lr[2] = {m_pHigh[0], m_pLow[1]};
    Point ul = Point(&c_ul[0], 2);
    Point lr = Point(&c_lr[0], 2);

    // Points/LineSegment for the segment
    Point p1 = Point(in.m_pStartPoint, 2);
    Point p2 = Point(in.m_pEndPoint, 2);
    

    //Check whether either or both the endpoints are within the region OR
    //whether any of the bounding segments of the Region intersect the segment
    return (containsPoint(p1) || containsPoint(p2) || 
            in.intersectsShape(LineSegment(ll, ul)) || in.intersectsShape(LineSegment(ul, ur)) ||
            in.intersectsShape(LineSegment(ur, lr)) || in.intersectsShape(LineSegment(lr, ll)));
    */
}

bool DTRegion::containsPoint(const Point& p) const
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::containsPoint: Point has different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (m_pLow[i].type() == VariableType::Continuous && ((double)m_pLow[i] > (double)p.getCoordinate(i) || (double)m_pHigh[i] < (double)p.getCoordinate(i)) || ((long long)m_pLow[i] > (long long)p.getCoordinate(i) || (long long)m_pHigh[i] < (long long)p.getCoordinate(i))) return false;
	}
	return true;
}

bool DTRegion::touchesPoint(const Point& p) const
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::touchesPoint: Point has different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (
            m_pLow[i].type() == VariableType::Continuous &&
			(((double)m_pLow[i] >= (double)p.getCoordinate(i) - std::numeric_limits<double>::epsilon() &&
			 (double)m_pLow[i] <= (double)p.getCoordinate(i) + std::numeric_limits<double>::epsilon()) ||
			((double)m_pHigh[i] >= (double)p.getCoordinate(i) - std::numeric_limits<double>::epsilon() &&
			 (double)m_pHigh[i] <= (double)p.getCoordinate(i) + std::numeric_limits<double>::epsilon())) ||
            ((long long)m_pLow[i] == (long long)p.getCoordinate(i) || (long long)m_pHigh[i] == (long long)p.getCoordinate(i))
        )
			return true;
	}
	return false;
}

double DTRegion::getMinimumDistance(const Point& p) const
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::getMinimumDistance: Point has different number of dimensions."
		);

	double ret = 0.0;

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
        if (m_pLow[i].type() == VariableType::Continuous) {
            if ((double)p.getCoordinate(i) < (double)m_pLow[i])
            {
                ret += std::pow((double)m_pLow[i] - (double)p.getCoordinate(i), 2.0);
            }
            else if ((double)p.getCoordinate(i) > (double)m_pHigh[i])
            {
                ret += std::pow((double)p.getCoordinate(i) - (double)m_pHigh[i], 2.0);
            }
        } else {
            // TODO think about manually casting the ((long long), (double)) ll
            // argument in case of some weird overflow/breakage.
            if ((long long)p.getCoordinate(i) < (long long)m_pLow[i])
            {
                ret += std::pow((long long)m_pLow[i] - (long long)p.getCoordinate(i), 2.0);
            }
            else if ((long long)p.getCoordinate(i) > (long long)m_pHigh[i])
            {
                ret += std::pow((long long)p.getCoordinate(i) - (long long)m_pHigh[i], 2.0);
            }
        }
    }

	return std::sqrt(ret);
}

Region DTRegion::getIntersectingRegion(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::getIntersectingRegion: Regions have different number of dimensions."
		);

	Region ret;
	ret.makeInfinite(m_dimension);

	// check for intersection.
	// marioh: avoid function call since this is called billions of times.
	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
		if (m_pLow[cDim].type() == VariableType::Continuous && ((double)m_pLow[cDim] > (double)r.m_pHigh[cDim] || (double)m_pHigh[cDim] < (double)r.m_pLow[cDim]) || ((long long)m_pLow[cDim] > (long long)r.m_pHigh[cDim] || (long long)m_pHigh[cDim] < (long long)r.m_pLow[cDim])) return ret;
	}

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
        if (m_pLow[cDim].type() == VariableType::Continuous) {
            ret.m_pLow[cDim] = std::max((double)m_pLow[cDim], (double)r.m_pLow[cDim]);
            ret.m_pHigh[cDim] = std::min((double)m_pHigh[cDim], (double)r.m_pHigh[cDim]);
        } else {
            ret.m_pLow[cDim] = std::max((long long)m_pLow[cDim], (long long)r.m_pLow[cDim]);
            ret.m_pHigh[cDim] = std::min((long long)m_pHigh[cDim], (long long)r.m_pHigh[cDim]);
        }
	}

	return ret;
}

double DTRegion::getIntersectingArea(const Region& r) const
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::getIntersectingArea: Regions have different number of dimensions."
		);

	double ret = 1.0;
	double f1, f2;

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
		if (m_pLow[cDim].type() == VariableType::Continuous && ((double)m_pLow[cDim] > (double)r.m_pHigh[cDim] || (double)m_pHigh[cDim] < (double)r.m_pLow[cDim]) || ((long long)m_pLow[cDim] > (long long)r.m_pHigh[cDim] || (long long)m_pHigh[cDim] < (long long)r.m_pLow[cDim])) return 0.0;

        if (m_pLow[cDim].type() == VariableType::Continuous) {
            f1 = std::max(m_pLow[cDim], r.m_pLow[cDim]);
            f2 = std::min(m_pHigh[cDim], r.m_pHigh[cDim]);
        } else {
            f1 = std::stod(std::to_string(std::max((long long)m_pLow[cDim], (long long)r.m_pLow[cDim])));
            f2 = std::stod(std::to_string(std::min((long long)m_pHigh[cDim], (long long)r.m_pHigh[cDim])));
        }
		ret *= f2 - f1;
	}

	return ret;
}

/*
 * Returns the margin of a region. It is calcuated as the sum of  2^(d-1) * width, in each dimension.
 * It is actually the sum of all edges, no matter what the dimensionality is.
*/
double DTRegion::getMargin() const
{
	double mul = std::pow(2.0, static_cast<double>(m_dimension) - 1.0);
	double margin = 0.0;

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (m_pLow[i].type() == VariableType::Continuous) margin += ((double)m_pHigh[i] - (double)m_pLow[i]) * mul;
		else margin += ((long long)m_pHigh[i] - (long long)m_pLow[i]) * mul;
	}

	return margin;
}

void DTRegion::combineRegion(const Region& r)
{
	if (m_dimension != r.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::combineRegion: Region has different number of dimensions."
		);

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
        if (m_pLow[cDim].type() == VariableType::Continuous) {
            m_pLow[cDim] = std::min((double)m_pLow[cDim], (double)r.m_pLow[cDim]);
            m_pHigh[cDim] = std::max((double)m_pHigh[cDim], (double)r.m_pHigh[cDim]);
        } else {
            m_pLow[cDim] = std::min((long long)m_pLow[cDim], (long long)r.m_pLow[cDim]);
            m_pHigh[cDim] = std::max((long long)m_pHigh[cDim], (long long)r.m_pHigh[cDim]);
        }
	}
}

void DTRegion::combinePoint(const Point& p)
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::combinePoint: Point has different number of dimensions."
		);

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
        if (m_pLow[cDim].type() == VariableType) {
            m_pLow[cDim] = std::min((double)m_pLow[cDim], (double)p.m_pCoords[cDim]);
            m_pHigh[cDim] = std::max((double)m_pHigh[cDim], (double)p.m_pCoords[cDim]);
        } else {
            m_pLow[cDim] = std::min((long long)m_pLow[cDim], (long long)p.m_pCoords[cDim]);
            m_pHigh[cDim] = std::max((long long)m_pHigh[cDim], (long long)p.m_pCoords[cDim]);
        }
	}
}

void DTRegion::getCombinedRegion(Region& out, const Region& in) const
{
	if (m_dimension != in.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTRegion::getCombinedRegion: Regions have different number of dimensions."
		);

	out = *this;
	out.combineRegion(in);
}

Variable DTRegion::getLow(uint32_t index) const
{
	if (index >= m_dimension)
		throw Tools::IndexOutOfBoundsException(index);

	return m_pLow[index];
}

Variable DTRegion::getHigh(uint32_t index) const
{
	if (index >= m_dimension)
		throw Tools::IndexOutOfBoundsException(index);

	return m_pHigh[index];
}

void DTRegion::makeInfinite(uint32_t dimension)
{
	makeDimension(dimension);
	for (uint32_t cIndex = 0; cIndex < m_dimension; ++cIndex)
	{
        if (m_pLow[cIndex].type() == VariableType::Continuous) {
            m_pLow[cIndex] = std::numeric_limits<double>::max();
            m_pHigh[cIndex] = -std::numeric_limits<double>::max();
        } else {
            m_pLow[cIndex] = std::numeric_limits<long long>::max();
            m_pHigh[cIndex] = -std::numeric_limits<long long>::max();
        }
	}
}

void DTRegion::makeDimension(uint32_t dimension)
{
	if (m_dimension != dimension)
	{
		delete[] m_pLow;
		delete[] m_pHigh;

		// remember that this is not a constructor. The object will be destructed normally if
		// something goes wrong (bad_alloc), so we must take care not to leave the object at an intermediate state.
		m_pLow = 0; m_pHigh = 0;

		m_dimension = dimension;
		m_pLow = new Variable[m_dimension];
		m_pHigh = new Variable[m_dimension];
	}
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const Region& r)
{
	uint32_t i;

	os << "Low: ";
	for (i = 0; i < r.m_dimension; ++i)
	{
		os << (r.m_pLow[i].type() == VariableType::Continuous ? (double)r.m_pLow[i] : (long long)r.m_pLow[i]) << " ";
	}

	os << ", High: ";

	for (i = 0; i < r.m_dimension; ++i)
	{
		os << (r.m_pHigh[i].type() == VariableType::Continuous ? (double)r.m_pHigh[i] : (long long)r.m_pHigh[i]) << " ";
	}

	return os;
}
