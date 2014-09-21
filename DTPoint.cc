#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

DTPoint::DTPoint()
	: m_dimension(0), m_pCoords(0)
{
}

DTPoint::DTPoint(const Variable* pCoords, uint32_t dimension)
	: m_dimension(dimension)
{
	// no need to initialize m_pCoords to 0 since if a bad_alloc is raised the destructor will not be called.

	m_pCoords = new Variable[m_dimension];
	memcpy(m_pCoords, pCoords, m_dimension * sizeof(Variable));
}

DTPoint::DTPoint(const DTPoint& p)
	: m_dimension(p.m_dimension)
{
	// no need to initialize m_pCoords to 0 since if a bad_alloc is raised the destructor will not be called.

	m_pCoords = new Variable[m_dimension];
	memcpy(m_pCoords, p.m_pCoords, m_dimension * sizeof(Variable));
}

DTPoint::~DTPoint()
{
	delete[] m_pCoords;
}

DTPoint& DTPoint::operator=(const DTPoint& p)
{
	if (this != &p)
	{
		makeDimension(p.m_dimension);
		memcpy(m_pCoords, p.m_pCoords, m_dimension * sizeof(Variable));
	}

	return *this;
}

bool DTPoint::operator==(const DTPoint& p) const
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTPoint::operator==: Points have different number of dimensions."
		);

	for (uint32_t i = 0; i < m_dimension; ++i)
	{
		if (
            m_pCoords[i].type() == VariableType::Continuous && (
			(double)m_pCoords[i] < (double)p.m_pCoords[i] - std::numeric_limits<double>::epsilon() ||
			(double)m_pCoords[i] > (double)p.m_pCoords[i] + std::numeric_limits<double>::epsilon()) ||
            (long long)m_pCoords[i] != (long long)p.m_pCoords[i]
        )
            return false;
	}

	return true;
}

//
// IObject interface
//
DTPoint* DTPoint::clone()
{
	return new DTPoint(*this);
}

//
// ISerializable interface
//
uint32_t DTPoint::getByteArraySize()
{
	return (sizeof(uint32_t) + m_dimension * sizeof(Variable));
}

void DTPoint::loadFromByteArray(const byte* ptr)
{
	uint32_t dimension;
	memcpy(&dimension, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	makeDimension(dimension);
	memcpy(m_pCoords, ptr, m_dimension * sizeof(Variable));
	//ptr += m_dimension * sizeof(Variable);
}

void DTPoint::storeToByteArray(byte** data, uint32_t& len)
{
	len = getByteArraySize();
	*data = new byte[len];
	byte* ptr = *data;

	memcpy(ptr, &m_dimension, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, m_pCoords, m_dimension * sizeof(Variable));
	//ptr += m_dimension * sizeof(Variable);
}

//
// IShape interface
//
bool DTPoint::intersectsShape(const IShape& s) const
{
	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0)
	{
		return pr->containsPoint(*this);
	}

	throw Tools::IllegalStateException(
		"DTPoint::intersectsShape: Not implemented yet!"
	);
}

bool DTPoint::containsShape(const IShape& s) const
{
	return false;
}

bool DTPoint::touchesShape(const IShape& s) const
{
	const DTPoint* ppt = dynamic_cast<const DTPoint*>(&s);
	if (ppt != 0)
	{
		if (*this == *ppt) return true;
		return false;
	}

	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0)
	{
		return pr->touchesPoint(*this);
	}

	throw Tools::IllegalStateException(
		"DTPoint::touchesShape: Not implemented yet!"
	);
}

void DTPoint::getCenter(DTPoint& out) const
{
	out = *this;
}

uint32_t DTPoint::getDimension() const
{
	return m_dimension;
}

void DTPoint::getMBR(Region& out) const
{
	out = Region(m_pCoords, m_pCoords, m_dimension);
}

double DTPoint::getArea() const
{
	return 0.0;
}

double DTPoint::getMinimumDistance(const IShape& s) const
{
	const DTPoint* ppt = dynamic_cast<const DTPoint*>(&s);
	if (ppt != 0)
	{
		return getMinimumDistance(*ppt);
	}

	const Region* pr = dynamic_cast<const Region*>(&s);
	if (pr != 0)
	{
		return pr->getMinimumDistance(*this);
	}

	throw Tools::IllegalStateException(
		"DTPoint::getMinimumDistance: Not implemented yet!"
	);
}

double DTPoint::getMinimumDistance(const DTPoint& p) const
{
	if (m_dimension != p.m_dimension)
		throw Tools::IllegalArgumentException(
			"DTPoint::getMinimumDistance: Shapes have different number of dimensions."
		);

	double ret = 0.0;

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
        if (m_pCoords[cDim].type() == VariableType::Continuous)
            ret += std::pow((long long)m_pCoords[cDim] - (long long)p.m_pCoords[cDim], 2.0);
        else
            ret += std::pow((long long)m_pCoords[cDim] - (long long)p.m_pCoords[cDim], 2.0);
	}

	return std::sqrt(ret);
}

Variable DTPoint::getCoordinate(uint32_t index) const
{
	if (index >= m_dimension)
		throw Tools::IndexOutOfBoundsException(index);

	return m_pCoords[index];
}

void DTPoint::makeInfinite(uint32_t dimension)
{
	makeDimension(dimension);
	for (uint32_t cIndex = 0; cIndex < m_dimension; ++cIndex)
	{
		if (m_pCoords[cIndex].type() == VariableType::Continuous) m_pCoords[cIndex] = std::numeric_limits<double>::max();
        else m_pCoords[cIndex] = std::numeric_limits<long long>::max();
	}
}

void DTPoint::makeDimension(uint32_t dimension)
{
	if (m_dimension != dimension)
	{
		delete[] m_pCoords;

		// remember that this is not a constructor. The object will be destructed normally if
		// something goes wrong (bad_alloc), so we must take care not to leave the object at an intermediate state.
		m_pCoords = 0;

		m_dimension = dimension;
		m_pCoords = new Variable[m_dimension];
	}
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const DTPoint& pt)
{
	for (uint32_t cDim = 0; cDim < pt.m_dimension; ++cDim)
	{
		os << (pt.m_pCoords[cDim].type() == VariableType::Continuous ? (double)pt.m_pCoords[cDim] : (long long)pt.m_pCoords[cDim]) << " ";
	}

	return os;
}
