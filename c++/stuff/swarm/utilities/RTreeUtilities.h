#ifndef __RTREE_UTILITIES_H__
#define __RTREE_UTILITIES_H__

class IntersectingVisitor : public SpatialIndex::IVisitor {
    private: 
        uint64_t intersections_;

    public:

        IntersectingVisitor() : intersections_(0) {};

        inline void prepareForQuery() { intersections_ = 0; }
        inline uint64_t intersections() { return intersections_; }

        void visitData(const SpatialIndex::IData&) { intersections_++; }
        // No interest in these functions.
        void visitNode(const SpatialIndex::INode&) {};
        void visitData(std::vector<const SpatialIndex::IData*>&) {};

};

#endif //__RTREE_UTILITIES_H__
