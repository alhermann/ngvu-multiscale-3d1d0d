// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/
#include <memory>
#include <algorithm>

#include <dune/grid/common/intersection.hh>

#include <dune/foamgrid/foamgrid/foamgridintersectioniterators.hh>
#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridgeometry.hh>
#include <dune/foamgrid/foamgrid/foamgridnulliteratorfactory.hh>

namespace Dune {

template <class GridImp>
class FoamGridLevelIntersectionIterator;

template <class GridImp>
class FoamGridLeafIntersectionIterator;

template <class GridImp>
class FoamGridLevelIntersection;

template <class GridImp>
class FoamGridLeafIntersection;


//! \brief Base class of all intersections within FoamGrid
//!
//! encapsulates common functionality of level and leaf intersections.
template<class GridImp>
class FoamGridIntersection
{
    enum {dimgrid  = GridImp::dimension};
    enum {dimworld = GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;

    // The iterators need access to all members
    friend class FoamGridLevelIntersectionIterator<GridImp>;
    friend class FoamGridLeafIntersectionIterator<GridImp>;

    // As weĺl as the implementations
    friend class FoamGridLevelIntersection<GridImp>;
    friend class FoamGridLeafIntersection<GridImp>;

    template<typename, typename>
    friend class Dune::Intersection;

    FoamGridIntersection()
      : center_(nullptr),
        neighbor_(FoamGridNullIteratorFactory<dimgrid, dimworld>::null())
    {}

    /**
     * \brief Initalizes an intersection.
     *
     * After initialization this object always represents the first intersection
     * related to an facet.
     * \param facet The index of the facet this intersection lives on.
     */
    FoamGridIntersection(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center,
                         int facet)
        : center_(center), facetIndex_(facet), neighbor_(FoamGridNullIteratorFactory<dimgrid, dimworld>::null())
    {}

public:

    typedef typename GridImp::template Codim<0>::Entity Entity;

        //! return Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        Entity inside() const {
            return Entity(FoamGridEntity<0, dimgrid, GridImp> (center_));
        }


        //! return Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        Entity outside(/*std::size_t neighborIndex = 0*/) const {
            // Return the 'other' element on the current facet
            return Entity(FoamGridEntity<0, dimgrid, GridImp> ((*neighbor_)));
        }

        //! equality
        bool equals(const FoamGridIntersection<GridImp>& i) const {
            return center_==i.center_ && neighbor_ == i.neighbor_ && facetIndex_ == i.facetIndex_;
        }

        /** \brief return true if intersection is with boundary.
        */
        bool boundary () const {
           return center_->facet_[facetIndex_]->elements_.size()==1;
        }

        //! return information about the Boundary
        int boundarySegmentIndex () const
        {
            return center_->facet_[facetIndex_]->boundarySegmentIndex();
        }

        //! Geometry type of an intersection
        GeometryType type () const
        {
            return GeometryType(GeometryType::simplex, dimgrid-1);
        }

        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const
        {
            assert(facetIndex_ < center_->corners());
            return facetIndex_;
        }

        virtual int indexInOutside(std::size_t neighborIndex = 0) const=0;

        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dimgrid-1>& local) const
        {
            // The intersection normal is a vector that is orthogonal to the element normal
            // and to the intersection itself.

            // only works for triangles and lines
            assert(center_->type().isTriangle() || center_->type().isLine());

            // Compute vertices
            const Dune::ReferenceElement<double, dimgrid>& refElement
                = Dune::ReferenceElements<double, dimgrid>::general(center_->type());

            // facet vertices (vertex in 1d, edge in 2d)
            int v0 = refElement.subEntity(facetIndex_, 1, 0, dimgrid);

            if(dimgrid == 2) {
                //second vertex only exists for an edge
                int v1 = refElement.subEntity(facetIndex_, 1, 1, dimgrid);
                // opposite vertex
                int v2 = (v1+1)%3;
                if (v2==v0)
                    v2 = (v0+1)%3;
                assert(v2!=v0 and v2!=v1);

                FieldVector<ctype, dimworld> facet = center_->vertex_[v0]->pos_ - center_->vertex_[v1]->pos_;
                FieldVector<ctype, dimworld> otherEdge = center_->vertex_[v2]->pos_ - center_->vertex_[v1]->pos_;

                //Cross product of facet and otherEdge is a scaled element normal
                FieldVector<ctype, dimworld> scaledElementNormal;

                if(dimworld == 3) //dimworld==3
                {
                    scaledElementNormal[0] = facet[1]*otherEdge[2] - facet[2]*otherEdge[1];
                    scaledElementNormal[1] = facet[2]*otherEdge[0] - facet[0]*otherEdge[2];
                    scaledElementNormal[2] = facet[0]*otherEdge[1] - facet[1]*otherEdge[0];
                    outerNormal_[0] = facet[1]*scaledElementNormal[2] - facet[2]*scaledElementNormal[1];
                    outerNormal_[1] = facet[2]*scaledElementNormal[0] - facet[0]*scaledElementNormal[2];
                    outerNormal_[2] = facet[0]*scaledElementNormal[1] - facet[1]*scaledElementNormal[0];
                }
                else //dimworld==2
                {
                    outerNormal_[0] = facet[1];
                    outerNormal_[1] = -facet[0];
                }

                //Check if scaled EdgeNormal is inner normal, if yes flip
                otherEdge = center_->vertex_[v0]->pos_ - center_->vertex_[v2]->pos_;
                if(otherEdge*outerNormal_ < 0)
                    outerNormal_ *= -1.0;

                return outerNormal_;
            }
            else //dimgrid ==1
            {
                int v1 = 1-v0;
                assert(v0!=v1);
                //automatically has right orientation
                outerNormal_ = center_->vertex_[v0]->pos_ - center_->vertex_[v1]->pos_;

                return outerNormal_;
            }
            DUNE_THROW(GridError, "Non-existing grid dimension requested! Has to be 1 or 2!");
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dimgrid-1>& local) const
        {

            if(dimgrid == 2)
            {
                const Dune::ReferenceElement<double, dimgrid>& refElement
                    = Dune::ReferenceElements<double, dimgrid>::general(center_->type());

                // facet vertices
                int v0 = refElement.subEntity(facetIndex_, 1, 0, dimgrid);
                int v1 = refElement.subEntity(facetIndex_, 1, 1, dimgrid);

                //facet length
                ctype facetLength = (center_->vertex_[v0]->pos_- center_->vertex_[v1]->pos_).two_norm();

                FieldVector<ctype, dimworld> integrationOuterNormal_ = unitOuterNormal(local);
                integrationOuterNormal_ *= facetLength;
                return integrationOuterNormal_;
            }
            else //dimgrid == 1
            {
                integrationOuterNormal_ = this->unitOuterNormal(local);
                return integrationOuterNormal_;
            }
            DUNE_THROW(GridError, "Non-existing grid dimension requested! Has to be 1 or 2!");
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dimgrid-1>& local) const
        {
            unitOuterNormal_ = this->outerNormal(local);
            unitOuterNormal_ /= unitOuterNormal_.two_norm();
            return unitOuterNormal_;
        }

        //! return unit outer normal at the intersection center
        FieldVector<ctype, dimworld> centerUnitOuterNormal () const
        {
            return unitOuterNormal(FieldVector<ctype,dimgrid-1>(0.5));
        }
    private:

        //! vector storing the outer normal
        mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;
        mutable FieldVector<typename GridImp::ctype, dimworld> unitOuterNormal_;
        mutable FieldVector<typename GridImp::ctype, dimworld> integrationOuterNormal_;

    protected:

        const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center_;

        /** \brief Count on which facet we are lookin' at.  */
        int facetIndex_;

        /** \brief Iterator to the other neighbor of the intersection. */
        typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::const_iterator neighbor_;

};

template<class GridImp>
class FoamGridLevelIntersection
    : public FoamGridIntersection<GridImp>
{
    friend class FoamGridLevelIntersectionIterator<GridImp>;

    template<typename, typename>
    friend class Dune::Intersection;

    FoamGridLevelIntersection() : FoamGridIntersection<GridImp>() {}

    enum {dimgrid = GridImp::dimension};
    enum{ dimworld = GridImp::dimensionworld };

    // Geometry is a CachedMultiLinearGeometry
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

    FoamGridLevelIntersection(const FoamGridEntityImp<dimgrid, dimgrid ,dimworld>* center, std::size_t facet)
                              : FoamGridIntersection<GridImp>(center, facet)
    {}

public:

    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid level intersections are always conforming
        return true;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside (std::size_t neighborIndex = 0) const {
        return std::find((*this->neighbor_)->facet_.begin(), (*this->neighbor_)->facet_.end(),
                         this->center_->facet_[this->facetIndex_])
            - (*this->neighbor_)->facet_.begin();
    }

    //! return true if across the facet a neighbor on this level exists
    bool neighbor () const {
      return this->neighbor_!=neighborEnd_;
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
        std::vector<FieldVector<double, dimgrid> > coordinates(dimgrid);

        // Get reference element
        const Dune::ReferenceElement<double, dimgrid>& refElement
            = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        for (int idx = 0; idx < dimgrid; ++idx)
            coordinates[idx] = refElement.position(refElement.subEntity(this->facetIndex_, 1, idx, dimgrid), dimgrid);

        geometryInInside_ = std::make_shared<LocalGeometryImpl>(this->type(), coordinates);

      return LocalGeometry(*geometryInInside_);
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    //! In the LevelIntersection we know that the intersection is conforming
    LocalGeometry geometryInOutside (std::size_t neighborIndex = 0) const {

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dimgrid>& refElement
           = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        std::array<FoamGridEntityImp<0, dimgrid, dimworld>*, dimgrid> vtx;

        for (int idx = 0; idx < dimgrid; ++idx)
            vtx[idx] = this->center_->vertex_[refElement.subEntity(this->facetIndex_, 1, idx, dimgrid)];

        std::vector<FieldVector<double, dimgrid> > coordinates(dimgrid);

        // Find the intersection vertices in local numbering of the outside element
        // That way we get the local orientation correctly.
        const Dune::ReferenceElement<double, dimgrid>& refElementOther
            = Dune::ReferenceElements<double, dimgrid>::general((*this->neighbor_)->type());

        for (std::size_t j=0; j<dimgrid; j++)
          for (int i=0; i<refElementOther.size(dimgrid); i++)
             if (vtx[j] == (*this->neighbor_)->vertex_[refElementOther.subEntity(0, 0, i, dimgrid)])
              coordinates[j] = refElement.position(refElement.subEntity(0, 0, i, dimgrid), dimgrid);

        geometryInOutside_ = std::make_shared<LocalGeometryImpl>(this->type(), coordinates);

        return LocalGeometry(*geometryInOutside_);
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const {

        std::vector<FieldVector<double, dimworld> > coordinates(dimgrid);

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double, dimgrid>& refElement
            = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        for (std::size_t idx = 0; idx < dimgrid; ++idx)
            coordinates[idx] = this->center_->vertex_[refElement.subEntity(this->facetIndex_, 1, idx, dimgrid)]->pos_;

        geometry_ = std::make_shared<GeometryImpl>(this->type(), coordinates);

        return Geometry(*geometry_);
    }


    private:
    /** \brief One-after-last iterator to the neighbor */
    typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid ,dimworld>*>::const_iterator neighborEnd_;
    //! pointer to global and local intersection geometries
    mutable std::shared_ptr<GeometryImpl> geometry_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInOutside_;

};




/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of an element!
*/
template<class GridImp>
class FoamGridLeafIntersection
    : public FoamGridIntersection<GridImp>
{

    friend class FoamGridLeafIntersectionIterator<GridImp>;

    template<typename, typename>
    friend class Dune::Intersection;

    FoamGridLeafIntersection() : FoamGridIntersection<GridImp>() {}

    enum {dimworld = GridImp::dimensionworld};
    enum {dimgrid  = GridImp::dimension};

    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef typename GridImp::Traits::template Codim<1>::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim<1>::LocalGeometryImpl LocalGeometryImpl;

    FoamGridLeafIntersection(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* center,
                             std::size_t facet)
        : FoamGridIntersection<GridImp>(center, facet)
    {}

public:

    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid leaf? intersections are always? conforming
        return true;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside (std::size_t neighborIndex = 0) const
    {
        // get our facet object
        auto* facet = this->center_->facet_[this->facetIndex_];

        // if the neighbor level is greater than the facet level
        // walk until our facet's father that has the same level as the neighbor
        // and look for the father facet in the neighbor to get the index
        while(facet->level() > (*this->neighbor_)->level())
        {
          assert(facet->father_ != nullptr);
          facet = facet->father_;
        }

        if (facet->level() == (*this->neighbor_)->level())
        {
          return std::distance((*this->neighbor_)->facet_.begin(),
                                std::find((*this->neighbor_)->facet_.begin(),
                                          (*this->neighbor_)->facet_.end(),
                                          facet)
                              );

        }

        // the neighbor level is smaller than the facet level
        // find the neighbor's facet that has the center element in it's element list
        else
        {
          return std::distance((*this->neighbor_)->facet_.begin(),
                                std::find_if((*this->neighbor_)->facet_.begin(),
                                             (*this->neighbor_)->facet_.end(),
                                             [this](const auto& facet) -> bool {
                                               for (const auto& e : facet->elements_)
                                                 if (this->center_ == e)
                                                   return true;
                                               return false;
                                             })
                              );
        }
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
        std::vector<FieldVector<double, dimgrid> > coordinates(dimgrid);

        // Get reference element
        const Dune::ReferenceElement<double, dimgrid>& refElement
            = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        for (std::size_t idx = 0; idx < dimgrid; ++idx)
            coordinates[idx] = refElement.position(refElement.subEntity(this->facetIndex_, 1, idx, dimgrid), dimgrid);

        geometryInInside_ = std::make_shared<LocalGeometryImpl>(this->type(), coordinates);

      return LocalGeometry(*geometryInInside_);
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    //! In the LevelIntersection we know that the intersection is conforming
    LocalGeometry geometryInOutside (std::size_t neighborIndex = 0) const {

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double,dimgrid>& refElement
           = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        std::array<FoamGridEntityImp<0, dimgrid, dimworld>*, dimgrid> vtx;

        for (std::size_t idx = 0; idx < dimgrid; ++idx)
            vtx[idx] = this->center_->vertex_[refElement.subEntity(this->facetIndex_, 1, idx, dimgrid)];

        std::vector<FieldVector<double, dimgrid> > coordinates(dimgrid);

        // Find the intersection vertices in local numbering of the outside element
        // That way we get the local orientation correctly.
        const Dune::ReferenceElement<double, dimgrid>& refElementOther
            = Dune::ReferenceElements<double, dimgrid>::general((*this->neighbor_)->type());

        for (std::size_t j=0; j< dimgrid; j++)
          for (int i=0; i<refElementOther.size(dimgrid); i++)
             if (vtx[j] == (*this->neighbor_)->vertex_[refElementOther.subEntity(0, 0, i, dimgrid)])
              coordinates[j] = refElement.position(refElement.subEntity(0, 0, i, dimgrid), dimgrid);

        geometryInOutside_ = std::make_shared<LocalGeometryImpl>(this->type(), coordinates);

        return LocalGeometry(*geometryInOutside_);
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const {

        std::vector<FieldVector<double, dimworld> > coordinates(dimgrid);

        // Get two vertices of the intersection
        const Dune::ReferenceElement<double, dimgrid>& refElement
            = Dune::ReferenceElements<double, dimgrid>::general(this->center_->type());

        for (std::size_t idx = 0; idx < dimgrid; ++idx)
            coordinates[idx] = this->center_->vertex_[refElement.subEntity(this->facetIndex_, 1, idx, dimgrid)]->pos_;

        geometry_ = std::make_shared<GeometryImpl>(this->type(), coordinates);

        return Geometry(*geometry_);
    }

    //! return true if across the facet a neighbor on this level exists
    bool neighbor () const {
        return this->neighbor_!=FoamGridNullIteratorFactory<dimgrid, dimworld>::null();
    }

private:

    FoamGridEntityImp<dimgrid-1, dimgrid ,dimworld>* facetPointer_;
    //! pointer to global and local intersection geometries
    mutable std::shared_ptr<GeometryImpl> geometry_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInInside_;
    mutable std::shared_ptr<LocalGeometryImpl> geometryInOutside_;
};


}  // namespace Dune

#endif
