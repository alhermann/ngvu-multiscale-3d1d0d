// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_FACTORY_HH
#define DUNE_FOAMGRID_FACTORY_HH

/** \file
    \brief The specialization of the generic GridFactory for FoamGrid
    \author Oliver Sander
 */

#include <vector>
#include <map>
#include <memory>

#include <dune/common/function.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

namespace Dune {

/** \brief Specialization of the generic GridFactory for FoamGrid<dimgrid, dimworld>
    */
template <int dimgrid, int dimworld>
    class GridFactoryBase
        : public GridFactoryInterface<FoamGrid<dimgrid, dimworld> >
    {
    /** \brief Type used by the grid for coordinates */
    typedef typename FoamGrid<dimgrid, dimworld>::ctype ctype;
    /** \brief Vertex iterator */
    typedef typename std::map<FieldVector<ctype,1>, unsigned int>::iterator VertexIterator;

    public:

        /** \brief Default constructor */
        GridFactoryBase()
            : factoryOwnsGrid_(true),
              vertexIndex_(0)
        {
            grid_ = new FoamGrid<dimgrid, dimworld>;
            grid_->entityImps_.resize(1);
        }

        /** \brief Constructor for a given grid object

        If you already have your grid object constructed you can
        hand it over using this constructor.

        If you construct your factory class using this constructor
        the pointer handed over to you by the method createGrid() is
        the one you supplied here.
         */
        GridFactoryBase(FoamGrid<dimgrid, dimworld>* grid)
            : factoryOwnsGrid_(false),
              vertexIndex_(0)
        {
            grid_ = grid;
            grid_->entityImps_.resize(1);
        }

        /** \brief Destructor */
        virtual ~GridFactoryBase() {
            if (grid_ && factoryOwnsGrid_)
                delete grid_;
        }

        /** \brief Insert a vertex into the coarse grid */
        virtual void insertVertex(const FieldVector<ctype,dimworld>& pos) {
            std::get<0>(grid_->entityImps_[0]).push_back(FoamGridEntityImp<0, dimgrid, dimworld> (0,   // level
                                                                         pos,  // position
                                                                         grid_->freeIdCounter_[0]++));
            vertexArray_.push_back(&*std::get<0>(grid_->entityImps_[0]).rbegin());
        }

        /** \brief Insert a boundary segment.

        This is only needed if you want to control the numbering of the boundary segments
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices) {
            DUNE_THROW(Dune::NotImplemented, "insertBoundarySegment not implemented yet!");
        }

        /** \brief Insert a boundary segment (== a line) and the boundary segment geometry
         *
            This influences the ordering of the boundary segments.
            Currently, the BoundarySegment object does not actually have any effect.
        */
        virtual void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                                           const std::shared_ptr<BoundarySegment<dimgrid, dimworld> >& boundarySegment)
        {
            insertBoundarySegment(vertices);
        }

        /** \brief Return the number of the element to insert parameters
         */
        virtual unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld>::Traits::template Codim<0>::Entity &entity) const
        {
            return grid_->getRealImplementation(entity).target_->leafIndex_;
        }

        /** \brief Return the number of the vertex to insert parameters
         */
        virtual unsigned int
        insertionIndex(const typename FoamGrid<dimgrid, dimworld>::Traits::template Codim<dimgrid>::Entity &vertex) const
        {
            return grid_->getRealImplementation(vertex).target_->leafIndex_;

        }

    protected:

        // Pointer to the grid being built
        FoamGrid<dimgrid, dimworld>* grid_;

        // True if the factory allocated the grid itself, false if the
        // grid was handed over from the outside
        bool factoryOwnsGrid_;

        /** \brief Counter that creates the vertex indices */
        unsigned int vertexIndex_;

        /** \brief Array containing all vertices */
        std::vector<FoamGridEntityImp<0, dimgrid, dimworld>*> vertexArray_;
    };

template <int dimgrid, int dimworld>
    class GridFactory<FoamGrid<dimgrid, dimworld> >
        : public GridFactoryBase<dimgrid, dimworld>
    {};

/** \brief Specialization of GridFactoryBase for 1D-FoamGrid<1, dimworld>
    */
template <int dimworld>
    class GridFactory<FoamGrid<1, dimworld> >
        : public GridFactoryBase<1, dimworld>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 1};
        typedef typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator FacetIterator;
        typedef typename FoamGrid<1, dimworld>::ctype ctype;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<1, dimworld>* grid):
            GridFactoryBase<1,dimworld>(grid)
        {}

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices) {
            assert(type.isLine());
            FoamGridEntityImp<1, dimgrid, dimworld> newElement(this->vertexArray_[vertices[0]],
                                                               this->vertexArray_[vertices[1]],
                                                               0,
                                                               this->grid_->freeIdCounter_[1]++);

            std::get<1>(this->grid_->entityImps_[0]).push_back(newElement);

        }

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization)
        {
            assert(type.isLine());
            FoamGridEntityImp<1, dimgrid, dimworld> newElement(this->vertexArray_[vertices[0]],
                                                               this->vertexArray_[vertices[1]],
                                                               0,
                                                               this->grid_->freeIdCounter_[1]++);
            // save the pointer to the element parametrization
            newElement.elementParametrization_ = elementParametrization;

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid

        The receiver takes responsibility of the memory allocated for the grid
        */
        virtual FoamGrid<1, dimworld>* createGrid() {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            typename std::list<FoamGridEntityImp<1, dimgrid, dimworld> >::iterator eIt    = std::get<1>(this->grid_->entityImps_[0]).begin();
            typename std::list<FoamGridEntityImp<1, dimgrid, dimworld> >::iterator eEndIt = std::get<1>(this->grid_->entityImps_[0]).end();

            for(;eIt!=eEndIt;eIt++) {

                // Get two vertices of the edge
                const FoamGridEntityImp<0, dimgrid, dimworld>* v0 = eIt->vertex_[0];
                const FoamGridEntityImp<0, dimgrid, dimworld>* v1 = eIt->vertex_[1];

                // make vertices know about edge
                // using const_cast because of the implementation of FoamGridEntityImp<1,dimgrid,dimworld>
                // the member variable vertex_ is an array with pointers to const vertices
                const_cast <FoamGridEntityImp<0, dimgrid, dimworld>*> (v0)->elements_.push_back(&*eIt);
                const_cast <FoamGridEntityImp<0, dimgrid, dimworld>*> (v1)->elements_.push_back(&*eIt);

            }

            // Create the index sets
            this->grid_->setIndices();

            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            unsigned int boundaryFacetCounter = 0;

            // Iterate over all facets (=vertices in 1d)
            FacetIterator fIt = std::get<0>(this->grid_->entityImps_[0]).begin();
            const FacetIterator fEndIt = std::get<0>(this->grid_->entityImps_[0]).end();
            for (; fIt != fEndIt; ++fIt)
                if(fIt->elements_.size()==1) // if boundary facet
                    fIt->boundarySegmentIndex_ = boundaryFacetCounter++;

            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld>* tmp = this->grid_;
            tmp->numBoundarySegments_ = boundaryFacetCounter;
            this->grid_ = nullptr;
            return tmp;
        }

    };

    /** \brief Specialization of GridFactoryBase for 2D-FoamGrid<2, dimworld>
    */
template <int dimworld>
    class GridFactory<FoamGrid<2, dimworld> >
        : public GridFactoryBase<2, dimworld>
    {
        /** \brief Grid dimension */
        enum {dimgrid = 2};
        typedef typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator FacetIterator;
        typedef typename FoamGrid<2, dimworld>::ctype ctype;

    public:

        GridFactory() {}

        GridFactory(FoamGrid<2, dimworld>* grid):
            GridFactoryBase<2,dimworld>(grid)
        {}

        /** \brief Insert an element into the coarse grid
            \param type The GeometryType of the new element
            \param vertices The vertices of the new element, using the DUNE numbering
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices) {

            assert(type.isTriangle());

            FoamGridEntityImp<dimgrid, dimgrid, dimworld> newElement(0,   // level
                                       this->grid_->freeIdCounter_[dimgrid]++);  // id
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Insert a parametrized element into the coarse grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        */
        virtual void insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization)
        {
            assert(type.isTriangle());
            FoamGridEntityImp<dimgrid, dimgrid, dimworld> newElement(0,   // level
                                       this->grid_->freeIdCounter_[dimgrid]++);  // id
            newElement.vertex_[0] = this->vertexArray_[vertices[0]];
            newElement.vertex_[1] = this->vertexArray_[vertices[1]];
            newElement.vertex_[2] = this->vertexArray_[vertices[2]];
            // save the pointer to the element parametrization
            newElement.elementParametrization_ = elementParametrization;

            std::get<dimgrid>(this->grid_->entityImps_[0]).push_back(newElement);
        }

        /** \brief Finalize grid creation and hand over the grid
        The receiver takes responsibility of the memory allocated for the grid
        */
        virtual FoamGrid<dimgrid, dimworld>* createGrid() {
            // Prevent a crash when this method is called twice in a row
            // You never know who may do this...
            if (this->grid_==nullptr)
                return nullptr;

            // ////////////////////////////////////////////////////
            //   Create the edges
            // ////////////////////////////////////////////////////

            // for convenience
            typedef FoamGridEntityImp<0, dimgrid, dimworld> FoamGridVertex;

            // For fast retrieval: a map from pairs of vertices to the edge that connects them
            std::map<std::pair<const FoamGridEntityImp<0, dimgrid, dimworld>*, const FoamGridEntityImp<0, dimgrid, dimworld>*>, FoamGridEntityImp<1, dimgrid, dimworld>*> edgeMap;

            typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator eIt    = std::get<dimgrid>(this->grid_->entityImps_[0]).begin();
            typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator eEndIt = std::get<dimgrid>(this->grid_->entityImps_[0]).end();

            for (; eIt!=eEndIt; ++eIt) {

                FoamGridEntityImp<dimgrid, dimgrid, dimworld>* element = &(*eIt);

                const Dune::ReferenceElement<double,dimgrid>& refElement
                    = Dune::ReferenceElements<double, dimgrid>::general(eIt->type());

                // Loop over all edges of this element
                for (size_t i=0; i<element->facet_.size(); ++i) {

                    // Get two vertices of the potential edge
                    const FoamGridVertex* v0 = element->vertex_[refElement.subEntity(i, 1, 0, 2)];
                    const FoamGridVertex* v1 = element->vertex_[refElement.subEntity(i, 1, 1, 2)];

                    FoamGridEntityImp<1, dimgrid, dimworld>* existingEdge = nullptr;
                    typename std::map<std::pair<const FoamGridEntityImp<0, dimgrid, dimworld>*, const FoamGridEntityImp<0, dimgrid, dimworld>*>, FoamGridEntityImp<1, dimgrid, dimworld>*>::const_iterator e = edgeMap.find(std::make_pair(v0,v1));

                    if (e != edgeMap.end()) {
                        existingEdge = e->second;
                    } else {
                        e = edgeMap.find(std::make_pair(v1,v0));
                        if (e != edgeMap.end())
                            existingEdge = e->second;
                    }

                    if (existingEdge == nullptr) {

                        // The current edge has not been inserted already.  We do that now.
                        std::get<1>(this->grid_->entityImps_[0]).push_back(FoamGridEntityImp<1, dimgrid, dimworld>(v0,
                                                                                                    v1,
                                                                                                    0,   // level
                                                                                                    this->grid_->freeIdCounter_[1]++    // id
                                                                                                    ));

                        existingEdge = &*std::get<1>(this->grid_->entityImps_[0]).rbegin();

                        edgeMap.insert(std::make_pair(std::make_pair(v0,v1), existingEdge));

                    }

                    // make element know about the edge
                    element->facet_[i] = existingEdge;

                    // make edge know about the element
                    existingEdge->elements_.push_back(element);
                }
            }

            // Create the index sets
            this->grid_->setIndices();


            // ////////////////////////////////////////////////
            //   Set the boundary ids
            // ////////////////////////////////////////////////

            unsigned int boundaryFacetCounter = 0;

            // Iterate over all facets (=edges in 2D)
            FacetIterator fIt = std::get<1>(this->grid_->entityImps_[0]).begin();
            const FacetIterator fEndIt = std::get<1>(this->grid_->entityImps_[0]).end();
            for (; fIt!=fEndIt; ++fIt)
                if(fIt->elements_.size()==1) //if boundary facet
                    fIt->boundarySegmentIndex_ = boundaryFacetCounter++;


            // ////////////////////////////////////////////////
            //   Hand over the new grid
            // ////////////////////////////////////////////////

            Dune::FoamGrid<dimgrid, dimworld>* tmp = this->grid_;
            tmp->numBoundarySegments_ = boundaryFacetCounter;
            this->grid_ = nullptr;
            return tmp;
        }
    };
}

#endif
