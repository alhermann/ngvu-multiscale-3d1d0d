// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_HH
#define DUNE_FOAMGRID_HH

/** \file
* \brief The FoamGrid class
*/

#include <list>
#include <set>
#include <map>
#include <tuple>

#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// Implementation classes
#include "foamgrid/foamgridvertex.hh"
#include "foamgrid/foamgridedge.hh"
#include "foamgrid/foamgridelements.hh"

// The components of the FoamGrid interface
#include "foamgrid/foamgridgeometry.hh"
#include "foamgrid/foamgridentity.hh"
#include "foamgrid/foamgridentitypointer.hh"
#include "foamgrid/foamgridentityseed.hh"
#include "foamgrid/foamgridintersectioniterators.hh"
#include "foamgrid/foamgridleveliterator.hh"
#include "foamgrid/foamgridleafiterator.hh"
#include "foamgrid/foamgridhierarchiciterator.hh"
#include "foamgrid/foamgridindexsets.hh"

namespace Dune {

// Forward declaration
template <int dimgrid, int dimworld>
class FoamGrid;


/** \brief Encapsulates loads of types exported by FoamGrid */
template<int dimgrid, int dimworld>
struct FoamGridFamily
{
    typedef GridTraits<
        dimgrid,   // dimgrid
        dimworld,   // dimworld
        Dune::FoamGrid<dimgrid, dimworld>,
        FoamGridGeometry,
        FoamGridEntity,
        FoamGridLevelIterator,
        FoamGridLeafIntersection,
        FoamGridLevelIntersection,
        FoamGridLeafIntersectionIterator,
        FoamGridLevelIntersectionIterator,
        FoamGridHierarchicIterator,
        FoamGridLeafIterator,
        FoamGridLevelIndexSet< const FoamGrid<dimgrid, dimworld> >,
        FoamGridLeafIndexSet< const FoamGrid<dimgrid, dimworld> >,
        FoamGridIdSet< const FoamGrid<dimgrid, dimworld> >,  // global IdSet
        unsigned int,   // global id type
        FoamGridIdSet< const FoamGrid<dimgrid, dimworld> >,  // local IdSet
        unsigned int,   // local id type
        CollectiveCommunication<Dune::FoamGrid<dimgrid, dimworld> > ,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        FoamGridEntitySeed
            > Traits;
};



/** \brief An implementation of the Dune grid interface: a 1- or 2-dimensional simplicial grid in an n-dimensional world
 *
 * \tparam dimgrid Dimension of the grid; must be either 1 or 2
 * \tparam dimworld Dimension of the world space
 */
template <int dimgrid, int dimworld>
class FoamGrid :
        public GridDefaultImplementation  <dimgrid, dimworld, double, FoamGridFamily<dimgrid, dimworld> >
{

    friend class FoamGridLevelIndexSet<const FoamGrid >;
    friend class FoamGridLeafIndexSet<const FoamGrid >;
    friend class FoamGridIdSet<const FoamGrid >;
    friend class FoamGridHierarchicIterator<const FoamGrid >;
    friend class FoamGridLevelIntersectionIterator<const FoamGrid >;
    friend class FoamGridLeafIntersectionIterator<const FoamGrid >;
    friend class FoamGridLevelIntersection<const FoamGrid >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLeafIterator;

    template <class GridType_>
    friend class GridFactory;

    template <int dimgrid_, int dimworld_>
    friend class GridFactoryBase;

    template<int codim_, int dim_, class GridImp_>
    friend class FoamGridEntity;

    template <int codim_, class GridImp_>
    friend class FoamGridEntityPointer;

    public:

    /** \brief FoamGrid is only implemented for 1 and 2 dimension */
    static_assert(dimgrid==1 || dimgrid==2, "Use FoamGrid only for 1d and 2d grids!");

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    typedef FoamGridFamily<dimgrid, dimworld>  GridFamily;

    //! Exports various types belonging to this grid class
    typedef typename FoamGridFamily<dimgrid, dimworld>::Traits Traits;

    //! The type used to store coordinates
    typedef double ctype;

    /** \brief Constructor, constructs an empty grid
     */
    FoamGrid()
        : leafIndexSet_(*this),
          globalRefined(0),
          numBoundarySegments_(0),
          growing_(false)
    {
        //static_assert(dimgrid == 2, "FoamGrid currently only works for 2D in nD");
        std::fill(freeIdCounter_.begin(), freeIdCounter_.end(), 0);
    }

        //! Destructor
        ~FoamGrid()
        {
            // Delete level index sets
            for (size_t i=0; i<levelIndexSets_.size(); i++)
                if (levelIndexSets_[i])
                    delete (levelIndexSets_[i]);
        }


        //! Return maximum level defined in this grid. Levels are numbered
        //! 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const {
            return entityImps_.size()-1;
        }


        //! Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimgrid, dimworld> >(std::get<dimgrid-codim>(entityImps_[level]).begin());
        }


        //! one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimgrid, dimworld> >(std::get<dimgrid-codim>(entityImps_[level]).end());
        }


        //! Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimgrid, dimworld> >(std::get<dimgrid-codim>(entityImps_[level]).begin());
        }


        //! one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimgrid, dimworld> >(std::get<dimgrid-codim>(entityImps_[level]).end());
        }


        //! Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >(*this);
        }


        //! one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >();
        }


        //! Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >(*this);
        }


        //! one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >();
        }


        /** \brief Number of grid entities per level and codim
         */
        int size (int level, int codim) const {

            // Turn dynamic index into static index
            if ((codim==2 && dimgrid==2) || (codim==1 && dimgrid==1))
                return std::get<0>(entityImps_[level]).size();
            if ((codim==1 && dimgrid==2))
                return std::get<1>(entityImps_[level]).size();
            if (codim==0)
                return std::get<dimgrid>(entityImps_[level]).size();

            return 0;
        }


        //! number of leaf entities per codim in this process
        int size (int codim) const{
            return leafIndexSet().size(codim);
        }


        //! number of entities per level, codim and geometry type in this process
        int size (int level, GeometryType type) const {
            return this->levelIndexSet(level).size(type);
        }


        //! number of leaf entities per codim and geometry type in this process
        int size (GeometryType type) const
        {
            return this->leafIndexSet().size(type);
        }

        /** \brief The number of boundary edges on the coarsest level */
        size_t numBoundarySegments() const
        {
            return numBoundarySegments_;
        }

        /** \brief Access to the GlobalIdSet */
        const typename Traits::GlobalIdSet& globalIdSet() const{
            return idSet_;
        }


        /** \brief Access to the LocalIdSet */
        const typename Traits::LocalIdSet& localIdSet() const{
            return idSet_;
        }


        /** \brief Access to the LevelIndexSets */
        const typename Traits::LevelIndexSet& levelIndexSet(int level) const
        {
          if (! levelIndexSets_[level])
          {
            levelIndexSets_[level] = new FoamGridLevelIndexSet<const FoamGrid>(*this, level);
            levelIndexSets_[level]->update();
          }
          return *levelIndexSets_[level];
        }


        /** \brief Access to the LeafIndexSet */
        const typename Traits::LeafIndexSet& leafIndexSet() const
        {
            return leafIndexSet_;
        }

        /** \brief Create an Entity from an EntitySeed */
        template <class EntitySeed>
        static typename Traits::template Codim<EntitySeed::codimension>::Entity
        entity(const EntitySeed& seed)
        {
          const int codim = EntitySeed::codimension;
          typedef typename Traits::template Codim<codim>::Entity Entity;
          return Entity(FoamGridEntity<codim, dimgrid, const FoamGrid>(FoamGrid::getRealImplementation(seed).getImplementationPointer()));
        }


        /** @name Grid Refinement Methods */
        /*@{*/


        /** \brief Refine the grid uniformly
         * \param refCount Number of times the grid is to be refined uniformly
        */
        void globalRefine (int refCount = 1);

        /** \brief Mark entity for refinement
        *
        * This only works for entities of codim 0.
        * The parameter is currently ignored
        *
        * \return <ul>
        * <li> true, if marking was successful </li>
        * <li> false, if marking was not possible </li>
        * </ul>
        */
        bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e)
        {
            if (!e.isLeaf())
                return false;

            /** \todo Why do I need those const_casts here? */
            if (refCount>=1)
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE;
            else if (refCount<0)
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN;
            else
                const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(e).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;

            return true;
        }

        /** \brief Return refinement mark for entity
        *
        * \return refinement mark (1,0,-1)
        */
        int getMark(const typename Traits::template Codim<0>::Entity & e) const
        {
            switch(this->getRealImplementation(e).target_->markState_)
            {
              case FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING:
              case FoamGridEntityImp<dimgrid, dimgrid, dimworld>::IS_COARSENED:
                return 0;
              case FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE:
                return 1;
              case FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN:
                return -1;
            }
            return 0;
        }

        //! \brief Book-keeping routine to be called before adaptation
        bool preAdapt();

        //! Triggers the grid refinement process
        bool adapt();

        /** \brief Clean up refinement markers */
        void postAdapt();

        /** \brief Sets a (leaf) vertex to a new position
         *
         *  \param e The codim dimgrid entity (vertex) to be moved. Note: The vertex
         *           must be a leaf vertex. The implementation stores copies of vertices
         *           for each level they exist on. Changing a vertex' position changes
         *           its position on all coarser grid levels, too! We could not think of
         *           an application for moving non-leaf vertices, write us if you need
         *           that feature.
         *  \param pos The new global position of the vertex
         */
        void setPosition(const typename Traits::template Codim<dimgrid>::Entity & e,
                         const FieldVector<ctype, dimworld>& pos);

        /** @name Grid Growth Methods */
        /*@{*/

        /** \brief Add new vertex to be added the grid
        * \param pos The position vector of the vertex
        * \return The index of the newly inserted vertex (to be able to insert elements with it)
        */
        unsigned int insertVertex(const FieldVector<ctype,dimworld>& pos)
        {
          if(!growing_) initializeGrowth_();

          // the final level of the vertex will be the minimum common vertex level of the new element
          verticesToInsert_.push_back(FoamGridEntityImp<0, dimgrid, dimworld>(0, pos, -verticesToInsert_.size()-1)); // initialize with some invalid id
          FoamGridEntityImp<0, dimgrid, dimworld>& newVertex = verticesToInsert_.back();
          newVertex.isNew_ = true;
          // new vertices are numbered consecutively starting from
          // the highest available index in the leaf index set +1
          return this->leafGridView().size(dimgrid) - 1 + verticesToInsert_.size();
        }

        /** \brief Add a new element to be added to the grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \return The growthInsertionIndex that can be used to attach user data to this element.
                It is valid between until calling postGrow.
        */
        unsigned int insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices)
        {
          // foamgrid only supports simplices until now
          assert(type.isTriangle() || type.isLine());
          assert(vertices.size() == dimgrid + 1);

          // the final level of the element will be the minimum common vertex level
          elementsToInsert_.push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(0, -elementsToInsert_.size()-1)); // initialize with some invalid id
          FoamGridEntityImp<dimgrid, dimgrid, dimworld>& newElement = elementsToInsert_.back();
          assert(vertices.size() == newElement.vertex_.size());

          for(std::size_t i = 0; i < vertices.size(); i++)
          {
            if(int(vertices[i]) >= this->leafGridView().size(dimgrid))
            {
              // initialize with pointer to vertex in verticesToInsert_ vector, later overwrite with actual pointer
              auto vIt = verticesToInsert_.begin();
              std::advance(vIt, vertices[i] - this->leafGridView().size(dimgrid));
              newElement.vertex_[i] = &*vIt;
            }
            else
            {
              // make sure the index to vertex map has been initialized
              if(!growing_) initializeGrowth_();
              // the vertex already exists in the grid, initialize with leaf vertex, later overwrite with lowest level father
              assert(indexToVertexMap_[vertices[i]]->isLeaf());
              newElement.vertex_[i] = indexToVertexMap_[vertices[i]];
            }
          }
          newElement.isNew_ = true;
          newElement.growthInsertionIndex_ = elementsToInsert_.size()-1;
          return elementsToInsert_.size()-1;
        }

        /** \brief Add a new element to be added to the grid
        \param type The GeometryType of the new element
        \param vertices The vertices of the new element, using the DUNE numbering
        \param elementParametrization A function prescribing the shape of this element
        \return The growthInsertionIndex that can be used to attach user data to this element.
                It is valid between until calling postGrow.
        */
        unsigned int insertElement(const GeometryType& type,
                                   const std::vector<unsigned int>& vertices,
                                   const std::shared_ptr<VirtualFunction<FieldVector<ctype,dimgrid>,FieldVector<ctype,dimworld> > >& elementParametrization)
        {
          auto growthInsertionIndex = insertElement(type, vertices);
          // save the pointer to the element parametrization
          elementsToInsert_.back().elementParametrization_ = elementParametrization;
          return growthInsertionIndex;
        }

        /** \brief Mark an element for removal from the grid
        \param e The codim 0 entity to be removed from the grid
        */
        void removeElement(const typename Traits::template Codim<0>::Entity & e)
        {
          // save entity for later, actual removal happens in grow()
          elementsToRemove_.push_back(const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*> (this->getRealImplementation(e).target_));
        }

        //! \brief Book-keeping routine to be called before growth
        bool preGrow();

        //! Triggers the grid growth process
        bool grow();

        /** \brief Clean up isNew markers */
        void postGrow();

        /**
         * \brief The index of insertion if the element was created in the current growth step.
         *         If this is the first element added to the growth queue by calling insertElement the index is 0 and so on.
         *         The index will be valid until postGrow is called.
         * \note   This is useful to attach user data to a created element. The data might only known at element creation time.
         *         As the final element index and id are not known yet when the element is added to the insertion queue, this
         *         index can be used to attach user data that can (after calling postGrow) be attached to the real element index/id.
         */
        unsigned int growthInsertionIndex(const typename Traits::template Codim<0>::Entity & e) const
        {
            int idx = this->getRealImplementation(e).target_->growthInsertionIndex_;
            assert(idx >= 0);
            return static_cast<unsigned int>(idx);
        }

        /*@}*/

        /** @name Methods for parallel computations */
        /*@{*/

        /** \brief Size of the overlap on the leaf level */
        unsigned int overlapSize(int codim) const {
            return 0;
        }


        /** \brief Size of the ghost cell layer on the leaf level */
        unsigned int ghostSize(int codim) const {
            return 0;
        }


        /** \brief Size of the overlap on a given level */
        unsigned int overlapSize(int level, int codim) const {
            return 0;
        }


        /** \brief Size of the ghost cell layer on a given level */
        unsigned int ghostSize(int level, int codim) const {
            return 0;
        }


#if 0
        /** \brief Distributes this grid over the available nodes in a distributed machine
        *
        * \param minlevel The coarsest grid level that gets distributed
        * \param maxlevel does currently get ignored
        */
        void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
            DUNE_THROW(NotImplemented, "FoamGrid::loadBalance()");
        }

        /** \brief The communication interface
        *  @param T: array class holding data associated with the entities
        *  @param P: type used to gather/scatter data in and out of the message buffer
        *  @param codim: communicate entites of given codim
        *  @param if: one of the predifined interface types, throws error if it is not implemented
        *  @param level: communicate for entities on the given level
        *
        *  Implements a generic communication function sending an object of type P for each entity
        *  in the intersection of two processors. P has two methods gather and scatter that implement
        *  the protocol. Therefore P is called the "protocol class".
        */
        template<class T, template<class> class P, int codim>
        void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);

        /*! The new communication interface

        communicate objects for all codims on a given level
        */
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
        {}

        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {}
#endif

        /** dummy collective communication */
        const typename Traits::CollectiveCommunication& comm () const
        {
            return ccobj_;
        }
        /*@}*/


        // **********************************************************
        // End of Interface Methods
        // **********************************************************

    private:

    //! \brief Prepares the grid for growth
    void initializeGrowth_()
    {
      // update the index to vertex map
      indexToVertexMap_.resize(this->leafGridView().size(dimgrid));
      for (const auto& vertex : vertices(this->leafGridView()))
      {
        std::size_t index = leafIndexSet().index(vertex);
        indexToVertexMap_[index] = const_cast<FoamGridEntityImp<0, dimgrid ,dimworld>*>(this->getRealImplementation(vertex).target_);
      }

      // tell the grid it's ready for growth
      growing_ = true;
    }

    //! \brief erases pointers in father elements to vanished entities of the element
    void erasePointersToEntities(std::list<FoamGridEntityImp<dimgrid, dimgrid ,dimworld> >& elements);

    //! \brief Erase Entities from memory that vanished due to coarsening.
    //! \warning This method has to be called first for i=0.
    //! \tparam i The dimension of the entities.
    //! \param  levelEntities The vector with the level entitied
    template<int i>
    void eraseVanishedEntities(std::list<FoamGridEntityImp<i, dimgrid, dimworld> >& levelEntities);

    //! \brief Coarsen an Element
    //! \param element The element to coarsen
    void coarsenSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element);

    //! \brief refine an Element
    //! \param element The element to refine
    //! \param refCount How many times to refine the element
    void refineSimplexElement(FoamGridEntityImp<2, dimgrid, dimworld>& element,
                       int refCount);
    //! Overloaded function for the 1d case
    void refineSimplexElement(FoamGridEntityImp<1, dimgrid, dimworld>& element,
                       int refCount);

    //! \brief remove this element resulting in grid shrinkage
    bool removeSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element);

    /**
     * \brief Overwrites the neighbours of this and descendant edges
     *
     * After returning all neighbours the previously pointed to the
     * father will point to the son element
     * \param edge The edge to start overwriting with.
     * \param son The son element to substitute the father with.
     * \param father Pointer to the father element that is to be substituted.
     */
    void overwriteFineLevelNeighbours(FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>& edge,
                                      const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* son,
                                      const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father);

    /**
     * \brief An an element to a facet (and all it's sons if it's not on the leaf)
     *
     * \param element The element to add to the facet's element vector
     * \param facet The facet that needs to add the element
     */
    void addElementForFacet(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* element,
                            FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* facet);

    //! Add a new facet for 1d grids (the facet already exists as vertex)
    void addNewFacet(FoamGridEntityImp<0, dimgrid, dimworld>* &facet,
                     std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray,
                     int level);

    //! Add a new facet for 2d grids
    void addNewFacet(FoamGridEntityImp<1, dimgrid, dimworld>* &facet,
                     std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray,
                     int level);

    //! compute the grid indices and ids
    void setIndices();

    //! Compute the codim 2-0 connectivity useful for removal of elements
    void computeTwoZeroConnectivity()
    {
      // the elements vector for vertices is already available as they are facets in 1d grids
      for(int level = 0; level <= maxLevel(); ++level)
      {
        // do it everytime freshly so that we don't have to take care of 2-0 connectivity in adaptivity
        for(auto&& vertex : std::get<0>(entityImps_[level]))
          vertex.elements_.clear();

        for(auto eIt = std::get<dimgrid>(entityImps_[level]).begin(); eIt != std::get<dimgrid>(entityImps_[level]).end(); ++eIt)
          for(auto&& vertex : eIt->vertex_)
            vertex->elements_.push_back(&*eIt);
      }
    }

    //! Collective communication interface
    typename Traits::CollectiveCommunication ccobj_;

    // Stores the lists of vertices, edges, elements for each level
    // std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
    // std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >,
    // std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> > > > entityImps_;

    // conditional typename depending on dimension of grid (1 or 2)
    typedef typename std::conditional<
                      dimgrid==2,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<2, dimgrid, dimworld> > > >,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> > > >
                   >::type EntityImps;

    EntityImps entityImps_;

    typedef typename std::conditional<
                        dimgrid==2,
                        std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<2, dimgrid, dimworld> > >,
                        std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                     std::list<FoamGridEntityImp<1, dimgrid, dimworld> > >
                   >::type EntityTuple;

    //! Our set of level indices
    mutable std::vector<FoamGridLevelIndexSet<const FoamGrid>*> levelIndexSets_;

    //! The leaf index set
    FoamGridLeafIndexSet<const FoamGrid > leafIndexSet_;

    //! The id set
    FoamGridIdSet<const FoamGrid > idSet_;

    /** \brief Counters that always provide the next free id for each dimension */
    array<unsigned int, dimgrid+1> freeIdCounter_;

    /** \brief How many times was the leaf level globally refined. */
    int globalRefined;

    /** \brief The number of boundary segements. */
    std::size_t numBoundarySegments_;

    // True if the last call to preadapt returned true
    bool willCoarsen;

    /** \brief A map from indices to leaf vertices. Gets updated when calling beginGrowth(). */
    std::vector<FoamGridEntityImp<0, dimgrid, dimworld>* > indexToVertexMap_;

    /** \brief The (temporary) vector of element(pointer)s to be deleted at runtime. Gets cleaned when calling grow(). */
    std::vector<FoamGridEntityImp<dimgrid, dimgrid, dimworld>* > elementsToRemove_;

    /** \brief The (temporary) vector of runtime inserted vertices. Gets cleaned when calling grow(). */
    std::list<FoamGridEntityImp<0, dimgrid, dimworld> > verticesToInsert_;

    /** \brief The (temporary) vector of runtime inserted elements. Gets cleaned when calling grow(). */
    std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> > elementsToInsert_;

    /** \brief If the grid is in a growing process (beginGrowth has been called). */
    bool growing_;

}; // end Class FoamGrid

#include <dune/foamgrid/foamgrid/foamgrid.cc>


namespace Capabilities
{
    /** \brief True if the grid implements entities of a given codim.
      *
      * FoamGrid implements all codimensions, hence this is always true
      */
    template<int dimgrid, int dimworld, int codim>
    struct hasEntity< FoamGrid<dimgrid, dimworld>, codim>
    {
        static const bool v = true;
    };

    //! \todo Please doc me !
    template <int dimgrid, int dimworld>
    struct isLevelwiseConforming< FoamGrid<dimgrid, dimworld> >
    {
        static const bool v = false;
    };

    //! \todo Please doc me !
    template <int dimgrid, int dimworld>
    struct isLeafwiseConforming< FoamGrid<dimgrid, dimworld> >
    {
        static const bool v = false;
    };
}

} // namespace Dune


// The factory should be present whenever the user includes foamgrid.hh.
// However since the factory needs to know the grid the include directive
// comes here at the end.
#include <dune/foamgrid/foamgrid/foamgridfactory.hh>

#endif
