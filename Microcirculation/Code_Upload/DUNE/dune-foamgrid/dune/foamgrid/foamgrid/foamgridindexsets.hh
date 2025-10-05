// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INDEXSETS_HH
#define DUNE_FOAMGRID_INDEXSETS_HH

/** \file
* \brief The index and id sets for the FoamGrid class
*/

#include <vector>
#include <list>

#include <dune/common/version.hh>

#include <dune/grid/common/indexidset.hh>

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>  // for FoamGridEntityImp

namespace Dune {

    /** \todo Take the index types from the host grid */
    template<class GridImp>
    class FoamGridLevelIndexSet :
        public IndexSet<GridImp,FoamGridLevelIndexSet<GridImp> >
    {

        /** \brief Dimension of the grid */
        enum {dimgrid  = GridImp::dimension};

        /** \brief Dimension of the space that the grid is embedded in */
        enum {dimworld = GridImp::dimensionworld};

    public:

        FoamGridLevelIndexSet(const GridImp& grid, int level)
        : grid_(&grid), level_(level), numQuads_(0), numTriangles_(0), numEdges_(0), numVertices_(0)
        {}

        //! get index of an entity
        template<int codim>
        int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
        {
            return GridImp::getRealImplementation(e).target_->levelIndex_;
        }

        //! get index of subentity of an entity
        template<int cc>
        int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e,
                      int i,
                      unsigned int codim) const
        {
            return GridImp::getRealImplementation(e).target_->subLevelIndex(i,codim);
        }

        //! get number of entities of given codim, type and on this level
        int size (int codim) const {
            if(dimgrid == 2) {
                switch (codim) {
                case 0:
                    return numTriangles_ + numQuads_;
                case 1:
                    return numEdges_;
                case 2:
                    return numVertices_;
                }
            } else { //dimgrid==1
                switch (codim) {
                case 0:
                    return numEdges_;
                case 1:
                    return numVertices_;
                }
            }

            return 0;
        }


        //! get number of entities of given codim, type and on this level
        int size (GeometryType type) const
        {
            if (type.isVertex())
                return numVertices_;
            if (type.isLine())
                return numEdges_;
            if (type.isTriangle())
                return numTriangles_;
            if (type.isQuadrilateral())
                return numQuads_;
            return 0;
        }


        /** \brief Deliver all geometry types used in this grid */
        const std::vector<GeometryType>& geomTypes (int codim) const
        {
            return myTypes_[codim];
        }

         /** \brief Deliver all geometry types used in this grid */
        std::vector<GeometryType> types (int codim) const
        {
            assert(codim<=dimgrid);
            return myTypes_[codim];
        }

        /** \brief Return true if the given entity is contained in the index set

        This checks only for the level.  We assume that e belongs to the correct grid
        */
        template<class EntityType>
        bool contains (const EntityType& e) const
        {
            return level_ == e.level();
        }

        /** \brief Set up the index set */
        void update()
        {
            numQuads_ = 0;
            numTriangles_ = 0;
            numEdges_ = 0;
            numVertices_ = 0;

            // //////////////////////////////
            //   Init the vertex indices
            // //////////////////////////////
            typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::const_iterator vIt;
            for (vIt =  std::get<0>(grid_->entityImps_[level_]).begin();
                 vIt != std::get<0>(grid_->entityImps_[level_]).end();
                 ++vIt)
                /** \todo Remove this const cast */
                *const_cast<unsigned int*>(&(vIt->levelIndex_)) = numVertices_++;

             // ///////////////////////////////
            //   Init the edges(2d) / element(1d) indices
            // ///////////////////////////////
            typename std::list<FoamGridEntityImp<1, dimgrid, dimworld> >::const_iterator edIt;
            for (edIt =  std::get<1>(grid_->entityImps_[level_]).begin();
                 edIt != std::get<1>(grid_->entityImps_[level_]).end();
                 ++edIt)
             /** \todo Remove this const cast */
                 *const_cast<unsigned int*>(&(edIt->levelIndex_)) = numEdges_++;


            // ///////////////////////////////
            //   Init the element (2d) indices
            // ///////////////////////////////
            if(dimgrid == 2) {

                typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::const_iterator eIt;
                for (eIt =  std::get<dimgrid>(grid_->entityImps_[level_]).begin();
                    eIt != std::get<dimgrid>(grid_->entityImps_[level_]).end();
                    ++eIt)
                    /** \todo Remove this const cast */
                    *const_cast<unsigned int*>(&(eIt->levelIndex_)) = (eIt->type().isTriangle()) ? numTriangles_++ : numQuads_++;
            }

            // ///////////////////////////////////////////////
            //   Update the list of geometry types present
            // ///////////////////////////////////////////////

            for (int i = 0; i < dimgrid+1; ++i)
                myTypes_[i].clear();

            if (numTriangles_>0 && dimgrid == 2)
                myTypes_[0].push_back(GeometryType(GeometryType::simplex, 2));

            if (numQuads_>0 && dimgrid == 2)
                myTypes_[0].push_back(GeometryType(GeometryType::cube, 2));

            if (numEdges_>0 && dimgrid == 2)
                myTypes_[1].push_back(GeometryType(GeometryType::simplex, 1));
            if (numEdges_>0 && dimgrid == 1)
                myTypes_[0].push_back(GeometryType(GeometryType::simplex, 1));

            if (numVertices_>0)
                myTypes_[dimgrid].push_back(GeometryType(GeometryType::simplex, 0));
        }

    private:
        const GridImp* grid_;
        int level_;

        int numQuads_;
        int numTriangles_;
        int numEdges_;
        int numVertices_;

        /** \brief The GeometryTypes present for each codim */
        array<std::vector<GeometryType>, dimgrid+1> myTypes_;

    };


template<class GridImp>
class FoamGridLeafIndexSet :
    public IndexSet<GridImp,FoamGridLeafIndexSet<GridImp> >
{

    // Grid dimension
    enum {dimgrid  = std::remove_const<GridImp>::type::dimension};
    // World dimension
    enum {dimworld = std::remove_const<GridImp>::type::dimensionworld};

public:

    /** \brief Copy constructor */
    FoamGridLeafIndexSet(const GridImp& g): grid_(g)
    {}

        //! get index of an entity
        /*
            We use the RemoveConst to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
        template<int codim>
        int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
        {
            return GridImp::getRealImplementation(e).target_->leafIndex_;
        }

        //! get index of subentity of an entity
        template<int cc>
        int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e,
                      int i,
                      unsigned int codim) const
        {
            return GridImp::getRealImplementation(e).target_->subLeafIndex(i,codim);
        }

        //! get number of entities of given codim, type and on this level
        int size (int codim) const {
            if(dimgrid == 2) {
                switch (codim) {
                case 0:
                    return numTriangles_ + numQuads_;
                case 1:
                    return numEdges_;
                case 2:
                    return numVertices_;
                }
            } else { //if dimgrid==1
                switch (codim) {
                case 0:
                    return numEdges_;
                case 1:
                    return numVertices_;
                }
            }

            return 0;
        }


        //! get number of entities of given codim, type and on this level
        int size (GeometryType type) const
        {
            if (type.isVertex())
                return numVertices_;
            if (type.isLine())
                return numEdges_;
            if (type.isTriangle())
                return numTriangles_;
            if (type.isQuadrilateral())
                return numQuads_;
            return 0;
        }


        /** \brief Deliver all geometry types used in this grid */
        const std::vector<GeometryType>& geomTypes (int codim) const
        {
            return myTypes_[codim];
        }

        /** \brief Deliver all geometry types used in this grid */
        std::vector<GeometryType> types (int codim) const
        {
            assert(codim<=dimgrid);
            return myTypes_[codim];
        }

        /** \brief Return true if the given entity is contained in the index set */
        template<class EntityType>
        bool contains (const EntityType& e) const
        {
            return GridImp::getRealImplementation(e).target_->isLeaf();
        }

        /** Recompute the leaf numbering */
        void update()
        {

        numQuads_ = 0;
        numTriangles_ = 0;
        numEdges_ = 0;
        numVertices_ = 0;

        // ///////////////////////////////
        //  Init codim 0 entity indices
        // ///////////////////////////////
        typename GridImp::Traits::template Codim<0>::LeafIterator eIt    = grid_.template leafbegin<0>();
        typename GridImp::Traits::template Codim<0>::LeafIterator eEndIt = grid_.template leafend<0>();

        for (; eIt!=eEndIt; ++eIt){
            if(eIt->type().isTriangle() && dimgrid == 2)
                *const_cast<unsigned int*>(&(GridImp::getRealImplementation(*eIt).target_->leafIndex_)) = numTriangles_++;
            if(eIt->type().isQuadrilateral() && dimgrid == 2)
                *const_cast<unsigned int*>(&(GridImp::getRealImplementation(*eIt).target_->leafIndex_)) = numQuads_++;
            if(eIt->type().isLine() && dimgrid == 1)
                *const_cast<unsigned int*>(&(GridImp::getRealImplementation(*eIt).target_->leafIndex_)) = numEdges_++;
        }

        // //////////////////////////////
        //   Init the codim 1 entitiy indices
        // //////////////////////////////

        for (int i=grid_.maxLevel(); i>=0; i--) {

            typename GridImp::Traits::template Codim<1>::LevelIterator edIt    = grid_.template lbegin<1>(i);
            typename GridImp::Traits::template Codim<1>::LevelIterator edEndIt = grid_.template lend<1>(i);

            for (; edIt!=edEndIt; ++edIt) {

            const FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* target = GridImp::getRealImplementation(*edIt).target_;

            if (target->isLeaf()){
                // The is a real leaf edge.
                if(edIt->type().isLine() && dimgrid==2)
                    *const_cast<unsigned int*>(&(target->leafIndex_)) = numEdges_++;
                if(edIt->type().isVertex() && dimgrid==1)
                    *const_cast<unsigned int*>(&(target->leafIndex_)) = numVertices_++;
            }
            else
            {
                if(target->nSons_==1)
                    // If there is green refinement an edge might only have
                    // one son. In this case son and father are identical and
                    // we inherit the leafIndex from the son.
                    *const_cast<unsigned int*>(&(target->leafIndex_)) = target->sons_[0]->leafIndex_;
                }
            }
        }

        // //////////////////////////////////////////////
        //   Init the codim 2 entity indices (only in 2d)
        // //////////////////////////////////////////////

        if(dimgrid==2){

            for (int i=grid_.maxLevel(); i>=0; i--) {

                typename GridImp::Traits::template Codim<dimgrid>::LevelIterator vIt    = grid_.template lbegin<dimgrid>(i);
                typename GridImp::Traits::template Codim<dimgrid>::LevelIterator vEndIt = grid_.template lend<dimgrid>(i);

                for (; vIt!=vEndIt; ++vIt) {

                    const FoamGridEntityImp<0, dimgrid, dimworld>* target = GridImp::getRealImplementation(*vIt).target_;

                    if (target->isLeaf())
                    {
                        if(vIt->type().isVertex() && dimgrid==2)
                            *const_cast<unsigned int*>(&(target->leafIndex_)) = numVertices_++;
                    }
                    else
                        *const_cast<unsigned int*>(&(target->leafIndex_)) = target->sons_[0]->leafIndex_;

                }
            }
        }

        // ///////////////////////////////////////////////
        //   Update the list of geometry types present
        // ///////////////////////////////////////////////
        for (int i = 0; i < dimgrid+1; ++i)
            myTypes_[i].clear();


        if (numTriangles_>0 && dimgrid == 2)
            myTypes_[0].push_back(GeometryType(GeometryType::simplex, 2));

        if (numQuads_>0 && dimgrid == 2)
            myTypes_[0].push_back(GeometryType(GeometryType::cube, 2));

        if (numEdges_>0 && dimgrid == 2)
            myTypes_[1].push_back(GeometryType(1));
        if (numEdges_>0 && dimgrid == 1)
            myTypes_[0].push_back(GeometryType(1));

        if (numVertices_>0)
            myTypes_[dimgrid].push_back(GeometryType(0));
        }

private:

    const GridImp& grid_;

    int numQuads_;
    int numTriangles_;
    int numEdges_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    array<std::vector<GeometryType>, dimgrid+1> myTypes_;

};




template <class GridImp>
class FoamGridIdSet :
    public IdSet<GridImp,FoamGridIdSet<GridImp>, unsigned int>
{

    public:
        //! define the type used for persistent indices
        typedef unsigned int IdType;


        //! get id of an entity
        /*
        We use the std::remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
        */
        template<int cd>
        IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
        {
            return GridImp::getRealImplementation(e).target_->id_;
        }


        //! get id of subEntity
        /*
            We use the std::remove_const to extract the Type from the mutable class,
            because the const class is not instantiated yet.
        */
        IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
        {
            return GridImp::getRealImplementation(e).subId(i,codim);
        }


        /** \todo Should be private */
        void update() {}

};

}  // namespace Dune


#endif
