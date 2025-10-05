// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:

#ifndef DUNE_FOAMGRID_EDGE_HH
#define DUNE_FOAMGRID_EDGE_HH

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>

namespace Dune {

    /** \brief Edge specialization of FoamGridEntityImp. Edges only exist for dimgrid == 2. */
    template <int dimworld>
    class FoamGridEntityImp<1, 2, dimworld>
        : public FoamGridEntityBase
    {
    public:
        /** \brief Dimension of the grid is always 2 if edges exist. */
        enum {dimgrid = 2};

        FoamGridEntityImp(const FoamGridEntityImp<0, dimgrid, dimworld>* v0,
                          const FoamGridEntityImp<0, dimgrid, dimworld>* v1,
                          int level, unsigned int id)
            : FoamGridEntityBase(level,id), elements_(), nSons_(0), father_(nullptr)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
            sons_[0] =sons_[1] = nullptr;
        }


        FoamGridEntityImp(const FoamGridEntityImp<0, dimgrid, dimworld>* v0,
                          const FoamGridEntityImp<0, dimgrid, dimworld>* v1,
                          int level, unsigned int id,
                          FoamGridEntityImp* father)
            : FoamGridEntityBase(level,id), elements_(), nSons_(0), father_(father)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
            sons_[0] =sons_[1] = nullptr;
        }

        /** \todo Implement this method! */
        bool isLeaf() const {
            return sons_[0]==nullptr;
        }

        //! This has no function yet in Foamgrid
        unsigned int boundarySegmentIndex() const {
            return boundarySegmentIndex_;
        }

        //! This has no function yet in Foamgrid
        unsigned int boundaryId() const {
            return boundaryId_;
        }

        GeometryType type() const {
            return GeometryType(1);
        }

        bool hasFather() const
        {
            return father_!=nullptr;
        }

        /** \brief Number of corners (==2) */
        int corners() const {
            return 2;
        }

        FieldVector<double, dimworld> corner(int i) const {
            return vertex_[i]->pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }

        /** \brief Return level index of sub entity with codim = cc and local number i
         */
        int subLevelIndex (int i, unsigned int codim) const {
            assert(1<=codim && codim<=2);
            switch (codim) {
            case 1:
                return this->levelIndex_;
            case 2:
                return vertex_[i]->levelIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        /** \brief Return leaf index of sub entity with codim = cc and local number i
         */
        int subLeafIndex (int i,unsigned int codim) const {
            assert(1<=codim && codim<=2);
            switch (codim) {
            case 1:
                return this->leafIndex_;
            case 2:
                return vertex_[i]->leafIndex_;
            }
            DUNE_THROW(GridError, "Non-existing codimension requested!");
        }

        std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*> elements_;

        std::array<const FoamGridEntityImp<0, dimgrid, dimworld>*, 2> vertex_;

        /** \brief The boundary id.  Only used if this edge is a boundary edge */
        unsigned int boundarySegmentIndex_;
        unsigned int boundaryId_;

        /** \brief links to refinements of this edge */
        array<FoamGridEntityImp<1, dimgrid, dimworld>*,2> sons_;

        /** \brief The number of refined edges (0 or 2). */
        unsigned int nSons_;

        /** \brief Pointer to father edge */
        FoamGridEntityImp<1, dimgrid, dimworld>* father_;

    };

}

#endif
