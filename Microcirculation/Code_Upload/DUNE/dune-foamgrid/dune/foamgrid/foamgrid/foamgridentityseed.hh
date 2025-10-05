#ifndef DUNE_FOAMGRID_ENTITY_SEED_HH
#define DUNE_FOAMGRID_ENTITY_SEED_HH

/**
 * \file
 * \brief Implementation of EntitySeed for the FoamGrid grid manager
 */

#include "foamgridentity.hh"

namespace Dune {


/**
 * \brief The EntitySeed class provides the minmal information needed to restore an Entity using the grid.
 * \ingroup FoamGrid
 *
 */
template<int codim, class GridImp>
class FoamGridEntitySeed
{
        template<int dimgrid, int dimworld>
        friend class FoamGrid;

    protected:

        enum {dimgrid = GridImp::dimension};
        enum {dimworld = GridImp::dimensionworld};
        enum {mydim = dimgrid-codim};

        // Entity type of the hostgrid
        typedef FoamGridEntityImp<mydim, dimgrid, dimworld> EntityImplType;

    public:

        enum {codimension = codim};

        //! default construct an invalid entity seed
        FoamGridEntitySeed() :
            target_(nullptr)
        {}

        //! construct entity seed from entity
        FoamGridEntitySeed(const FoamGridEntity<codim, dimgrid, GridImp>& entity) :
            target_(entity.target_)
        {}

        FoamGridEntitySeed(const FoamGridEntity<codim, dimgrid, GridImp>* target) :
            target_(target)
        {}

        /** \brief check whether it is safe to create an Entity from this Seed */
        bool isValid() const
        {
          return target_ != nullptr;
        }

    protected:

        const EntityImplType* getImplementationPointer() const
        {
            return target_;
        }

    private:

        const EntityImplType* target_;
};

} // namespace Dune


#endif
