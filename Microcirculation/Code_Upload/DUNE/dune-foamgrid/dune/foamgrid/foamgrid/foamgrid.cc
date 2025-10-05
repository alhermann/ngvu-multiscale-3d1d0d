// Refine the grid uniformly
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::globalRefine (int refCount)
{
  willCoarsen=false;

  if (maxLevel()+refCount<0)
    DUNE_THROW(GridError, "Grid has only " << maxLevel() << " levels. Cannot do "
                           << " globalRefine(" << refCount << ")");

  // The leafiterator is simply successively visiting all levels from fine to coarse
  // and just checking whether the isLeaf flag is set. Using it to identify
  // elements that need refinement would produce an endless loop, as the newly
  // added elements
  // would later be identified as elements that still need refinement.
  std::size_t oldLevels =entityImps_.size();
  //conditional typename depending on dimension of grid (1 or 2)
  typedef typename std::conditional<
                      dimgrid==2,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<2, dimgrid, dimworld> >
                                          > >::reverse_iterator,
                      typename std::vector<std::tuple<std::list<FoamGridEntityImp<0, dimgrid, dimworld> >,
                                                 std::list<FoamGridEntityImp<1, dimgrid, dimworld> >
                                          > >::reverse_iterator
                   >::type LevelIterator;

  // Allocate space for the new levels. Thus the rend iterator will not get
  // invalid due to newly added levels.
  entityImps_.reserve(oldLevels+refCount);
  levelIndexSets_.reserve(oldLevels+refCount);
  LevelIterator level = entityImps_.rbegin();

  // Add tuples for the new levels.
  for (int i=0; i < refCount; ++i)
  {
    entityImps_.push_back(EntityTuple());
    // add space for new LevelIndexSets. They are not created until requested
    levelIndexSets_.push_back( (FoamGridLevelIndexSet<const FoamGrid > *) 0 );
  }

  if (refCount < 0)
  {
    if (globalRefined+refCount<0)
    {
      for (int i=refCount; i<0; ++i)
      {
        // Mark each leaf element for coarsening.
        for (const auto& element : elements(this->leafGridView()))
        {
          mark(-1, element);
        }

        // do the adaptation
        preAdapt();
        adapt();
        postAdapt();
      }

      globalRefined=0;
      return;
    }
    else
    {
      for (int i=refCount; i<0; ++i)
      {
        delete levelIndexSets_.back();
        levelIndexSets_.pop_back();
      }

      entityImps_.resize(oldLevels+refCount);

      // To be able to create the leaf level we need to set
      // the sons of the entities of maxlevel to null
      typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::iterator vIt
        = std::get<0>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<0, dimgrid, dimworld> >::iterator vEndIt
        = std::get<0>(entityImps_[maxLevel()]).end();
      for (; vIt!=vEndIt; ++vIt)
      {
        vIt->sons_[0]=nullptr;
        if(dimgrid == 1)
          vIt->nSons_=0;
      }

      if(dimgrid == 2)
      {
        typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator edIt
          = std::get<dimgrid-1>(entityImps_[maxLevel()]).begin();
        typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator edEndIt
          = std::get<dimgrid-1>(entityImps_[maxLevel()]).end();
        for (; edIt!=edEndIt; ++edIt)
        {
          edIt->sons_[0]=nullptr;
          edIt->sons_[1]=nullptr;
          edIt->nSons_=0;
        }
      }

      typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator elIt
        = std::get<dimgrid>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<dimgrid, dimgrid ,dimworld> >::iterator elEndIt
        = std::get<dimgrid>(entityImps_[maxLevel()]).end();
      for (; elIt!=elEndIt; ++elIt)
      {
        for(int sonIdx=0; sonIdx < (1<<dimgrid); ++sonIdx)
          elIt->sons_[sonIdx]=nullptr;
        elIt->nSons_=0;
      }
    }

  }
  else if (refCount > 0) // for refCount = 0 do nothing
  {
    // sanity check whether the above asssumption really holds.
    assert(&entityImps_[oldLevels-1]==&(*level));

    // Do the actual refinement.
    // We start with the finest level
    std::size_t levelIndex;
    for (levelIndex=oldLevels-1; level!=entityImps_.rend(); ++level, --levelIndex)
    {
      typedef typename std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >::iterator ElementIterator;
      bool foundLeaf=false;

      for (ElementIterator element=std::get<dimgrid>(*level).begin(); element != std::get<dimgrid>(*level).end(); ++element)
        if(element->isLeaf())
        {
          foundLeaf = true;
          if (element->type().isTriangle() || element->type().isLine())
            refineSimplexElement(*element, refCount);
          else
            DUNE_THROW(NotImplemented, "Refinement only supported for simplices!");
        }

      if (!foundLeaf)
        break;
    }
  }

  // Update the leaf indices
  leafIndexSet_.update();

  globalRefined=std::max(globalRefined+refCount,0);
  postAdapt();
}


//f Book-keeping routine to be called before adaptation
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::preAdapt()
{
  // Loop over all leaf entities and check whether they might be
  // coarsened. If there is one return true.
  int addLevels = 0;
  willCoarsen = false;

  for (const auto& element : elements(this->leafGridView()))
  {
    int mark=getMark(element);
    addLevels=std::max(addLevels, element.level()+mark-maxLevel());

    if (mark<0)
    {
      // If the element is marked for coarsening but it is not allowed because of growth
      // i.e. it contains a junction facet with has no father, we reset to DO_NOTHING
      if(this->getRealImplementation(element).target_->coarseningBlocked_)
        const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(element).target_)->markState_ = FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;

      // If this element is marked for coarsening, but another child
      // of this element's father is marked for refinement or has children, then we
      // need to reset the marker to doNothing
      bool otherChildRefined=false;
      FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father = *this->getRealImplementation(element).target_->father_;
      typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::iterator ChildrenIter;
      for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
        otherChildRefined = otherChildRefined ||
                            (*child)->markState_==FoamGridEntityImp<dimgrid, dimgrid, dimworld>::REFINE ||
                            !(*child)->isLeaf();

      if (otherChildRefined)
      {
        for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
          if ((*child)->markState_==FoamGridEntityImp<dimgrid, dimgrid, dimworld>::COARSEN)
            (*child)->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;
      }
      else
        willCoarsen = willCoarsen || mark<0;
    }
  }

  if (addLevels)
  {
    entityImps_.push_back(EntityTuple());
    // add space for new LevelIndexSets. They are not created until requested
    levelIndexSets_.push_back( (FoamGridLevelIndexSet<const FoamGrid > *) 0);
  }

  return willCoarsen;
}


// Triggers the grid refinement process
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::adapt()
{
  bool haveRefined=false;

  // Loop over all leaf elements and refine/coarsen those that marked for it.
  for (const auto& element : elements(this->leafGridView()))
  {
    int mark=getMark(element);
    if (mark>0)
    {
      // Refine simplices
      if (element.type().isTriangle() || element.type().isLine())
      {
        refineSimplexElement(*const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(element).target_), 1);
        haveRefined=true;
      }
      else
        DUNE_THROW(NotImplemented, "Refinement only supported for simplices!");
    }

    if (mark<0) // If simplex was already treated by coarsenSimplex mark is 0
    {
      // Coarsen simplices
      if (element.type().isTriangle() || element.type().isLine())
      {
        assert(element.level());
        coarsenSimplexElement(*const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(element).target_));
      }
      else
        DUNE_THROW(NotImplemented, "Coarsening only supported for simplices!");
    }
  }

  if (!willCoarsen)
  {
    if(haveRefined)
    {
      // Update the leaf indices
      leafIndexSet_.update();
    }
    return haveRefined;
  }

  for (int level = maxLevel(); level >= 0; --level)
  {
    // First delete the pointers to vanishing entities
    erasePointersToEntities(std::get<dimgrid>(entityImps_[level]));

    // Now delete the vertices marked with willVanish_ == true
    if (dimgrid > 1)
      eraseVanishedEntities(std::get<0>(entityImps_[level]));

    // erase vanished facets
    // the erased element were replaced by the father, so we erase all facets that don't have elements
    // on the same level
    typedef typename std::list<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld> >::iterator FacetIter;
    for (FacetIter facet=std::get<dimgrid-1>(entityImps_[level]).begin();
         facet != std::get<dimgrid-1>(entityImps_[level]).end();)
    {
      for (auto&& element : facet->elements_)
      {
        if(element->willVanish_)
        {
          element = element->father_;
        }
      }

      // check if we have same level elements, if so this facet stays
      bool hasSameLevelElements=false;
      for (auto&& element : facet->elements_)
        hasSameLevelElements = hasSameLevelElements || element->level()==level;

      assert(facet->willVanish_!=hasSameLevelElements);

      if (!hasSameLevelElements && facet->willVanish_)
      {
        // erase returns next position
        facet=std::get<dimgrid-1>(entityImps_[level]).erase(facet);
      }
      else
      {
        // increment
        (*facet).willVanish_ = false;
        ++facet;
      }
    }

    // And the elements
    eraseVanishedEntities(std::get<dimgrid>(entityImps_[level]));
  }

  // delete uppermost level if there are no entities in this level anymore
  if(!std::get<0>(entityImps_.back()).size())
  {
    assert(!std::get<dimgrid-1>(entityImps_.back()).size() &&
           !std::get<dimgrid>(entityImps_.back()).size());
    entityImps_.pop_back();
  }

  // Renumber the entities
  setIndices();

  globalRefined=0;

  return haveRefined;
}



// Clean up refinement markers
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::postAdapt()
{
  willCoarsen=false;

  // Loop over all leaf entities and remove the isNew Marker.
  for (int level = maxLevel(); level >= 0; --level)
  {
    for(auto eIt = std::get<dimgrid>(entityImps_[level]).begin(); eIt != std::get<dimgrid>(entityImps_[level]).end(); ++eIt)
    {
      eIt->isNew_=false;
      assert(!eIt->willVanish_);
      if (eIt->father_)
        eIt->father_->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::DO_NOTHING;
    }
  }
}


// Erases pointers in father elements to vanished entities of the element
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::erasePointersToEntities(std::list<FoamGridEntityImp<dimgrid, dimgrid, dimworld> >& elements)
{
  for(auto&& element : elements)
  {
    if(element.willVanish_)
    {
      FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father=*element.father_;
      for (unsigned int i=0; i<father.nSons_; i++)
        father.sons_[i]=nullptr;
      // reset the number of sons for the father
      father.nSons_ = 0;
      for (int i=0; i<father.corners(); i++)
        if (father.vertex_[i]->sons_[0]!=nullptr && father.vertex_[i]->sons_[0]->willVanish_)
        {
          father.vertex_[i]->sons_[0]=nullptr;
          --father.vertex_[i]->nSons_;
        }
      for (int i=0; i<father.corners(); i++)
        if (father.facet_[i]->sons_[0]!=nullptr
            && father.facet_[i]->sons_[0]->willVanish_)
          for (unsigned int j=0; j<dimgrid; j++)
          {
            assert(father.facet_[i]->sons_[j]!=nullptr);
            assert(father.facet_[i]->sons_[j]->willVanish_);
            father.facet_[i]->sons_[j]=nullptr;
            --father.facet_[i]->nSons_;
          }
      assert(father.isLeaf());
    }
  }
}

// Erase Entities from memory that vanished due to coarsening.
template <int dimgrid, int dimworld>
template <int i>
void Dune::FoamGrid<dimgrid, dimworld>::eraseVanishedEntities(std::list<FoamGridEntityImp<i, dimgrid, dimworld> >& levelEntities)
{
  typedef typename std::list<FoamGridEntityImp<i, dimgrid, dimworld> >::iterator EntityIterator;
  for (EntityIterator entity=levelEntities.begin();
      entity != levelEntities.end();)
  {
    if(entity->willVanish_)
      entity=levelEntities.erase(entity);
    else
      ++entity;
  }
}

// Coarsen a simplex element
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::coarsenSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element)
{
  // If we coarsen an element, this means that we erase all chidren of its father
  // to prevent inconsistencies.
  const FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father = *(element.father_);

  // The facets that might be erased
  std::set<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*> childFacets;

  // The vertices that might be erased
  std::set<FoamGridEntityImp<0, dimgrid, dimworld>*> childVertices;

  // Iterate over all children of the father -> elements to be erased
  typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::const_iterator ChildrenIter;
  for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
  {
    // Remember element for the actual deletion taking place later
    (*child)->markState_=FoamGridEntityImp<dimgrid, dimgrid, dimworld>::IS_COARSENED;
    (*child)->willVanish_=true;

    // Iterate over the facets of this vanishing element
    typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;
    for (FacetIter facet=(*child)->facet_.begin(); facet != (*child)->facet_.end(); ++facet)
    {
      // Remove references to elements that will be erased
      typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator
                    ElementIter;
      for (ElementIter element = (*facet)->elements_.begin();
           element != (*facet)->elements_.end();
           ++element)
      {
        for (ChildrenIter child1=father.sons_.begin(); child1 !=
             father.sons_.end(); ++child1)
          if(*element==*child1)
          {
            // To prevent removing an element from a vector
            // we just overwrite the element with its father.
            // later we will check the facets and erase all
            // of them that do not have any elements on the same
            // level.
            *element=&father;
          }
      }
    }
    // Save potential sub entities that might be erased
    if(dimgrid > 1)
    {
      typedef typename array<FoamGridEntityImp<0, dimgrid, dimworld>*, dimgrid+1>::iterator VertexIter;
      for (VertexIter vertex=(*child)->vertex_.begin(); vertex != (*child)->vertex_.end(); ++vertex)
        childVertices.insert(*vertex);
    }
    typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;
    for (FacetIter facet=(*child)->facet_.begin(); facet!= (*child)->facet_.end();
        ++facet)
      childFacets.insert(*facet);
  }

  // Check whether those guys are really erased.
  // We remove all vertices and facets that are part of a child element of one
  // of the neighbours of the father element from the "to be erased" lists
  typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::const_iterator FacetIter;
  for(FacetIter facet=father.facet_.begin(); facet != father.facet_.end(); ++facet)
  {
    typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator NeighborIter;
    for (NeighborIter neighbor = (*facet)->elements_.begin();
         neighbor != (*facet)->elements_.end();
         ++neighbor)
    {
      assert((*neighbor)->level()<=father.level());
      if(*neighbor == &father)
        continue;

      if((*neighbor)->level()==father.level() && !((*neighbor)->isLeaf()))
      {
        // This is a real neighbor element on the same level
        // Check whether one of its children is marked for coarsening
        bool coarsened=false;
        for (ChildrenIter child=(*neighbor)->sons_.begin();
             child != (*neighbor)->sons_.end(); ++child)
          if((*child)->mightVanish())
          {
            coarsened=true;
            break;
          }

        if(!coarsened)
        {
          // Remove all entities that exist in the children of the neighbor
          // from the removal list
          for (ChildrenIter child=(*neighbor)->sons_.begin();
               child != (*neighbor)->sons_.end(); ++child)
          {
            typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIter;
            for (FacetIter facet=(*child)->facet_.begin(); facet != (*child)->facet_.end(); ++facet)
              childFacets.erase(*facet);
            typedef typename array<FoamGridEntityImp<0, dimgrid, dimworld>*, dimgrid+1>::iterator VertexIter;
            for (VertexIter vertex=(*child)->vertex_.begin();
                 vertex != (*child)->vertex_.end();
                 ++vertex)
              childVertices.erase(*vertex);
          }
        }
      }
    }
  }
  // Mark the sub entities that are left for removal
  typedef typename std::set<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*>::iterator SFacetIter;
  for (SFacetIter e=childFacets.begin(); e!=childFacets.end(); ++e)
    (*e)->willVanish_=true;
  if(dimgrid > 1)
  {
    typedef typename std::set<FoamGridEntityImp<0, dimgrid, dimworld>*>::iterator VertexIter;
    for (VertexIter v=childVertices.begin(); v!=childVertices.end(); ++v)
      (*v)->willVanish_=true;
  }
}

// Refine one simplex element (2D simplex)
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::refineSimplexElement(FoamGridEntityImp<2, dimgrid, dimworld>& element,
                                                    int refCount)
{
  if(refCount<=0)
  {
    DUNE_THROW(NotImplemented, "Called refineSimplexElement with refCount <= 0");
    return;
  }

  unsigned int nextLevel=element.level()+1;

  array<FoamGridEntityImp<0, dimgrid, dimworld>*, 3*dimgrid> nextLevelVertices;
  std::size_t vertexIndex=0;

  // create copies of the vertices of the element
  for(int c=0; c<element.corners(); ++c)
  {
    if(element.vertex_[c]->sons_[0]==nullptr){
      // Vertex doesn't exist yet on the next level
      std::get<0>(entityImps_[nextLevel])
      .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel,
                                               element.vertex_[c]->pos_,
                                               element.vertex_[c]->id_));
      FoamGridEntityImp<0, dimgrid, dimworld>& newVertex =
      std::get<0>(entityImps_[nextLevel]).back();
      element.vertex_[c]->sons_[0]=&newVertex;
      newVertex.father_ = element.vertex_[c];
    }
    nextLevelVertices[vertexIndex++]=element.vertex_[c]->sons_[0];
  }

  // create new vertices from facet-midpoints together with the new facets that
  // have a father
  typedef typename array<FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>*, dimgrid+1>::iterator FacetIterator;

  array<FoamGridEntityImp<1, dimgrid, dimworld>*, 9> nextLevelFacets;
  std::size_t facetIndex=0;
#ifndef NDEBUG
  const auto& refElement = Dune::ReferenceElements<double, dimgrid>::general(element.type());
#endif

  // I am just to dumb for a general facet to vertice mapping.
  // Therefore we just store it here
  array<std::pair<unsigned int,unsigned int>,3 > facetVertexMapping;
  facetVertexMapping[0]=std::make_pair(0,1);
  facetVertexMapping[1]=std::make_pair(2,0);
  facetVertexMapping[2]=std::make_pair(1,2);

  for(FacetIterator facet=element.facet_.begin(); facet != element.facet_.end(); ++facet)
  {
#ifndef NDEBUG
    const auto* v0 = element.vertex_[refElement.subEntity(facetIndex/2, 1, 0, 2)];
    const auto* v1 = element.vertex_[refElement.subEntity(facetIndex/2, 1, 1, 2)];
#endif

    if(!(*facet)->nSons_)
    {
      // Not refined yet
      // Compute facet midpoint
      FieldVector<double, dimworld> midPoint;
      for(int dim=0; dim<dimworld;++dim)
        midPoint[dim]=((*facet)->vertex_[0]->pos_[dim]
        + (*facet)->vertex_[1]->pos_[dim]) /2.0;

      // if the element is parametrized we obtain the new point by the parametrization
      if(element.elementParametrization_)
      {
        // we know the local coordinates of the midpoint
        FoamGridEntity<0, dimgrid, Dune::FoamGrid<dimgrid, dimworld> > e(&element);
        FieldVector<double, dimgrid> localMidPoint(e.geometry().local(midPoint));
        while(e.hasFather())
        {
          localMidPoint = e.geometryInFather().global(localMidPoint);
          e.target_ = e.target_->father_;
        }
        // overwrite the midpoint with the coordinates mapped by the parametrization
        element.elementParametrization_->evaluate(localMidPoint, midPoint);
      }

      //create midpoint
      std::get<0>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel, midPoint,
                                               freeIdCounter_[0]++));
      FoamGridEntityImp<0, dimgrid, dimworld>& midVertex =
        std::get<0>(entityImps_[nextLevel]).back();
      nextLevelVertices[vertexIndex++]=&midVertex;

      // sanity check for DUNE numbering
      assert(v0->sons_[0]!=nullptr);
      assert(v1->sons_[0]!=nullptr);
      assert(v0->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].first] ||
        v0->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].second]);

      assert(v1->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].first] ||
        v1->sons_[0] == nextLevelVertices[facetVertexMapping[facetIndex/2].second]);

      // create the facets and publish them in the father
      std::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[facetVertexMapping[facetIndex/2].first], &midVertex,
                                                 nextLevel, freeIdCounter_[1]++, *facet));
      (*facet)->sons_[0] = &std::get<1>(entityImps_[nextLevel]).back();
      ++((*facet)->nSons_);
      // Inherit the boundaryId and SegmentIndex (are treated as the same for now)
      (*facet)->sons_[0]->boundaryId_=(*facet)->boundaryId_;
      (*facet)->sons_[0]->boundarySegmentIndex_=(*facet)->boundarySegmentIndex_;
      nextLevelFacets[facetIndex++]= (*facet)->sons_[0];

      // Initialize the elements_ vector of the new facet
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*facet)->sons_[0]->elements_=(*facet)->elements_;

      assert((*facet)->vertex_[1]->sons_[0]!=nullptr);
      std::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(&midVertex, nextLevelVertices[facetVertexMapping[facetIndex/2].second],
                                                 nextLevel, freeIdCounter_[1]++, *facet));
      (*facet)->sons_[1] = &std::get<1>(entityImps_[nextLevel]).back();
      ++((*facet)->nSons_);
      // Inherit the boundaryId and SegmentIndex (are treated as the same for now)
      (*facet)->sons_[1]->boundaryId_=(*facet)->boundaryId_;
      (*facet)->sons_[1]->boundarySegmentIndex_=(*facet)->boundarySegmentIndex_;
      nextLevelFacets[facetIndex++]= (*facet)->sons_[1];

      // Initialize the elements_ vector of the new facet
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*facet)->sons_[1]->elements_=(*facet)->elements_;

      (*facet)->nSons_=2;
    } else {
      // Facets do already exist. Just add its sons to nextLevelFacets
      // but make sure that the one containing vertex facetIndex comes first
      if((*facet)->sons_[0]->vertex_[0]->id_ ==
        nextLevelVertices[facetVertexMapping[facetIndex/2].first]->id_ ||
        (*facet)->sons_[0]->vertex_[1]->id_ ==
        nextLevelVertices[facetVertexMapping[facetIndex/2].first]->id_)
      {
        nextLevelFacets[facetIndex++]=(*facet)->sons_[0];
      nextLevelFacets[facetIndex++]=(*facet)->sons_[1];
      }else{
        nextLevelFacets[facetIndex++]=(*facet)->sons_[1];
        nextLevelFacets[facetIndex++]=(*facet)->sons_[0];
      }
      if((*facet)->sons_[0]->vertex_[0]->id_!=(*facet)->vertex_[0]->id_ &&
         (*facet)->sons_[0]->vertex_[0]->id_!=(*facet)->vertex_[1]->id_)
      {
        //vertex 0 is the midpoint
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>((*facet)->sons_[0]->vertex_[0]);
      } else {
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>((*facet)->sons_[0]->vertex_[1]);
      }
    }
  }
  assert(facetIndex==6);
  // Create the facets that lie within the father element
  // first the one that lies opposite to the vertex 0 in the father
  std::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[4], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&std::get<1>(entityImps_[nextLevel]).back();

  // the one opposite to father vertex 1
  std::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&std::get<1>(entityImps_[nextLevel]).back();

  // and the one opposite to father vertex 2
  std::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1, dimgrid, dimworld>(nextLevelVertices[4],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelFacets[facetIndex++]=&std::get<1>(entityImps_[nextLevel]).back();

  assert(facetIndex==nextLevelFacets.size());
  assert(vertexIndex==nextLevelVertices.size());

  array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 4> nextLevelElements;
  // create the new triangles that lie in the corners
  // First the one that contains vertex 0 of the father.
  std::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));

  FoamGridEntityImp<dimgrid, dimgrid, dimworld>* newElement = &(std::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[0];
  newElement->facet_[1]=nextLevelFacets[3];
  newElement->facet_[2]=nextLevelFacets[6];
  newElement->vertex_[0]=nextLevelVertices[0];
  newElement->vertex_[1]=nextLevelVertices[3];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=0;
  nextLevelElements[0]=newElement;
  element.sons_[0]=newElement;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  // Next the one that contains vertex 1 of the father.
  std::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(std::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[4];
  newElement->facet_[1]=nextLevelFacets[1];
  newElement->facet_[2]=nextLevelFacets[7];
  newElement->vertex_[0]=nextLevelVertices[1];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[3];
  newElement->refinementIndex_=1;
  nextLevelElements[1]=newElement;
  element.sons_[1]=newElement;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  // Last the one that contains vertex 2 of the father.
  std::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid ,dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(std::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[2];
  newElement->facet_[1]=nextLevelFacets[5];
  newElement->facet_[2]=nextLevelFacets[8];
  newElement->vertex_[0]=nextLevelVertices[2];
  newElement->vertex_[1]=nextLevelVertices[4];
  newElement->vertex_[2]=nextLevelVertices[5];
  newElement->refinementIndex_=2;
  nextLevelElements[2]=newElement;
  element.sons_[2]=newElement;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  // create the triangle in the center
  std::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(std::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->facet_[0]=nextLevelFacets[7];
  newElement->facet_[1]=nextLevelFacets[6];
  newElement->facet_[2]=nextLevelFacets[8];
  newElement->vertex_[0]=nextLevelVertices[3];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=3;
  nextLevelElements[3]=newElement;
  element.sons_[3]=newElement;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  // Now that all the triangle are created, we can update the elements attached
  // to the facets.
  // The new (inside) neighbors of the facets lying on facets of the father element.
  std::size_t neighbors[6] = {0, 1, 2, 0, 1, 2};
  for(std::size_t i=0; i<6; ++i){
    // Overwrite the father element by the newly created elements.
    overwriteFineLevelNeighbours(*nextLevelFacets[i], nextLevelElements[neighbors[i]],
                                 &element);
  }

  // Update the neighbours of the inner facets
  nextLevelFacets[6]->elements_.push_back(nextLevelElements[0]);
  nextLevelFacets[6]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[7]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[7]->elements_.push_back(nextLevelElements[1]);
  nextLevelFacets[8]->elements_.push_back(nextLevelElements[3]);
  nextLevelFacets[8]->elements_.push_back(nextLevelElements[2]);

  element.nSons_=4;

  if((refCount--)>1)
  {
    typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid>::iterator ElementIterator;
    for(ElementIterator elem=nextLevelElements.begin();
        elem != nextLevelElements.end(); ++elem)
    {
      refineSimplexElement(**elem, refCount);
    }
  }
}

// Refine one simplex element (1D simplex element)
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::refineSimplexElement(FoamGridEntityImp<1, dimgrid, dimworld>& element,
                                                    int refCount)
{
  if(refCount<=0)
  {
    DUNE_THROW(NotImplemented, "Called refineSimplexElement with refCount <= 0");
    return;
  }

  unsigned int nextLevel=element.level()+1;

  array<FoamGridEntityImp<0, dimgrid, dimworld>*, 3*dimgrid> nextLevelVertices;
  std::size_t vertexIndex=0;

  // create copies of the vertices of the element
  for(int c=0; c<element.corners(); ++c)
  {
    if(element.vertex_[c]->sons_[0]==nullptr){
      // Vertex doesn't exist yet on the next level
      std::get<0>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel,
                                               element.vertex_[c]->pos_,
                                               element.vertex_[c]->id_));
      FoamGridEntityImp<0, dimgrid, dimworld>& newVertex =
        std::get<0>(entityImps_[nextLevel]).back();

      // publish vertex in the father
      element.vertex_[c]->sons_[0] = &newVertex;
      element.vertex_[c]->nSons_++;
      assert(element.vertex_[c]->nSons_==1); // Vertex can't have more than one son
      // publish father in the new vertex
      newVertex.father_ = element.vertex_[c];
      assert(newVertex.hasFather());

      // Inherit the boundaryId_ and the boundarySegmentIndex
      element.vertex_[c]->sons_[0]->boundaryId_= element.vertex_[c]->boundaryId_;
      element.vertex_[c]->sons_[0]->boundarySegmentIndex_= element.vertex_[c]->boundarySegmentIndex_;

      // Initialize the elements_ vector if the new facet with that of the father. Later,
      // we will overwrite it with the correct values
      element.vertex_[c]->sons_[0]->elements_= element.vertex_[c]->elements_;
    }
    //add vertex to nextLevelVertices
    nextLevelVertices[vertexIndex++]=element.vertex_[c]->sons_[0];
  }
  assert(vertexIndex==2);

  // Create the facet/vertex which lies within the father element
  // Compute element midpoint
  typedef FoamGridEntityImp<0, dimgrid, dimworld> FoamGridVertex;
  const FoamGridVertex* v0 = element.vertex_[0];
  const FoamGridVertex* v1 = element.vertex_[1];
  FieldVector<double, dimworld> midPoint;
  for(int dim=0; dim<dimworld;++dim)
    midPoint[dim]=(v0->pos_[dim] + v1->pos_[dim])*0.5;

  // if the element is parametrized we obtain the new point by the parametrization
  if(element.elementParametrization_)
  {
    // we know the local coordinates of the midpoint
    FoamGridEntity<0, dimgrid, Dune::FoamGrid<dimgrid, dimworld> > e(&element);
    FieldVector<double, dimgrid> localMidPoint(0.5);
    while(e.hasFather())
    {
      localMidPoint = e.geometryInFather().global(localMidPoint);
      e.target_ = e.target_->father_;
    }
    // overwrite the midpoint with the coordinates mapped by the parametrization
    element.elementParametrization_->evaluate(localMidPoint, midPoint);
  }
  // Create element midpoint
  std::get<0>(entityImps_[nextLevel]).push_back(FoamGridEntityImp<0, dimgrid, dimworld>(nextLevel, midPoint, freeIdCounter_[0]++));
  FoamGridEntityImp<0, dimgrid, dimworld>& midVertex = std::get<0>(entityImps_[nextLevel]).back();
  nextLevelVertices[vertexIndex++]=&midVertex;

  assert(vertexIndex==nextLevelVertices.size()); //==3

  // Create next level elements
  array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid > nextLevelElements;
  // create the elements
  // First the one that contains vertex 0 of the father.
  std::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[dimgrid]++));

  FoamGridEntityImp<dimgrid, dimgrid, dimworld>* newElement = &(std::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->vertex_[0]=nextLevelVertices[0];
  newElement->vertex_[1]=nextLevelVertices[2];
  newElement->facet_[0]=nextLevelVertices[0]; // are equal to vertices but are
  newElement->facet_[1]=nextLevelVertices[2]; // set for consistent interface
  newElement->refinementIndex_=0;
  nextLevelElements[0]=newElement;
  element.sons_[0]=newElement;
  element.nSons_++;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  // Next the one that contains vertex 1 of the father.
  std::get<dimgrid>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<dimgrid, dimgrid, dimworld>(nextLevel, freeIdCounter_[dimgrid]++));
  newElement = &(std::get<dimgrid>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->vertex_[0]=nextLevelVertices[2];
  newElement->vertex_[1]=nextLevelVertices[1];
  newElement->facet_[0]=nextLevelVertices[2]; // are equal to vertices but are
  newElement->facet_[1]=nextLevelVertices[1]; // set for consistent interface
  newElement->refinementIndex_=1;
  nextLevelElements[1]=newElement;
  element.sons_[1]=newElement;
  element.nSons_++;

  // if the father is parametrized the son will be parametrized too
  if(element.elementParametrization_)
    newElement->elementParametrization_ = element.elementParametrization_;

  assert(element.nSons_== 1<<dimgrid); //==2

  // Now that all the elements are created, we can update the elements attached
  // to the facets.
  // The new (inside) neighbors of the facets lying on facets of the father element.
  std::size_t neighbors[2] = {0, 1};
  for(std::size_t i=0; i<2; ++i)
  {
    // Overwrite the father element by the newly created elements.
    overwriteFineLevelNeighbours(*nextLevelVertices[i], nextLevelElements[neighbors[i]],
                                    &element);
  }
  // Update the neighbours of the inner vertex
  nextLevelVertices[2]->elements_.push_back(nextLevelElements[0]);
  nextLevelVertices[2]->elements_.push_back(nextLevelElements[1]);

  if((refCount--)>1)
  {
    typedef typename array<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, 1<<dimgrid >::iterator ElementIterator;
    for(ElementIterator elem=nextLevelElements.begin();
        elem != nextLevelElements.end(); ++elem)
    {
      refineSimplexElement(**elem, refCount);
    }
  }
}

// Overwrites the neighbours of this and descendant facets
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::overwriteFineLevelNeighbours(FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>& facet,
                                                            const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* son,
                                                            const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* father)
{
  typedef typename std::vector<const FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>::iterator ElementIterator;

  for(ElementIterator elem=facet.elements_.begin();
      elem != facet.elements_.end();
      ++elem)
  {
    // father is replaced by the son in the elements_ vector of the vertex
    if (*elem == father)
    {
      *elem = son;
    }
  }

  for(std::size_t i=0; i<facet.nSons_; ++i)
    overwriteFineLevelNeighbours(*facet.sons_[i], son, father);
}


// Recompute the grid indices after the grid has changed
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::setIndices()
{
  // //////////////////////////////////////////
  //   Create the index sets
  // //////////////////////////////////////////
  for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
    // add space for new LevelIndexSets. They are not created until requested
    levelIndexSets_.push_back((FoamGridLevelIndexSet< const FoamGrid > *) 0);
  }

  // Delete old LevelIndexSets if the grid hierarchy got lower
  int excess = levelIndexSets_.size() - (maxLevel() + 1);
  for (int i=0; i<excess; i++) {
    if (levelIndexSets_.back())
      delete(levelIndexSets_.back());
    levelIndexSets_.pop_back();
  }

  for (int i=0; i<=maxLevel(); i++)
    if (levelIndexSets_[i])
      levelIndexSets_[i]->update();

  // Update the leaf indices
  leafIndexSet_.update();
}

//f Book-keeping routine to be called before growth
// Returns true if an element will be removed
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::preGrow()
{ return !elementsToRemove_.empty(); }


// Triggers the grid growth
// Returns true if new elements were created
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::grow()
{
  bool newEntities = false;
  bool removedEntities = false;

  // find a possible level for new elements and vertices, if there is more than one
  // possibility to insert the element/vertex we choose the one with the lowest level
  // first, find the min and max level for each element
  std::map<FoamGridEntityImp<0, dimgrid, dimworld>*, int> minVertexLevel, maxVertexLevel;
  std::map<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*, int> minElementLevel, maxElementLevel;
  for (auto eIt = elementsToInsert_.begin(); eIt != elementsToInsert_.end(); ++eIt)
  {
    // initialize the maps with numerical limits
    minElementLevel[&(*eIt)] = std::numeric_limits<int>::min();
    maxElementLevel[&(*eIt)] = std::numeric_limits<int>::max();
    for(auto vIt = eIt->vertex_.begin(); vIt != eIt->vertex_.end(); ++vIt)
    {
      // get pointer from the iterator
      FoamGridEntityImp<0, dimgrid, dimworld>* vertex = *vIt;

      // only leaf vertices can be chosen for growth so this vertex has
      // the maximum level if it is an already existing vertex
      int level = vertex->level_;
      maxVertexLevel[vertex] = vertex->isNew_ ? std::numeric_limits<int>::max() : level;

      // search for the lowest level this vertex exists on
      while(vertex->hasFather())
      {
        --level;
        vertex = vertex->father_;
      }
      minVertexLevel[vertex] = vertex->isNew_ ? std::numeric_limits<int>::min() : level;

      // the min and max element level is the intersection of possible element and vertex levels
      minElementLevel[&(*eIt)] = std::max(minElementLevel[&(*eIt)], minVertexLevel[vertex]);
      maxElementLevel[&(*eIt)] = std::min(maxElementLevel[&(*eIt)], maxVertexLevel[vertex]);
    }
    // if the element is not possible (to insert) we want to erase it (not add it)
    // in order to avoid erasing it now we just change the isNew_ marker
    // later only elements with the marker get added, all others discarded
    if(minElementLevel[&(*eIt)] > maxElementLevel[&(*eIt)])
    {
      eIt->isNew_ = false;
    }
  }

  // Now check if new vertices are shared (otherwise element level is minElementLevel)
  for(auto vIt = verticesToInsert_.begin(); vIt != verticesToInsert_.end(); ++vIt)
  {
    for (auto eIt = elementsToInsert_.begin(); eIt != elementsToInsert_.end(); ++eIt)
      for(auto otherVIt = eIt->vertex_.begin(); otherVIt != eIt->vertex_.end(); ++otherVIt)
        if(vIt->id_ == (*otherVIt)->id_ && vIt->level_ == (*otherVIt)->level_)
        {
          assert((*otherVIt)->isNew_);
          minVertexLevel[&(*vIt)] = std::max(minElementLevel[&(*eIt)], minVertexLevel[&(*vIt)]);
          maxVertexLevel[&(*vIt)] = std::min(maxElementLevel[&(*eIt)], maxVertexLevel[&(*vIt)]);
        }

    // if it's not possible to insert the vertex the connecting elements are also not possible
    if(minVertexLevel[&(*vIt)] > maxVertexLevel[&(*vIt)])
    {
      // erase elements with pointers to the impossible new vertex
      for (auto eIt = elementsToInsert_.begin(); eIt != elementsToInsert_.end(); ++eIt)
      {
        for(auto otherVIt = eIt->vertex_.begin(); otherVIt != eIt->vertex_.end(); ++otherVIt)
        {
          if(vIt->id_ == (*otherVIt)->id_ && vIt->level_ == (*otherVIt)->level_)
          {
            assert((*otherVIt)->isNew_);
            eIt->isNew_ = false;
            break;
          }
        }
      }
    }
  }

  // All elements that are left now can be validly inserted in the maximum of minElementLevel and
  // minVertexLevel of their respective new vertices
  using VertexPointer = FoamGridEntityImp<0, dimgrid, dimworld>*;
  std::map<VertexPointer, VertexPointer> insertedMap;
  for (auto eIt = elementsToInsert_.begin(); eIt != elementsToInsert_.end(); ++eIt)
  {
    // skip impossible elements determined above
    if(!eIt->isNew_)
      continue;

    // find the element level
    for(auto vIt = eIt->vertex_.begin(); vIt != eIt->vertex_.end(); ++vIt)
      if((*vIt)->isNew_)
        minElementLevel[&(*eIt)] = std::max(minElementLevel[&(*eIt)], minVertexLevel[*vIt]);

    // set the element level
    const int level = std::max(minElementLevel[&(*eIt)], 0);
    eIt->level_ = level;

    // set the vertex levels
    for(auto vIt = eIt->vertex_.begin(); vIt != eIt->vertex_.end(); ++vIt)
    {
      if((*vIt)->isNew_)
      {
        // set the vertex level
        (*vIt)->level_ = level;
        // Insert the new vertex into the grid if it wasn't already inserted
        if(!insertedMap.count(*vIt))
        {
          // add the vertex to the grid
          (*vIt)->id_ = freeIdCounter_[0]++;
          std::get<0>(entityImps_[level]).push_back(*(*vIt));

          // add a map entry so we know this new vertex was already inserted
          insertedMap[*vIt] = &*std::get<0>(entityImps_[level]).rbegin();

          // publish actual vertex pointer in element
          (*vIt) = &*std::get<0>(entityImps_[level]).rbegin();
        }
        else
        {
          // the vertex was already inserted, only publish the vertex pointer in element
          (*vIt) = insertedMap[*vIt];
        }
      }
      // if the vertex was already there, replace it with the vertex on the right level
      else
      {
        while((*vIt)->level_ != level)
          (*vIt) = (*vIt)->father_;
      }
    }

    // We are ready to insert the element into the grid
    eIt->id_ = freeIdCounter_[dimgrid]++;
    std::get<dimgrid>(entityImps_[level]).push_back(*eIt);
    newEntities = true;
  }

  // cleanup
  growing_ = false;
  verticesToInsert_.clear();
  elementsToInsert_.clear();
  indexToVertexMap_.clear();

  if(newEntities)
  {
    // Now we deal with the facets
    // Existing facets have to get knowledge of the new element, non-existing ones have to be added
    // Construct a map from vertex arrays to facets
    std::vector<std::map<std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid>,
                                 FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* > > facetMap;
    for(int level = 0; level <= maxLevel(); level++)
    {
      std::map<std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid>, FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* > tmp;
      facetMap.push_back(tmp);
      for(auto fIt = std::get<dimgrid-1>(entityImps_[level]).begin(); fIt != std::get<dimgrid-1>(entityImps_[level]).end(); ++fIt)
      {
        // directly using the facet's vertex information is not possible because it doesn't exist in 1d
        // we have to get the vertex pointers from a neighbouring element
        // every old vertex has at least one element
        if(fIt->elements_.size() > 0)
        {
          auto element = fIt->elements_[0];
          const Dune::ReferenceElement<double,dimgrid>& refElement
                      = Dune::ReferenceElements<double,dimgrid>::general(element->type());

          // obtain the local index of this facet
          std::size_t localFacetIndex = 0;
          for (std::size_t fIdx = 0; fIdx < element->facet_.size(); ++fIdx)
            if(&*fIt == element->facet_[fIdx])
              localFacetIndex = fIdx;

          // construct the vertex array of the facet
          std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray;
          for (std::size_t vIdx = 0; vIdx < dimgrid; ++vIdx)
            vertexArray[vIdx] = const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>
                                  (element->vertex_[refElement.subEntity(localFacetIndex, 1, vIdx, dimgrid)]);
          // add map entry
          facetMap[level][vertexArray] = &*fIt;
        }
      }
    }

    // facet map for each level, we could only do the levels that new elements get inserted on...
    for(int level = 0; level <= maxLevel(); level++)
    {
      // loop over all elements on this level
      for(auto eIt = std::get<dimgrid>(entityImps_[level]).begin(); eIt != std::get<dimgrid>(entityImps_[level]).end(); ++eIt)
      {
        const Dune::ReferenceElement<double,dimgrid>& refElement
                  = Dune::ReferenceElements<double,dimgrid>::general(eIt->type());
        // loop over all facets of this element
        for (std::size_t fIdx = 0; fIdx < eIt->facet_.size(); ++fIdx)
        {
          // get the vertices of the facet
          std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray;
          std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> flippedVertexArray;
          for (std::size_t vIdx = 0; vIdx < dimgrid; ++vIdx)
          {
            if(dimgrid == 2)
            {
              vertexArray[vIdx] = eIt->vertex_[refElement.subEntity(fIdx, dimgrid-1, vIdx, dimgrid)];
              flippedVertexArray[1-vIdx] = eIt->vertex_[refElement.subEntity(fIdx, dimgrid-1, vIdx, dimgrid)];
            }
            if(dimgrid == 1)
              vertexArray[vIdx] = eIt->vertex_[fIdx];
          }

          FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* existingFacet = nullptr;
          auto f = facetMap[level].find(vertexArray);

          if(f != facetMap[level].end())
          {
            // facet is already in map
            existingFacet = f->second;
          }
          else
          {
            // in 2d also check the permutation of the array
            if(dimgrid == 2)
            {
              f = facetMap[level].find(flippedVertexArray);
              if(f != facetMap[level].end()) existingFacet = f->second;
            }
          }
          if(existingFacet == nullptr)
          {
            // I'm too dumb to find a general algorithm working in 1d and 2d so we use function overloads here
            addNewFacet(existingFacet, vertexArray, level);
            // put the facet in the map
            facetMap[level][vertexArray] = existingFacet;
          }

          // make element know about the facet
          eIt->facet_[fIdx] = existingFacet;

          // make facet know about the element if it's new
          // if the facet is not a leaf facet we have to make sure all sons of it
          // until the leaf level have the new element in their element vector
          if(eIt->isNew_)
            addElementForFacet(&*eIt, existingFacet);
        }
      }
    }
  }

  // Set the element's (to be removed) flag to willVanish_ = true
  // Remark: Vertices to delete can only be deleted efficiently with codim 2-0 connectivity information (touching point grid problem)
  // This is why we calculate the connectivity first
  if(dimgrid == 2 && (!elementsToRemove_.empty()))
    computeTwoZeroConnectivity();

  for(auto&& ep : elementsToRemove_)
    removedEntities = removeSimplexElement(*ep) || removedEntities;

  // cleanup
  elementsToRemove_.clear();

  if(removedEntities)
  {
    // actually remove the entities
    for(int level = 0; level <= maxLevel(); level++)
    {
      eraseVanishedEntities(std::get<0>(entityImps_[level]));
      if(dimgrid == 2)
      {
        eraseVanishedEntities(std::get<dimgrid-1>(entityImps_[level]));
      }
      eraseVanishedEntities(std::get<dimgrid>(entityImps_[level]));
    }
  }

  if(removedEntities || newEntities)
  {
    // set boundary indices (it is a completely new set -> need for boundary IDs?? (see FS#1369) to transfer boundary data?)
    unsigned int boundaryFacetCounter = 0;
    for(int level = 0; level <= maxLevel(); level++)
      for (auto&& facet : std::get<dimgrid-1>(entityImps_[level]))
        if(facet.isLeaf() && facet.elements_.size()==1) //if boundary facet
          facet.boundarySegmentIndex_ = boundaryFacetCounter++;
    numBoundarySegments_ = boundaryFacetCounter;
  }

  // update the leaf index if something happened
  if(removedEntities || newEntities)
    leafIndexSet_.update();

  return newEntities;
}

// helper function to add an element to the element vectors of a facet's sons
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::addElementForFacet(const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* element,
                                                           FoamGridEntityImp<dimgrid-1, dimgrid, dimworld>* facet)
{
  // publish element in the facet
  facet->elements_.push_back(element);
  // and it's sons (recursively) if it's not on the leaf
  if(!facet->isLeaf())
    for(auto&& son : facet->sons_)
      addElementForFacet(element, son);
}

// helper function to add new facets in 1d (they already exist as vertices)
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::addNewFacet(FoamGridEntityImp<0, dimgrid, dimworld>* &facet,
                                                    std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray,
                                                    int level)
{ facet = vertexArray[0]; }

// helper function to add new facets in 2d
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::addNewFacet(FoamGridEntityImp<1, dimgrid, dimworld>* &facet,
                                                    std::array<FoamGridEntityImp<0, dimgrid, dimworld>*,dimgrid> vertexArray,
                                                    int level)
{
  std::get<1>(entityImps_[level]).push_back(FoamGridEntityImp<1, 2, dimworld>(vertexArray[0], vertexArray[1], level, freeIdCounter_[1]++));
  facet = &*std::get<1>(entityImps_[level]).rbegin();
}

// Clean up refinement markers
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::postGrow()
{
  // Loop over all leaf entities and remove the isNew Marker
  // and invalidate the growth inserstion index
  // and set the coarseing blocker in case we created a T-junction
  for (auto&& element : elements(this->leafGridView()))
  {
    FoamGridEntityImp<dimgrid, dimgrid, dimworld>& e = *const_cast<FoamGridEntityImp<dimgrid, dimgrid, dimworld>*>(this->getRealImplementation(element).target_);
    e.isNew_=false;
    e.growthInsertionIndex_=-1;

    for(auto&& vertex : e.vertex_)
      vertex->isNew_ = false;

    //Block elements that do now have a facet being a junction and not having a father for coarsening
    for(auto&& facet : e.facet_)
      if(facet->elements_.size() > 2 && !(facet->hasFather()))
        e.coarseningBlocked_ = true;

    // Also block all children of the father as one child marked for coarsening is enough to coarsen all
    if(e.coarseningBlocked_ && e.hasFather())
    {
      FoamGridEntityImp<dimgrid, dimgrid, dimworld>& father=*(e.father_);
      for (auto&& child : father.sons_)
        child->coarseningBlocked_ = true;
    }
  }
}

// Remove an element from the grid at runtime
template <int dimgrid, int dimworld>
bool Dune::FoamGrid<dimgrid, dimworld>::removeSimplexElement(FoamGridEntityImp<dimgrid, dimgrid, dimworld>& element)
{
  element.willVanish_ = true;
  // check which facets of the element only belong to this element and have to be removed too
  for(auto&& facet : element.facet_)
  {
    if(facet->elements_.size() < 2)
    {
      facet->willVanish_ = true;
      if(facet->hasFather())
      {
        --(facet->father_->nSons_);
        for(auto&& son : facet->father_->sons_)
          if(son->willVanish_)
            son = nullptr;
      }
    }
    else
    {
      // this facet will become a boundary through the removal if it only has two associated elements
      if(facet->elements_.size() == 2)
      {
        facet->boundaryId_= numBoundarySegments_;
        facet->boundarySegmentIndex_= numBoundarySegments_;
        ++numBoundarySegments_;
      }
      // remove the entity pointers of the remaining subentities and the father
      for(auto e = facet->elements_.begin(); e != facet->elements_.end();)
      {
        if((*e)->willVanish_)
          e = facet->elements_.erase(e);
        else
          ++e;
      }
    }
  }
  if(element.hasFather())
  {
    --(element.father_->nSons_);
    for(auto&& son : element.father_->sons_)
      if(son->willVanish_)
        son = nullptr;
  }
  // in 2d check for vertices to be deleted
  // vertices to delete can only be deleted efficiently with codim 2-0 connectivity information (touching point grid problem)
  if(dimgrid == 2)
  {
    for(auto&& vertex : element.vertex_)
    {
      if(vertex->elements_.size() < 2)
      {
        vertex->willVanish_ = true;
        if(vertex->hasFather())
        {
          vertex->father_->nSons_ = 0;
          vertex->father_->sons_[0] = nullptr;
        }
      }
    }
  }

  return true;
}

// change a vertex' position on all levels
template <int dimgrid, int dimworld>
void Dune::FoamGrid<dimgrid, dimworld>::setPosition(const typename Traits::template Codim<dimgrid>::Entity & e,
                                                    const FieldVector<ctype, dimworld>& pos)
{
  auto vertex = const_cast<FoamGridEntityImp<0, dimgrid, dimworld>*>(this->getRealImplementation(e).target_);

  if (!vertex->isLeaf())
    DUNE_THROW(Dune::NotImplemented, "Moving vertices that are not on the leaf!");

  vertex->pos_ = pos;

  // set position of possible fathers on coarser levels
  auto father = vertex;
  while (father->hasFather())
  {
    father = father->father_;
    father->pos_ = pos;
  }
}
