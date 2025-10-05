#pragma once
namespace Dune {

/*
 * A simple AMG block diagonal preconditioner
 */
template <class M, class X, class Y, int blockLevel = 2>
class BlockDiagAMGPreconditioner : public Dune::Preconditioner<X, Y> {
    template <std::size_t i>
    using DiagBlockType =
        std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template <std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template <std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template <std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template <std::size_t i>
    using ScalarProduct = Dune::SeqScalarProduct<VecBlockType<i>>;

    template <std::size_t i>
    using BlockAMG = Dune::Amg::AMG<LinearOperator<i>, VecBlockType<i>, Smoother<i>>;

    using AMGTuple = typename makeFromIndexedType<std::tuple, BlockAMG, std::make_index_sequence<M::size()>>::type;

   public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of.
        category = SolverCategory::sequential
    };

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The (multi type block) matrix to operate on.
       \param w The relaxation factor.
     */
    template <class LOP, class Criterion, class SmootherArgs>
    BlockDiagAMGPreconditioner(const LOP& lop, const Criterion& c, const SmootherArgs& sa)
        : BlockDiagAMGPreconditioner(lop, c, sa, std::make_index_sequence<M::size()>{}) {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre(X& v, Y& d) {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)),
                [&](const auto i) { std::get<decltype(i)::value>(amg_).pre(v[i], d[i]); });
    }

    /*!
     * \brief Apply the preconditoner.
     * \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply(X& v, const Y& d) {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)),
                [&](const auto i) { std::get<decltype(i)::value>(amg_).apply(v[i], d[i]); });
    }

    /*!
     * \brief Clean up.
     * \copydoc Preconditioner::post(X&)
     */
    virtual void post(X& v) {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)),
                [&](const auto i) { std::get<decltype(i)::value>(amg_).post(v[i]); });
    }

   private:
    template <class LOP, class Criterion, class SmootherArgs, std::size_t... Is>
    BlockDiagAMGPreconditioner(const LOP& lop, const Criterion& c, const SmootherArgs& sa,
                               std::index_sequence<Is...> is)
        : amg_(std::make_tuple(BlockAMG<Is>(*std::get<Is>(lop), *std::get<Is>(c), *std::get<Is>(sa))...)) {}

    AMGTuple amg_;
};

/*
 * \brief A simple AMG block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagAMGBiCGSTABSolver {
    template <class M, std::size_t i>
    using DiagBlockType =
        std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template <class X, std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template <class M, class X, std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

    template <class M, class X, std::size_t i>
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother<M, X, i>>::Arguments;

    template <class M, std::size_t i>
    using Criterion =
        Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DiagBlockType<M, i>, Dune::Amg::FirstDiagonal>>;

    template <class M, class X, std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

   public:
    // Solve saddle-point problem using a Schur complement based preconditioner
    template <int precondBlockLevel = 2, class Matrix, class Vector>
    bool solve(const Matrix& m, Vector& x, const Vector& b, double residualReduction, int maxIter, int verbosity) {
        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        //! \todo make parameters changeable at runtime from input file / parameter tree
        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(3);
        params.setDebugLevel(verbosity);

        auto criterion = makeCriterion_<Criterion, Matrix>(params, std::make_index_sequence<Matrix::size()>{});
        auto smootherArgs = makeSmootherArgs_<SmootherArgs, Matrix, Vector>(std::make_index_sequence<Matrix::size()>{});

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(m)), [&](const auto i) {
            auto& args = std::get<decltype(i)::value>(smootherArgs);
            args->iterations = 1;
            args->relaxationFactor = 1;
        });

        auto linearOperator =
            makeLinearOperator_<LinearOperator, Matrix, Vector>(m, std::make_index_sequence<Matrix::size()>{});

        BlockDiagAMGPreconditioner<Matrix, Vector, Vector> preconditioner(linearOperator, criterion, smootherArgs);

        Dune::MatrixAdapter<Matrix, Vector, Vector> op(m);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, residualReduction, maxIter, verbosity);
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const { return result_; }

    std::string name() const { return "block-diagonal AMG-preconditioned BiCGSTAB solver"; }

   private:
    template <template <class M, std::size_t i> class Criterion, class Matrix, class Params, std::size_t... Is>
    auto makeCriterion_(const Params& p, std::index_sequence<Is...>) {
        return std::make_tuple(std::make_shared<Criterion<Matrix, Is>>(p)...);
    }

    template <template <class M, class X, std::size_t i> class SmootherArgs, class Matrix, class Vector,
              std::size_t... Is>
    auto makeSmootherArgs_(std::index_sequence<Is...>) {
        return std::make_tuple(std::make_shared<SmootherArgs<Matrix, Vector, Is>>()...);
    }

    template <template <class M, class X, std::size_t i> class LinearOperator, class Matrix, class Vector,
              std::size_t... Is>
    auto makeLinearOperator_(const Matrix& m, std::index_sequence<Is...>) {
        return std::make_tuple(std::make_shared<LinearOperator<Matrix, Vector, Is>>(
            m[Dune::index_constant<Is>{}][Dune::index_constant<Is>{}])...);
    }

    Dune::InverseOperatorResult result_;
};

}  // end namespace Dune
