#ifndef MAXIMUM_LIKELIHOOD_SOLVER_HPP
#define MAXIMUM_LIKELIHOOD_SOLVER_HPP

#include "solver.hpp"

class MaximumLikelihoodSolverInterface {
	public:
		//Virtual destructor, as this is an abstract class.
		virtual ~MaximumLikelihoodSolverInterface() {};
		
		virtual bool updateGrouping(GroupingType groupingType, bool (Grouping::*updateFunc)()) = 0;
		virtual void setGroupings(const GroupingSet &groupings) = 0;
		virtual void setGrouping(GroupingType groupingType, const Grouping &grouping) = 0;
		virtual const GroupingSet &getGroupings() const = 0;
		virtual const Grouping &getGrouping(GroupingType groupingType) const = 0;
		
		virtual Parameters getSolutionParameters() = 0;
		virtual double getSolutionAdditionalParameter(size_t i) = 0;
		
		virtual Eigen::VectorXd getSolutionPredictions() = 0;
		virtual Eigen::VectorXd getSolutionResiduals() = 0;
		
		virtual double getDeviance() = 0;
		virtual size_t getNumParameters() = 0;
		virtual double getAIC() = 0;
		virtual double getAICc() = 0;
		virtual double getR2() = 0;
		
		//For parallelisation purposes, we sometimes need copies of the entire solver.
		//We need this as a virtual function, so it can correctly copy each specialisation of the solver.
		virtual MaximumLikelihoodSolverInterface *getCopy() const = 0;
};

template<typename SolverT> class MaximumLikelihoodSolver : public SolverT, public MaximumLikelihoodSolverInterface {
	static_assert(std::is_base_of_v<Solver, SolverT>);
	protected:
		GroupingSet groupings;
	public:
		typedef SolverT::ParametersT ParametersT;
		
		MaximumLikelihoodSolver(Model model, Data data):
			SolverT(model, data), groupings({Grouping(data.getNumRowSpecies()), Grouping(data.getNumRowSpecies()), Grouping(data.getNumColSpecies())}) {};
		MaximumLikelihoodSolver(Model model, Data data, GroupingSet groupings):
			SolverT(model, data),
			groupings(groupings)
		{
			for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
				assert(groupings[i].numSpecies == data.getNumSpecies((GroupingType)i));
			}
		};
	protected:
		ParametersT solution;
		bool isDirtySolution = true;
		bool isDirtyNullSolution = true; //Not all specailisations of this class use this, but it must be maintained here nonetheless.
	public:
		//A function to update groupings.
		//Can use reset(), separate(), or advance().
		bool updateGrouping(GroupingType groupingType, bool (Grouping::*updateFunc)()) override {
			isDirtySolution = true;
			isDirtyNullSolution = true;
			return (groupings[groupingType].*updateFunc)();
		}
		void setGroupings(const GroupingSet &gs) override {
			isDirtySolution = true;
			isDirtyNullSolution = true;
			groupings = gs;
		}
		void setGrouping(GroupingType groupingType, const Grouping &g) override {
			isDirtySolution = true;
			isDirtyNullSolution = true;
			groupings[groupingType] = g;
		}
		
		const GroupingSet &getGroupings() const override {return groupings;}
		const Grouping &getGrouping(GroupingType groupingType) const override {return groupings[groupingType];}
	protected:
		virtual void calculateSolution() = 0;
		ParametersT getSolution();
		
		Eigen::VectorXd parametersToVector(const ParametersT &parameters) const;
		ParametersT vectorToParameters(const Eigen::VectorXd &vector) const;
	public:
		Parameters getSolutionParameters() override;
		double getSolutionAdditionalParameter(size_t i) override;
		
		Eigen::VectorXd getSolutionPredictions() override;
		Eigen::VectorXd getSolutionResiduals() override;
		
		virtual double getDeviance() = 0;
		size_t getNumParameters() override;
		double getAIC() override;
		double getAICc() override;
		virtual double getR2() = 0;
		
		virtual MaximumLikelihoodSolverInterface *getCopy() const = 0;
};

class GaussNewtonSolver : public MaximumLikelihoodSolver<Solver> {
	public:
		using MaximumLikelihoodSolver<Solver>::MaximumLikelihoodSolver;
	protected:
		Eigen::VectorXd getResidualsFromVector(const Eigen::VectorXd &parameterVector) const;
		Jacobian getResidualsJacobianFromVector(const Eigen::VectorXd &parameterVector) const;
		
		void calculateSolution() override;
	public:
		double getR2() override;
		double getDeviance() override;
		
		MaximumLikelihoodSolverInterface *getCopy() const override {return new GaussNewtonSolver(*this);};
};

template<typename ErrDistT> class NLoptSolver : public MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>> {
	public:
		using MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>>::MaximumLikelihoodSolver;
		
		typedef MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>>::ParametersT ParametersT;
	protected:
		ParametersT nullSolution;
		
		//NLopt requires a function pointer, not the std::function produced by std::bind.
		//A std::function cannot be converted to a function pointer; function pointers must be to free functions, not to methods.
		//So we define a static (free) function as a wrapper, which takes this object as a void pointer.
		double getLogLikelihoodFromVector(const std::vector<double> &parametersVector);
		static double optimisationFunc(const std::vector<double>& parametersVector, std::vector<double>& grad, void* solver);
		
		ParametersT solve(bool isNull) const;
		void calculateSolution() override;
		void calculateNullSolution();
		ParametersT getNullSolution();
	public:
		double getR2() override;
		double getDeviance() override;
		
		MaximumLikelihoodSolverInterface *getCopy() const override {return new NLoptSolver<ErrDistT>(*this);};
};

#endif
