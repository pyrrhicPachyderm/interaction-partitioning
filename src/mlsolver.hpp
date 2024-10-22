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
		
		virtual Parameters getSolutionParameters(bool isNull) = 0;
		virtual double getSolutionAdditionalParameter(bool isNull, size_t i) = 0;
		
		virtual Eigen::VectorXd getSolutionPredictions(bool isNull) = 0;
		virtual Eigen::VectorXd getSolutionResiduals(bool isNull) = 0;
		
		virtual double getDeviance(bool isNull) = 0;
		virtual size_t getNumParameters(bool isNull) = 0;
		virtual double getAIC(bool isNull) = 0;
		virtual double getAICc(bool isNull) = 0;
		virtual double getR2(bool isNull) = 0;
		virtual double getMcFaddenR2() = 0;
		
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
		ParametersT nullSolution;
		bool isDirtySolution = true;
		bool isDirtyNullSolution = true;
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
		virtual void calculateSolution(bool isNull) = 0;
		ParametersT getSolution(bool isNull);
		
		Eigen::VectorXd parametersToVector(const ParametersT &parameters) const;
		ParametersT vectorToParameters(const Eigen::VectorXd &vector) const;
	public:
		Parameters getSolutionParameters(bool isNull) override;
		double getSolutionAdditionalParameter(bool isNull, size_t i) override;
		
		Eigen::VectorXd getSolutionPredictions(bool isNull) override;
		Eigen::VectorXd getSolutionResiduals(bool isNull) override;
		
		virtual double getDeviance(bool isNull) = 0;
		size_t getNumParameters(bool isNull) override;
		double getAIC(bool isNull) override;
		double getAICc(bool isNull) override;
		double getR2(bool isNull) override;
		double getMcFaddenR2() override;
		
		virtual MaximumLikelihoodSolverInterface *getCopy() const = 0;
};

class GaussNewtonSolver : public MaximumLikelihoodSolver<Solver> {
	public:
		using MaximumLikelihoodSolver<Solver>::MaximumLikelihoodSolver;
	protected:
		Eigen::VectorXd getResidualsFromVector(const Eigen::VectorXd &parameterVector) const;
		Jacobian getResidualsJacobianFromVector(const Eigen::VectorXd &parameterVector) const;
		
		void calculateSolution(bool isNull) override;
	public:
		double getDeviance(bool isNull) override;
		
		MaximumLikelihoodSolverInterface *getCopy() const override {return new GaussNewtonSolver(*this);};
};

template<typename ErrDistT> class NLoptSolver : public MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>> {
	public:
		using MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>>::MaximumLikelihoodSolver;
		
		typedef MaximumLikelihoodSolver<GeneralisedSolver<ErrDistT>>::ParametersT ParametersT;
	protected:
		
		//NLopt requires a function pointer, not the std::function produced by std::bind.
		//A std::function cannot be converted to a function pointer; function pointers must be to free functions, not to methods.
		//So we define a static (free) function as a wrapper, which takes this object as a void pointer.
		double getLogLikelihoodFromVector(const std::vector<double> &parametersVector);
		static double optimisationFunc(const std::vector<double>& parametersVector, std::vector<double>& grad, void* solver);
		
		void calculateSolution(bool isNull) override;
	public:
		double getDeviance(bool isNull) override;
		
		MaximumLikelihoodSolverInterface *getCopy() const override {return new NLoptSolver<ErrDistT>(*this);};
};

#endif
