#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"
#include "parameters.hpp"

class Solver {
	protected:
		Data data;
	public:
		Solver(Data data): data(data) {};
	protected:
		Eigen::MatrixXd colGroupedDesign;
		bool isDirtyColGroupedDesign = true;
		
		//We will have functions to say that particular elements have changed and mark appropriate things as dirty.
		//There needs to be a set for the superclass, Solver, as well as ones for subclasses.
		//The subclass ones will be define pure virtual here, and will need overriding.
		//The superclass ones will dirty elements of the superclass, then call the subclass ones.
		virtual void dirtyDataSubclass() = 0;
		virtual void dirtyGroupingSubclass(GroupingType groupingType) = 0;
		
		void dirtyData() {
			isDirtyColGroupedDesign = true;
			dirtyDataSubclass();
		}
		void dirtyGrouping(GroupingType groupingType) {
			if(groupingType == COL) isDirtyColGroupedDesign = true;
			dirtyGroupingSubclass(groupingType);
		}
		
		//Again, this should be templated on the enum, but can't be.
		//So we add a level of indirection to expose a templated interface, for consistency.
		virtual const Grouping &getGroupingSubclass(GroupingType groupingType) const = 0;
	public:
		template<GroupingType groupingType> const Grouping &getGrouping() const {
			return getGroupingSubclass(groupingType);
		}
	protected:
		void calculateColGroupedDesign();
		Eigen::MatrixXd getColGroupedDesign();
		
		Eigen::VectorXd getPredictions(const Parameters &parameters);
		Eigen::VectorXd getResiduals(const Parameters &parameters);
};

#endif
