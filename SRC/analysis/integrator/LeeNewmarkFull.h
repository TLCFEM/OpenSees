/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#ifndef LeeNewmarkFull_h
#define LeeNewmarkFull_h

#include "Newmark.h"
#include <vector>
#include <sparseGEN/sp_sparse/triplet_form.hpp>

class LeeNewmarkFull final : public Newmark {
public:
	// declare which damping model type
	enum class Type { T0, T1, T2, T3, T4 };

	// Define a structure Basis to store model parameters
	struct Basis {
		Type t;
		std::vector<double> p;
		double wp, zp;
	};

private:
	// initial vector (after Basis) to store each damping basis
	std::vector<Basis> damping_basis;

	// declare three sparse matrices ub a triplet form
	triplet_form<double> rabbit; // temporary matrix to store the top-left matrix of original n-dof size
	const triplet_form<double> current_stiffness;
	const triplet_form<double> current_mass;

	// ---------------------------------------------------------------------------------------------------

	unsigned get_amplifier() const;  // this is to obtain the size?
	unsigned get_total_size() const; // this is to obtain the total big matrix size
	void update_residual() const;    // this is to update the residual (RHS) based on the big matrix

	// The following lines are functions to assemble the big matrix based on damping model type.
	void assemble_Bcx(triplet_form<double>&, unsigned&, double, double, int&) const;
	void assemble_by_type_zero(triplet_form<double>&, unsigned&, double, double) const;
	void assemble_by_type_one(triplet_form<double>&, unsigned&, double, double, int) const;
	void assemble_by_type_two(triplet_form<double>&, unsigned&, double, double, int, int) const;
	void assemble_by_type_three(triplet_form<double>&, unsigned&, double, double, double) const;
	void assemble_by_type_four(triplet_form<double>&, unsigned&, double, double, int, int, int, int, double) const;

	// an enumerated list (MatType) to check for which matrix to compute.  OpenSees only has one matrix.  This implementation needs to have 5 matrices: K_initial at t_0, K_current at t_n, and K_tangent at t_{n+1}, M, and C
	enum class MatType : unsigned { None, InitialStiffness, CurrentStiffness, TangentStiffness, Mass, Damping };

	// Initialize a variable (which_matrix) as none first.
	MatType which_matrix = MatType::None;

	// a constant (protected from being changed) variable (stiffness_type), initialized at K_current
	const MatType stiffness_type = MatType::CurrentStiffness;

	// initialize a variable as 0
	unsigned n_block = 0;

	// a logical variable to check if it is first_iteration
	bool first_iteration = true;

	// declare two vectors
	Vector current_internal, trial_internal;

public:
	LeeNewmarkFull(double,               // alpha
	               double,               // beta
	               std::vector<Basis>&&, // damping basis vector
	               unsigned              // flag to indicate which stiffness to be used
	);

	int formTangent(int) override;

	int formEleTangent(FE_Element*) override;

	int formNodTangent(DOF_Group*) override;

	int commit() override;
	int revertToLastStep() override;
	int revertToStart() override;

	int domainChanged() override;

	int update(const Vector&) override;
};

void* OPS_LeeNewmarkFull(unsigned);

template<typename T> std::enable_if_t<!std::numeric_limits<T>::is_integer, bool> approx_equal(T x, T y, int ulp = 2) { return fabs(x - y) <= std::numeric_limits<T>::epsilon() * fabs(x + y) * ulp || fabs(x - y) < std::numeric_limits<T>::min(); }

#endif
