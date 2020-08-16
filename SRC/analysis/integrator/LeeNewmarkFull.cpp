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

// ReSharper disable CppClangTidyClangDiagnosticShorten64To32
// ReSharper disable CppClangTidyClangDiagnosticSignConversion
#include "LeeNewmarkFull.h"
#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <elementAPI.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <iostream>
#include <LinearSOE.h>
#include <sparseGEN/LeeSparse.h>
#include <string.h>

/* Paper references
C&S1: Lee (2021) Adjustable bandwidth
C&S2: Lee (2022) Softening response
*/

// Read user inputs for damping model parameters and Newmark coefficients
void* OPS_LeeNewmarkFull(const unsigned stiffness_type) {
	opserr << "Draft implementation of LeeNewmarkFull, version 1\n";

	double gamma = -1., beta = -1.;

	int num_to_read = 1;

	OPS_GetDoubleInput(&num_to_read, &gamma);

	OPS_GetDoubleInput(&num_to_read, &beta);

	std::vector<LeeNewmarkFull::Basis> damping_basis;

	double wp, zp;

	while(OPS_GetNumRemainingInputArgs() > 0) {
		// get type identifier
		// if failed break;
		const char* type_raw_string = OPS_GetString();

		OPS_GetDoubleInput(&num_to_read, &wp);
		OPS_GetDoubleInput(&num_to_read, &zp);

		// get wp and zp
		if(strcmp(type_raw_string, "-type0") == 0) {
			// type0

			damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T0, std::vector<double>{}, wp, zp});
		}
		else if(strcmp(type_raw_string, "-type1") == 0) {
			// type1

			int para_a;
			OPS_GetIntInput(&num_to_read, &para_a); // n_p

			if(0 == para_a) {
				opserr << "warning: Type 1 changed to Type 0 due to trivial np.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T0, std::vector<double>{}, wp, zp});
			}
			else damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T1, std::vector<double>{static_cast<double>(para_a)}, wp, zp});
		}
		else if(strcmp(type_raw_string, "-type2") == 0) {
			// type2

			int para_a, para_b;
			OPS_GetIntInput(&num_to_read, &para_a); // n_pr
			OPS_GetIntInput(&num_to_read, &para_b); // n_pl

			if(0 == para_a + para_b) {
				opserr << "warning: Type 2 changed to Type 0 due to trivial npr, npl.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T0, std::vector<double>{}, wp, zp});
			}
			else damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T2, std::vector<double>{static_cast<double>(para_a), static_cast<double>(para_b)}, wp, zp});
		}
		else if(strcmp(type_raw_string, "-type3") == 0) {
			// type3

			double para_a;
			OPS_GetDoubleInput(&num_to_read, &para_a); // gamma

			if(para_a < -1.) {
				opserr << "warning: gamma shall be strictly greater than -1.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T3, std::vector<double>{para_a}, wp, zp});
			}
			else if(approx_equal(para_a, 0.)) {
				opserr << "warning: Type 3 changed to Type 0 due to trivial gamma.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T0, std::vector<double>{}, wp, zp});
			}
			else damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T3, std::vector<double>{para_a}, wp, zp});
		}
		else if(strcmp(type_raw_string, "-type4") == 0) {
			// type4

			int para_a, para_b, para_c, para_d;
			double para_e;
			OPS_GetIntInput(&num_to_read, &para_a);    // n_pr
			OPS_GetIntInput(&num_to_read, &para_b);    // n_pl
			OPS_GetIntInput(&num_to_read, &para_c);    // n_pk
			OPS_GetIntInput(&num_to_read, &para_d);    // n_pm
			OPS_GetDoubleInput(&num_to_read, &para_e); // gamma

			if(para_e <= -1.) {
				opserr << "warning: gamma shall be strictly greater than -1.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T4, std::vector<double>{static_cast<double>(para_a), static_cast<double>(para_b), static_cast<double>(para_c), static_cast<double>(para_d), para_e}, wp, zp});
			}
			else if(approx_equal(para_e, 0.)) {
				if(0 == para_a + para_b) {
					opserr << "warning: Type 4 changed to Type 0 due to trivial gamma and npr, npl.\n";
					damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T0, std::vector<double>{}, wp, zp});
				}
				else {
					opserr << "warning: Type 4 changed to Type 2 due to trivial gamma.\n";
					damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T2, std::vector<double>{static_cast<double>(para_a), static_cast<double>(para_b)}, wp, zp});
				}
			}
			else if(0 == para_a + para_b + para_c + para_d) {
				opserr << "warning: Type 4 changed to Type 3 due to trivial npr, npl, npk, npm.\n";
				damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T3, std::vector<double>{para_e}, wp, zp});
			}
			else damping_basis.emplace_back(LeeNewmarkFull::Basis{LeeNewmarkFull::Type::T4, std::vector<double>{static_cast<double>(para_a), static_cast<double>(para_b), static_cast<double>(para_c), static_cast<double>(para_d), para_e}, wp, zp});
		}
		else {
			opserr << "unknown flag\n";
			return nullptr;
		}
	}

	if(damping_basis.empty()) return nullptr;

	return new LeeNewmarkFull(gamma, beta, std::move(damping_basis), stiffness_type);
}

// Estimate the number of Ks and Ms in each damping basis matrix CbX
unsigned LeeNewmarkFull::get_amplifier() const {
	auto n_size = 2u;

	for(auto& I : damping_basis)
		if(Type::T0 == I.t) n_size += 5u;
		else if(Type::T1 == I.t) n_size += 5u + 6u * static_cast<unsigned>(I.p[0]);
		else if(Type::T2 == I.t) n_size += 4u + 5u * static_cast<unsigned>(.5 * (I.p[0] + I.p[1] - 1.));
		else if(Type::T3 == I.t) n_size += 9u;
		else if(Type::T4 == I.t) n_size += 2u * static_cast<unsigned>(I.p[0] + I.p[1] + I.p[2] + I.p[3] + 4.);

	return n_size;
}

// Estimate the matrix size (number of rows / columns) in the total big system matrix
unsigned LeeNewmarkFull::get_total_size() const {
	auto n_size = 1u;

	// add additional size, NOT the size of each mode
	for(auto& I : damping_basis)
		if(Type::T0 == I.t) n_size += 1u;
		else if(Type::T1 == I.t) n_size += 2u * static_cast<unsigned>(I.p[0]) + 1u;
		else if(Type::T2 == I.t) n_size += static_cast<unsigned>(I.p[0] + I.p[1]) + 1u;
		else if(Type::T3 == I.t) n_size += 2u;
		else if(Type::T4 == I.t) n_size += static_cast<unsigned>(I.p[0] + I.p[1] + I.p[2] + I.p[3] + 2.);

	return n_size * n_block;
}

// Update the residual considering the presence of big damping C
void LeeNewmarkFull::update_residual() const {
	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	auto& residual = t_soe->residual;

	// 1. get damping force due to big C from matrix-vector multiplication
	auto trial_vel = (-1) * trial_internal;
	for(auto I = 0u; I < n_block; ++I) trial_vel(I) = (*Udot)(I) / -c2;

	// now residual holds the damping force solely due to global damping model big C
	residual = t_soe->global_stiffness * trial_vel;

	// 2. get original residual due to K, M and C.
	auto& original_residual = t_soe->getB();

	for(auto I = 0u; I < n_block; ++I) residual(I) += original_residual(I);
}

void LeeNewmarkFull::assemble_Bcx(triplet_form<double>& stiffness, unsigned& current_pos, const double mass_coef, const double stiffness_coef, int& order) const {
	const auto ic1 = current_pos;   // first block
	const auto ic2 = ic1 + n_block; // second block
	const auto& s_row = current_stiffness.row_idx;
	const auto& s_col = current_stiffness.col_idx;
	const auto& s_val = current_stiffness.val_idx;

	if(order > 0) // for standard cases with order = 1, 2, etc.
	{
		// fill the big matrix (called stiffness) with two off-diagonal Kc matrices (first and second blocks)
		for(index_t Ki = 0; Ki < current_stiffness.c_size; ++Ki) {
			const index_t i1 = ic1 + s_row[Ki], j1 = ic1 + s_col[Ki]; // first block in the extended
			const index_t i2 = i1 + n_block, j2 = j1 + n_block;       // second block in the extended
			stiffness.at(i1, j2) = stiffness.at(i2, j1) = stiffness_coef * s_val[Ki];
		}

		// fill the big matrix with Mc
		const auto& m_row = current_mass.row_idx;
		const auto& m_col = current_mass.col_idx;
		const auto& m_val = current_mass.val_idx;

		if(order > 1) // for order = 2, 3, etc
		{
			// fill the big matrix (called stiffness) with two off-diagonal Mc matrices (second and third blocks)
			for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
				const index_t i2 = ic2 + m_row[Mi], j2 = ic2 + m_col[Mi]; // second block in the extended
				const index_t i3 = i2 + n_block, j3 = j2 + n_block;       // third block in the extended
				stiffness.at(i2, j3) = stiffness.at(i3, j2) = mass_coef * m_val[Mi];
			}
		}
		else {
			// for order = 1

			// fill the big matrix (called stiffness) with only one negative Mc matrix in the second block diagonal
			for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
				const index_t i2 = ic2 + m_row[Mi], j2 = ic2 + m_col[Mi]; // second block in the extended
				stiffness.at(i2, j2) = -mass_coef * m_val[Mi];
			}
		}

		// update the current position, and update order.
		current_pos += 2 * n_block;
		order -= 2;
	}
	else {
		// for order = 0, note Bcx(0) = Kc

		// fill the big matrix (called stiffness) with one Kc matrix in the diagonal (first block)
		for(index_t Ki = 0; Ki < current_stiffness.c_size; ++Ki) {
			const index_t i1 = ic1 + s_row[Ki], j1 = ic1 + s_col[Ki]; // first block in the extended
			stiffness.at(i1, j1) = stiffness_coef * s_val[Ki];
		}
		// update the current position with only one block shift, and update order.
		current_pos += n_block;
		order -= 1;
	}
}

// Assemble the new system matrix for C based on type 0 model
void LeeNewmarkFull::assemble_by_type_zero(triplet_form<double>& stiffness, unsigned& current_pos, double mass_coef, double stiffness_coef) const {
	const auto I = current_pos;
	mass_coef = 4. * mass_coef;
	stiffness_coef = 4. * stiffness_coef;

	auto row = current_mass.row_idx;
	auto col = current_mass.col_idx;
	auto val = current_mass.val_idx;
	for(index_t J = 0; J < current_mass.c_size; ++J) {
		const index_t K = row[J], L = col[J], M = K + I, N = L + I;
		stiffness.at(K, N) = stiffness.at(M, L) = -(stiffness.at(K, L) = stiffness.at(M, N) = mass_coef * val[J]);
	}
	row = current_stiffness.row_idx;
	col = current_stiffness.col_idx;
	val = current_stiffness.val_idx;
	for(index_t J = 0; J < current_stiffness.c_size; ++J) stiffness.at(row[J] + I, col[J] + I) = stiffness_coef * val[J];

	current_pos += n_block;
}

// Assemble the new system matrix for C based on type 1 model
void LeeNewmarkFull::assemble_by_type_one(triplet_form<double>& stiffness, unsigned& current_pos, const double mass_coef, const double stiffness_coef, int order) const {
	const auto mass_coef1 = 4. * mass_coef;
	const auto stiffness_coef1 = 4. * stiffness_coef;
	const auto mass_coefs = .5 * mass_coef1;           // eq. 10
	const auto stiffness_coefs = .5 * stiffness_coef1; // eq. 10

	// current mass and stiffness are now stored as objects so use dot (.) instead of ->
	const auto& m_row = current_mass.row_idx;
	const auto& m_col = current_mass.col_idx;
	const auto& m_val = current_mass.val_idx;
	const auto& s_row = current_stiffness.row_idx;
	const auto& s_col = current_stiffness.col_idx;
	const auto& s_val = current_stiffness.val_idx;

	index_t I = 0;
	index_t J = current_pos;
	index_t K = current_pos += n_block;
	index_t L = current_pos += n_block;
	index_t M = current_pos += n_block;

	while(order > 1) {
		// eq. 61
		for(index_t N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
			stiffness.at(O + I, R) = stiffness.at(Q, P + I) = mass_coef1 * m_val[N];
			stiffness.at(Q, P + K) = stiffness.at(O + K, R) = stiffness.at(O + L, P + M) = stiffness.at(O + M, P + L) = mass_coefs * m_val[N];
		}

		for(index_t N = 0; N < current_stiffness.c_size; ++N) {
			const auto O = s_row[N], P = s_col[N], Q = O + K, R = O + L, S = P + K, T = P + L;
			stiffness.at(Q, T) = stiffness.at(R, S) = stiffness_coef1 * s_val[N];
			stiffness.at(O + J, S) = stiffness.at(Q, P + J) = stiffness.at(R, P + M) = stiffness.at(O + M, T) = stiffness_coefs * s_val[N];
		}

		I = current_pos;
		J = current_pos += n_block;
		K = current_pos += n_block;
		L = current_pos += n_block;
		M = current_pos += n_block;
		order -= 2;
	}

	if(order < 1) {
		for(index_t N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + I, R = P + I, S = O + J, T = P + J;
			stiffness.at(Q, R) = stiffness.at(S, T) = stiffness.at(Q, T) = stiffness.at(S, R) = mass_coef1 * m_val[N];
		}

		for(index_t N = 0; N < current_stiffness.c_size; ++N) stiffness.at(s_row[N] + J, s_col[N] + J) = stiffness_coef1 * s_val[N];

		current_pos = K;

		return;
	}

	// eq. 53
	for(index_t N = 0; N < current_mass.c_size; ++N) {
		const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
		stiffness.at(O + L, P + L) = -(stiffness.at(O + I, R) = stiffness.at(Q, P + I) = mass_coef1 * m_val[N]);
		stiffness.at(Q, P + K) = stiffness.at(O + K, R) = mass_coefs * m_val[N];
	}

	for(index_t N = 0; N < current_stiffness.c_size; ++N) {
		const auto O = s_row[N], P = s_col[N], Q = O + K, R = P + K, S = O + L, T = P + L;
		stiffness.at(S, T) = -(stiffness.at(Q, T) = stiffness.at(S, R) = stiffness_coef1 * s_val[N]);
		stiffness.at(O + J, R) = stiffness.at(Q, P + J) = stiffness_coefs * s_val[N];
	}
}

// Assemble the new system matrix for C based on type 2 model
void LeeNewmarkFull::assemble_by_type_two(triplet_form<double>& stiffness, unsigned& current_pos, const double mass_coef, const double stiffness_coef, const int npr, int npl) const {
	const auto r = (2. * static_cast<double>(npl) + 1.) / (2. * static_cast<double>(npr) + 1.);
	const auto nps = static_cast<double>(npr) + static_cast<double>(npl) + 1.;
	const auto alpha = 2. * (1. + r) / pow(r, (static_cast<double>(npl) + 1.) / nps);
	const auto beta = alpha * pow(r, 1. / nps);
	const auto mass_coef2 = alpha * mass_coef;
	const auto stiffness_coef2 = beta * stiffness_coef;

	// current mass and stiffness are now stored as objects so use dot (.) instead of ->
	const auto& m_row = current_mass.row_idx;
	const auto& m_col = current_mass.col_idx;
	const auto& m_val = current_mass.val_idx;

	const auto ic1 = current_pos; // first block

	if(npr > 0) // standard, refers to eq. 81 of C&S1
	{
		const auto Bcxrsize = npr * n_block;

		// Add Mc2 in off-diagonal
		for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
			const auto i0 = m_row[Mi], j0 = m_col[Mi];         // for the zeroth block in the original equation
			const auto i1 = ic1 + i0, j1 = ic1 + j0;           // for the first block in the extended
			const auto i2 = i1 + Bcxrsize, j2 = j1 + Bcxrsize; // for the second block in the extended
			stiffness.at(i0, j1) = stiffness.at(i1, j0) = stiffness.at(i1, j2) = stiffness.at(i2, j1) = mass_coef2 * m_val[Mi];
		}

		// Add -Bcx(npr-1) in the diagonal
		auto npr1 = npr - 1;
		while(npr1 > -1) assemble_Bcx(stiffness, current_pos, -mass_coef2, -stiffness_coef2, npr1);

		// Add Bcx(npl) in the diagonal
		while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef2, stiffness_coef2, npl);
	}
	else // npr == 0, refer eq. 100 of C&S2, with the last block row and column removed.
	{
		//top left corner 2x2 mass block
		// Add Mc2
		for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
			const auto i0 = m_row[Mi], j0 = m_col[Mi]; // for the zeroth block in the original equation
			const auto i1 = ic1 + i0, j1 = ic1 + j0;   // for the first block in the extended
			stiffness.at(i0, j0) = stiffness.at(i0, j1) = stiffness.at(i1, j0) = stiffness.at(i1, j1) = mass_coef2 * m_val[Mi];
		}

		// Add Bcx(npl)
		while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef2, stiffness_coef2, npl);
	}
}

// Assemble the new system matrix for C based on type 3 model, refer to eq. 100 of C&S2
void LeeNewmarkFull::assemble_by_type_three(triplet_form<double>& stiffness, unsigned& current_pos, const double mass_coef, const double stiffness_coef, const double gammac) const {
	const auto gammab = 4. * gammac; // eq. 10
	const auto mass_coef3 = 4. * (1. + gammac) * mass_coef;
	const auto stiffness_coef3 = 4. * (1. + gammac) * stiffness_coef;
	const auto mass_coef3t = mass_coef3 * (1. + 1. / gammab);
	const auto stiffness_coef3t = stiffness_coef3 / gammab;

	// current mass and stiffness are now stored as objects so use dot (.) instead of ->
	const auto& m_row = current_mass.row_idx;
	const auto& m_col = current_mass.col_idx;
	const auto& m_val = current_mass.val_idx;
	const auto& s_row = current_stiffness.row_idx;
	const auto& s_col = current_stiffness.col_idx;
	const auto& s_val = current_stiffness.val_idx;

	const index_t ic1 = current_pos; // first block in the extended

	// Add Mc3
	for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
		index_t i0 = m_row[Mi], j0 = m_col[Mi], i1 = ic1 + i0, j1 = ic1 + j0, i2 = i1 + n_block, j2 = j1 + n_block;
		stiffness.at(i0, j0) = stiffness.at(i0, j1) = stiffness.at(i0, j2) = stiffness.at(i1, j0) = stiffness.at(i1, j1) = stiffness.at(i1, j2) = stiffness.at(i2, j0) = stiffness.at(i2, j1) = mass_coef3 * m_val[Mi];
		stiffness.at(i2, j2) = mass_coef3t * m_val[Mi];
	}

	// Add Kc3
	for(index_t Ki = 0; Ki < current_stiffness.c_size; ++Ki) {
		index_t i1 = ic1 + s_row[Ki], j1 = ic1 + s_col[Ki], i2 = i1 + n_block, j2 = j1 + n_block;
		stiffness.at(i1, j1) = stiffness_coef3 * s_val[Ki];
		stiffness.at(i2, j2) = stiffness_coef3t * s_val[Ki];
	}

	// update the current position index
	current_pos += 2 * n_block;
}

// Assemble the new system matrix for C based on type 4 model
void LeeNewmarkFull::assemble_by_type_four(triplet_form<double>& stiffness, unsigned& current_pos, double mass_coef, double stiffness_coef, int npr, int npl, int npk, int npm, double gammac) const {
	// Type 4 C&S2
	const auto& m_row = current_mass.row_idx;
	const auto& m_col = current_mass.col_idx;
	const auto& m_val = current_mass.val_idx;

	// serial
	const auto rs = (2. * static_cast<double>(npl) + 1.) / (2. * static_cast<double>(npr) + 1.);
	const auto nps = static_cast<double>(npr) + static_cast<double>(npl) + 1.;
	const auto alphas = 2. * (1. + rs) / pow(rs, (static_cast<double>(npl) + 1.) / nps);
	const auto betas = alphas * pow(rs, 1. / nps);
	const auto mass_coef4s = alphas * (1. + gammac) * mass_coef;
	const auto stiffness_coef4s = betas * (1. + gammac) * stiffness_coef;
	// parallel
	const auto rp = (2. * static_cast<double>(npm) + 1.) / (2. * static_cast<double>(npk) + 1.);
	const auto npt = static_cast<double>(npk) + static_cast<double>(npm) + 1.;
	const auto alphap = 2. * (1. + rp) / pow(rp, (static_cast<double>(npm) + 1.) / npt);
	const auto betap = alphap * pow(rp, 1. / npt);
	const auto mass_coef4p = alphap * (1. + gammac) * mass_coef;
	const auto stiffness_coef4p = betap * (1. + gammac) * stiffness_coef;

	const auto gammab = .25 * alphap * betap * gammac;

	const auto ic1 = current_pos; // first block in the extended

	if(npr + npm == 0) // eq. 100
	{
		const auto ic2 = ic1 + (npl + 1) * n_block; // second block in the extended

		// Add Mc4s and Mc4p
		const auto mass_coef4t = mass_coef4s + mass_coef4p / gammab;
		for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
			index_t i0 = m_row[Mi], j0 = m_col[Mi], i1 = ic1 + i0, j1 = ic1 + j0, i2 = ic2 + i0, j2 = ic2 + j0;
			stiffness.at(i0, j0) = stiffness.at(i0, j1) = stiffness.at(i0, j2) = stiffness.at(i1, j0) = stiffness.at(i1, j1) = stiffness.at(i1, j2) = stiffness.at(i2, j0) = stiffness.at(i2, j1) = mass_coef4s * m_val[Mi];
			stiffness.at(i2, j2) = mass_coef4t * m_val[Mi];
		}

		// Add Bcx(npl)
		while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef4s, stiffness_coef4s, npl);

		// Add Bcx(npk)/gammab
		while(npk > -1) assemble_Bcx(stiffness, current_pos, mass_coef4p / gammab, stiffness_coef4p / gammab, npk);

		return;
	}

	if(npr == 0) // eq. 98 
	{
		const auto ic2 = ic1 + (npl + 1) * n_block; // second block in the extended
		const auto ic3 = ic2 + (npk + 1) * n_block; // third block in the extended

		// Add Mc4s and Mc4p
		for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
			index_t i0 = m_row[Mi], j0 = m_col[Mi], i1 = ic1 + i0, j1 = ic1 + j0, i2 = ic2 + i0, j2 = ic2 + j0, i3 = ic3 + i0, j3 = ic3 + j0;
			stiffness.at(i0, j0) = stiffness.at(i0, j1) = stiffness.at(i0, j2) = stiffness.at(i1, j0) = stiffness.at(i1, j1) = stiffness.at(i1, j2) = stiffness.at(i2, j0) = stiffness.at(i2, j1) = stiffness.at(i2, j2) = mass_coef4s * m_val[Mi];
			stiffness.at(i2, j3) = stiffness.at(i3, j2) = mass_coef4p * m_val[Mi];
		}

		// Add Bcx(npl)
		while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef4s, stiffness_coef4s, npl);

		// Add Bcx(npk)/gammab
		while(npk > -1) assemble_Bcx(stiffness, current_pos, mass_coef4p / gammab, stiffness_coef4p / gammab, npk);

		// Add -gammab*Bcx(npm-1)
		auto npm1 = npm - 1;
		while (npm1 > -1) assemble_Bcx(stiffness, current_pos, -gammab * mass_coef4p, -gammab * stiffness_coef4p, npm1);

		return;
	}

	if(npm == 0) // eq. 97
	{
		const auto ic2 = ic1 + npr * n_block;       // second block in the extended
		const auto ic3 = ic2 + (npl + 1) * n_block; // third block in the extended

		// Add Mc4s and Mc4p
		const auto mass_coef4pm = mass_coef4p / gammab;
		for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
			index_t i0 = m_row[Mi], j0 = m_col[Mi], i1 = ic1 + i0, j1 = ic1 + j0, i2 = ic2 + i0, j2 = ic2 + j0, i3 = ic3 + i0, j3 = ic3 + j0;
			stiffness.at(i0, j1) = stiffness.at(i1, j0) = stiffness.at(i1, j2) = stiffness.at(i2, j1) = stiffness.at(i1, j3) = stiffness.at(i3, j1) = mass_coef4s * m_val[Mi];
			stiffness.at(i3, j3) = mass_coef4pm * m_val[Mi];
		}

		// Add -Bcx(npr-1)
		auto npr1 = npr - 1;
		while(npr1 > -1) assemble_Bcx(stiffness, current_pos, -mass_coef4s, -stiffness_coef4s, npr1);

		// Add Bcx(npl)
		while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef4s, stiffness_coef4s, npl);

		// Add Bcx(npk)/gammab
		while(npk > -1) assemble_Bcx(stiffness, current_pos, mass_coef4p / gammab, stiffness_coef4p / gammab, npk);

		return;
	}

	// eq. 84
	const auto ic2 = ic1 + npr * n_block;       // second block in the extended
	const auto ic3 = ic2 + (npl + 1) * n_block; // third block in the extended
	const auto ic4 = ic3 + (npk + 1) * n_block; // fourth block in the extended

	// Add Mc4s and Mc4p
	for(index_t Mi = 0; Mi < current_mass.c_size; ++Mi) {
		index_t i0 = m_row[Mi], j0 = m_col[Mi], i1 = ic1 + i0, j1 = ic1 + j0, i2 = ic2 + i0, j2 = ic2 + j0, i3 = ic3 + i0, j3 = ic3 + j0, i4 = ic4 + i0, j4 = ic4 + j0;
		stiffness.at(i0, j1) = stiffness.at(i1, j0) = stiffness.at(i1, j2) = stiffness.at(i2, j1) = stiffness.at(i1, j3) = stiffness.at(i3, j1) = mass_coef4s * m_val[Mi];
		stiffness.at(i3, j4) = stiffness.at(i4, j3) = mass_coef4p * m_val[Mi];
	}

	// Add -Bcx(npr-1)
	auto npr1 = npr - 1;
	while(npr1 > -1) assemble_Bcx(stiffness, current_pos, -mass_coef4s, -stiffness_coef4s, npr1);

	// Add Bcx(npl)
	while(npl > -1) assemble_Bcx(stiffness, current_pos, mass_coef4s, stiffness_coef4s, npl);

	// Add Bcx(npk)/gammab
	while(npk > -1) assemble_Bcx(stiffness, current_pos, mass_coef4p / gammab, stiffness_coef4p / gammab, npk);

	// Add -gammab*Bcx(npm-1)
	auto npm1 = npm - 1;
	while(npm1 > -1) assemble_Bcx(stiffness, current_pos, -gammab * mass_coef4p, -gammab * stiffness_coef4p, npm1);
}

// constructor
// std::vector<Basis>&& is a r-value
LeeNewmarkFull::LeeNewmarkFull(const double _gamma, const double _beta, std::vector<Basis>&& _basis, const unsigned _stiffness_type)
	: Newmark(_gamma, _beta)
	, damping_basis(std::forward<std::vector<Basis>>(_basis))
	, stiffness_type(static_cast<MatType>(_stiffness_type)) {}

int LeeNewmarkFull::formTangent(const int F) {
	statusFlag = F;

	auto result = 0;

	auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());
	auto t_model = getAnalysisModel();

	if(t_soe == nullptr || t_model == nullptr) {
		opserr << "WARNING LeeNewmarkFull::formTangent() no LeeSparse or AnalysisModel has been set\n";
		return -1;
	}

	if(first_iteration || stiffness_type == MatType::TangentStiffness) {
		t_soe->zero_current_stiffness();
		t_soe->zero_mass();
	}
	t_soe->zero_stiffness();
	t_soe->zero_damping();

	auto& t_dof = t_model->getDOFs();
	DOF_Group* t_dofptr;

	while((t_dofptr = t_dof()) != nullptr) {
		if(first_iteration || stiffness_type == MatType::TangentStiffness) {
			which_matrix = stiffness_type;
			if(t_soe->add_current_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:dof\n";
				result = -1;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_mass:dof\n";
				result = -1;
			}
		}
		which_matrix = MatType::TangentStiffness;
		if(t_soe->add_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:dof\n";
			result = -1;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_damping:dof\n";
			result = -1;
		}
	}

	auto& t_ele = t_model->getFEs();
	FE_Element* t_eleptr;
	while((t_eleptr = t_ele()) != nullptr) {
		if(first_iteration || stiffness_type == MatType::TangentStiffness) {
			which_matrix = stiffness_type;
			if(t_soe->add_current_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:ele\n";
				result = -2;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_mass:ele\n";
				result = -2;
			}
		}
		which_matrix = MatType::TangentStiffness;
		if(t_soe->add_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:ele\n";
			result = -2;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_damping:ele\n";
			result = -2;
		}
	}

	/* After this line, all matrices (of original size) are assembled. Elemental work is done.
	 * Now we deal with the formulation of the big matrix.
	 */

	auto& stiffness = t_soe->global_stiffness;

	// global stiffness is now populated with global damping matrix
	// need to add K+M+C to top left corner

	if(first_iteration || stiffness_type == MatType::TangentStiffness) {
		// preallocate memory
		stiffness = triplet_form<double>(get_total_size(), get_total_size(), get_amplifier() * t_soe->stiffness.c_size);

		t_soe->mass.csc_condense();
		t_soe->current_stiffness.csc_condense();

		access::rw(current_mass) = t_soe->mass;
		access::rw(current_stiffness) = t_soe->current_stiffness;

		// now current mass and stiffness are formulated
		// assemble unrolled damping matrix and the corresponding damping force
		// assemble stiffness

		// populating global stiffness matrix with global damping matrix
		auto IDX = n_block;

		for(auto& I : damping_basis) {
			const auto mass_coef = I.zp * I.wp * c2;
			const auto stiffness_coef = I.zp / I.wp * c2;
			if(Type::T0 == I.t) assemble_by_type_zero(stiffness, IDX, mass_coef, stiffness_coef);
			else if(Type::T1 == I.t) assemble_by_type_one(stiffness, IDX, mass_coef, stiffness_coef, I.p[0]);
			else if(Type::T2 == I.t) assemble_by_type_two(stiffness, IDX, mass_coef, stiffness_coef, I.p[0], I.p[1]);
			else if(Type::T3 == I.t) assemble_by_type_three(stiffness, IDX, mass_coef, stiffness_coef, I.p[0]);
			else if(Type::T4 == I.t) assemble_by_type_four(stiffness, IDX, mass_coef, stiffness_coef, I.p[0], I.p[1], I.p[2], I.p[3], I.p[4]);
		}

		// now stiffness holds global damping matrix big C
		stiffness.csc_condense();

		// update residual according to global damping matrix
		update_residual();

		const auto& row = stiffness.row_idx;
		const auto& col = stiffness.col_idx;
		const auto& val = stiffness.val_idx;

		rabbit = triplet_form<double>(n_block, n_block, t_soe->stiffness.c_size);
		for(unsigned I = 0; I < stiffness.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// check in left top block of unrolled damping matrix to be used in subsequent iterations
			if(row[I] < n_block) rabbit.at(row[I], col[I]) = val[I];
		}

		// after this line, rabbit holds the top left corner of global damping matrix may or may not be zeros
		// depending on types

		// now add K, M, C to global effective stiffness matrix
		stiffness += t_soe->stiffness + t_soe->mass * c3 + t_soe->damping * c2;

		// no need to condense again
		// it will be done before solving automatically

		first_iteration = false;
	}
	else {
		// if not first iteration
		// erase the tangent stiffness entries
		if(!stiffness.csc_sort()) return -1;

		const auto& row = stiffness.row_idx;
		const auto& col = stiffness.col_idx;
		const auto& val = stiffness.val_idx;

		for(size_t I = 0; I < stiffness.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// erase existing entries if fall in intact stiffness matrix
			if(row[I] < n_block) val[I] = 0.;
		}

		// check in original nonzero entries in unrolled damping matrix
		stiffness += rabbit;

		// ! now global damping matrix is recovered and stored in stiffness
		update_residual();

		// formulate global effective stiffness matrix
		stiffness += t_soe->stiffness + t_soe->mass * c3 + t_soe->damping * c2;
	}

	// ! SuperLU does not accept zero diagonal.

	return result;
}

int LeeNewmarkFull::formEleTangent(FE_Element* theEle) {
	theEle->zeroTangent();

	if(MatType::InitialStiffness == which_matrix) theEle->addKiToTang();
	else if(MatType::TangentStiffness == which_matrix || MatType::CurrentStiffness == which_matrix)
		if(statusFlag == CURRENT_TANGENT) theEle->addKtToTang();
		else if(statusFlag == INITIAL_TANGENT) theEle->addKiToTang();
		else if(statusFlag == HALL_TANGENT) {
			theEle->addKtToTang(cFactor);
			theEle->addKiToTang(iFactor);
		}
		else { opserr << "Newmark::formEleTangent - unknown FLAG\n"; }
	else if(MatType::Mass == which_matrix) theEle->addMtoTang();
	else if(MatType::Damping == which_matrix) theEle->addCtoTang();

	return 0;
}

int LeeNewmarkFull::formNodTangent(DOF_Group* theDof) {
	theDof->zeroTangent();

	if(MatType::Mass == which_matrix) theDof->addMtoTang();
	else if(MatType::Damping == which_matrix) theDof->addCtoTang();

	return 0;
}

int LeeNewmarkFull::commit() {
	first_iteration = true;

	current_internal = trial_internal;

	return IncrementalIntegrator::commit();
}

int LeeNewmarkFull::revertToLastStep() {
	first_iteration = true;

	trial_internal = current_internal;

	return Newmark::revertToLastStep();
}

int LeeNewmarkFull::revertToStart() {
	first_iteration = true;

	trial_internal.Zero();
	current_internal.Zero();

	return Newmark::revertToStart();
}

// A function to update the size parameters whenever there is a domain changed.
int LeeNewmarkFull::domainChanged() {
	const auto flag = Newmark::domainChanged();

	if(flag != 0) return flag;

	n_block = getLinearSOE()->getX().Size();

	const auto n_size = get_total_size();

	first_iteration = true;

	current_internal = Vector(n_size);
	trial_internal = Vector(n_size);

	return 0;
}

// A function to update the increment DU and other internal variables
int LeeNewmarkFull::update(const Vector& DU) {
	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	trial_internal += t_soe->increment;

	return Newmark::update(DU);
}
