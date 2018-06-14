#include "../GMGAvgRstr.h"
#include "../GMGDrctIntp.h"
#include "../InterLevelComm.h"
#include "catch.hpp"
using namespace std;
const int n = 8;
// generate 2 level simple test
void generateTwoLevel(shared_ptr<DomainCollection> &coarse, shared_ptr<DomainCollection> &fine)
{
	// generate simple tree
	OctTree t("3uni.bin");
	coarse.reset(new DomainCollection(t, 2, n));
	fine.reset(new DomainCollection(t, n));
	fine->zoltanBalance();
	coarse->zoltanBalanceWithLower(*fine);
}
void octFill(double *vec, int oct, double val)
{
	switch (oct) {
		case 0:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi) + (yi) *n + (zi) *n * n] = val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 1:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi + n / 2) + (yi) *n + (zi) *n * n] = val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 2:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi) + (yi + n / 2) * n + (zi) *n * n] = val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 3:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi + n / 2) + (yi + n / 2) * n + (zi) *n * n]
						= val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 4:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi) + (yi) *n + (zi + n / 2) * n * n] = val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 5:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi + n / 2) + (yi) *n + (zi + n / 2) * n * n]
						= val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 6:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi) + (yi + n / 2) * n + (zi + n / 2) * n * n]
						= val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		case 7:
			for (int zi = 0; zi < n / 2; zi++) {
				for (int yi = 0; yi < n / 2; yi++) {
					for (int xi = 0; xi < n / 2; xi++) {
						vec[(xi + n / 2) + (yi + n / 2) * n + (zi + n / 2) * n * n]
						= val + xi + yi * n + zi * n * n;
					}
				}
			}
			break;
		default:
			break;
	}
}
TEST_CASE("InterLevelComm scatter works", "[GMG]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<InterLevelComm> comm(new InterLevelComm(coarse, fine));
		PW<Vec>                    coarse_expected       = coarse->getNewDomainVec();
		PW<Vec>                    coarse_local_expected = comm->getNewCoarseDistVec();
		double *                   ce_vec;
		VecGetArray(coarse_expected, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int i = 0; i < n * n * n; i++) {
				ce_vec[d.id_local * n * n * n + i] = d.id + i;
			}
		}
		VecRestoreArray(coarse_expected, &ce_vec);
		double coarse_expected_norm;
		VecNorm(coarse_expected, NORM_2, &coarse_expected_norm);

		VecGetArray(coarse_local_expected, &ce_vec);
		for (auto data : comm->getFineDomains()) {
			Domain &d = *data.d;
			for (int i = 0; i < n * n * n; i++) {
				ce_vec[data.local_index * n * n * n + i] = d.parent_id + i;
			}
		}
		VecRestoreArray(coarse_local_expected, &ce_vec);
		double coarse_local_expected_norm;
		VecNorm(coarse_local_expected, NORM_2, &coarse_local_expected_norm);

		// check that forward scatter works as expected
		PW<VecScatter> scatter             = comm->getScatter();
		PW<Vec>        coarse_local_result = comm->getNewCoarseDistVec();
		VecScatterBegin(scatter, coarse_expected, coarse_local_result, ADD_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter, coarse_expected, coarse_local_result, ADD_VALUES, SCATTER_FORWARD);

		PetscBool vec_scatter_forward_works;
		VecEqual(coarse_local_result, coarse_local_expected, &vec_scatter_forward_works);
		REQUIRE(vec_scatter_forward_works == PETSC_TRUE);

		double coarse_local_result_norm;
		VecNorm(coarse_local_result, NORM_2, &coarse_local_result_norm);
		REQUIRE(coarse_local_result_norm == coarse_local_expected_norm);

		// check that reverse scatter works as expected
		PW<Vec> coarse_result = coarse->getNewDomainVec();
		VecScatterBegin(scatter, coarse_local_expected, coarse_result, INSERT_VALUES,
		                SCATTER_REVERSE);
		VecScatterEnd(scatter, coarse_local_expected, coarse_result, INSERT_VALUES,
		              SCATTER_REVERSE);

		PetscBool vec_scatter_reverse_works;
		VecEqual(coarse_result, coarse_expected, &vec_scatter_reverse_works);
		REQUIRE(vec_scatter_forward_works == PETSC_TRUE);

		double coarse_result_norm;
		VecNorm(coarse_result, NORM_2, &coarse_result_norm);
		REQUIRE(coarse_result_norm == coarse_expected_norm);
	}
}
TEST_CASE("GMGAvgRstr works", "[GMG]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<InterLevelComm> comm(new InterLevelComm(coarse, fine));
		shared_ptr<GMGRestrictor>  op(new GMGAvgRstr(coarse, fine, comm));
		PW<Vec>                    coarse_expected = coarse->getNewDomainVec();
		PW<Vec>                    fine_start      = fine->getNewDomainVec();
		double *                   ce_vec;
		VecGetArray(coarse_expected, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int oct = 0; oct < 8; oct++) {
				octFill(ce_vec + d.id_local * n * n * n, oct, d.child_id[oct]);
			}
		}
		VecRestoreArray(coarse_expected, &ce_vec);
		double coarse_expected_norm;
		VecNorm(coarse_expected, NORM_2, &coarse_expected_norm);

		VecGetArray(fine_start, &ce_vec);
		for (auto data : fine->domains) {
			Domain &d = *data.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n]
						= d.id + xi / 2 + yi / 2 * n + zi / 2 * n * n;
					}
				}
			}
		}
		VecRestoreArray(fine_start, &ce_vec);

		// check that restrictor works
		PW<Vec> coarse_result = coarse->getNewDomainVec();
		op->restrict(coarse_result, fine_start);

		PetscBool restrictor_works;
		VecEqual(coarse_result, coarse_expected, &restrictor_works);
		REQUIRE(restrictor_works == PETSC_TRUE);

		double coarse_result_norm;
		VecNorm(coarse_result, NORM_2, &coarse_result_norm);
		REQUIRE(coarse_result_norm == coarse_expected_norm);
		REQUIRE(coarse_result_norm != 0);
	}
}
TEST_CASE("GMGDrctIntp works", "[GMG]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<InterLevelComm>  comm(new InterLevelComm(coarse, fine));
		shared_ptr<GMGInterpolator> op(new GMGDrctIntp(coarse, fine, comm));
		PW<Vec>                     coarse_start  = coarse->getNewDomainVec();
		PW<Vec>                     fine_expected = fine->getNewDomainVec();
		double *                    ce_vec;
		VecGetArray(coarse_start, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int oct = 0; oct < 8; oct++) {
				octFill(ce_vec + d.id_local * n * n * n, oct, d.child_id[oct]);
			}
		}
		VecRestoreArray(coarse_start, &ce_vec);

		VecGetArray(fine_expected, &ce_vec);
		for (auto data : fine->domains) {
			Domain &d = *data.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n]
						= d.id + xi / 2 + yi / 2 * n + zi / 2 * n * n;
					}
				}
			}
		}
		VecRestoreArray(fine_expected, &ce_vec);
		double fine_expected_norm;
		VecNorm(fine_expected, NORM_2, &fine_expected_norm);

		// check that interpolator works
		PW<Vec> fine_result = fine->getNewDomainVec();
		op->interpolate(coarse_start, fine_result);

		PetscBool restrictor_works;
		VecEqual(fine_result, fine_expected, &restrictor_works);
		REQUIRE(restrictor_works == PETSC_TRUE);

		double fine_result_norm;
		VecNorm(fine_result, NORM_2, &fine_result_norm);
		REQUIRE(fine_result_norm == fine_expected_norm);
		REQUIRE(fine_result_norm != 0);
	}
	PetscFinalize();
}
