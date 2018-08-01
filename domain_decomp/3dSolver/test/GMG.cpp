#include "../GMG/AvgRstr.h"
#include "../GMG/DrctIntp.h"
#include "../GMG/InterLevelComm.h"
#include "../GMG/TriLinIntp.h"
#include "catch.hpp"
#ifdef HAVE_VTK
#include "../Writers/VtkWriter.h"
#endif
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
void generateTwoLevelRefined(shared_ptr<DomainCollection> &coarse,
                             shared_ptr<DomainCollection> &fine)
{
	// generate simple tree
	OctTree t("2refine.bin");
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
	SECTION("uniform mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		PW<Vec>                         coarse_expected       = coarse->getNewDomainVec();
		PW<Vec>                         coarse_local_expected = comm->getNewCoarseDistVec();
		double *                        ce_vec;
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
	SECTION("refined mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevelRefined(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		PW<Vec>                         coarse_expected       = coarse->getNewDomainVec();
		PW<Vec>                         coarse_local_expected = comm->getNewCoarseDistVec();
		double *                        ce_vec;
		VecGetArray(coarse_expected, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int i = 0; i < n * n * n; i++) {
				ce_vec[d.id_local * n * n * n + i] = d.id + i;
			}
		}
		cerr << "SIZES: " << coarse->domains.size() << ", " << comm->getFineDomains().size()
		     << endl;
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

#ifdef HAVE_VTK
		// print?
		VtkWriter writer(*coarse, "refined_inter_level_scatter");
		writer.add(coarse_local_expected, "Local Expected");
		writer.add(coarse_expected, "Expected");
		writer.write();
#endif

		// check that forward scatter works as expected
		PW<VecScatter> scatter             = comm->getScatter();
		PW<Vec>        coarse_local_result = comm->getNewCoarseDistVec();
		VecScatterBegin(scatter, coarse_expected, coarse_local_result, ADD_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter, coarse_expected, coarse_local_result, ADD_VALUES, SCATTER_FORWARD);

		PetscBool vec_scatter_forward_works;
		VecEqual(coarse_local_result, coarse_local_expected, &vec_scatter_forward_works);
		CHECK(vec_scatter_forward_works == PETSC_TRUE);

#ifdef HAVE_VTK
		// print?
		writer.add(coarse_local_result, "Local Result");
		writer.write();
#endif

		double coarse_local_result_norm;
		VecNorm(coarse_local_result, NORM_2, &coarse_local_result_norm);
		CHECK(coarse_local_result_norm == coarse_local_expected_norm);

		// check that reverse scatter works as expected
		PW<Vec> coarse_result = coarse->getNewDomainVec();
		VecScatterBegin(scatter, coarse_local_expected, coarse_result, INSERT_VALUES,
		                SCATTER_REVERSE);
		VecScatterEnd(scatter, coarse_local_expected, coarse_result, INSERT_VALUES,
		              SCATTER_REVERSE);

#ifdef HAVE_VTK
		// print?
		writer.add(coarse_result, "Result");
		writer.write();
#endif

		PetscBool vec_scatter_reverse_works;
		VecEqual(coarse_result, coarse_expected, &vec_scatter_reverse_works);
		CHECK(vec_scatter_forward_works == PETSC_TRUE);

		double coarse_result_norm;
		VecNorm(coarse_result, NORM_2, &coarse_result_norm);
		CHECK(coarse_result_norm == coarse_expected_norm);
	}
}
TEST_CASE("GMGAvgRstr works", "[GMG]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	SECTION("uniform mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		shared_ptr<GMG::Restrictor>     op(new GMG::AvgRstr(coarse, fine, comm));
		PW<Vec>                         coarse_expected = coarse->getNewDomainVec();
		PW<Vec>                         fine_start      = fine->getNewDomainVec();
		double *                        ce_vec;
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
	SECTION("refined mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevelRefined(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		shared_ptr<GMG::Restrictor>     op(new GMG::AvgRstr(coarse, fine, comm));
		PW<Vec>                         coarse_expected = coarse->getNewDomainVec();
		PW<Vec>                         fine_start      = fine->getNewDomainVec();
		double *                        ce_vec;
		VecGetArray(coarse_expected, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			if (d.hasChildren()) {
				for (int oct = 0; oct < 8; oct++) {
					octFill(ce_vec + d.id_local * n * n * n, oct, d.child_id[oct]);
				}
			} else {
				for (int zi = 0; zi < n; zi++) {
					for (int yi = 0; yi < n; yi++) {
						for (int xi = 0; xi < n; xi++) {
							ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n]
							= d.id + xi + yi * n + zi * n * n;
						}
					}
				}
			}
		}
		VecRestoreArray(coarse_expected, &ce_vec);
		double coarse_expected_norm;
		VecNorm(coarse_expected, NORM_2, &coarse_expected_norm);

		VecGetArray(fine_start, &ce_vec);
		for (auto data : fine->domains) {
			Domain &d = *data.second;
			if (d.id != d.parent_id) {
				for (int zi = 0; zi < n; zi++) {
					for (int yi = 0; yi < n; yi++) {
						for (int xi = 0; xi < n; xi++) {
							ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n]
							= d.id + xi / 2 + yi / 2 * n + zi / 2 * n * n;
						}
					}
				}
			} else {
				for (int zi = 0; zi < n; zi++) {
					for (int yi = 0; yi < n; yi++) {
						for (int xi = 0; xi < n; xi++) {
							ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n]
							= d.id + xi + yi * n + zi * n * n;
						}
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

	SECTION("uniform mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		shared_ptr<GMG::Interpolator>   op(new GMG::DrctIntp(coarse, fine, comm));
		PW<Vec>                         coarse_start  = coarse->getNewDomainVec();
		PW<Vec>                         fine_expected = fine->getNewDomainVec();
		double *                        ce_vec;
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
}
void getXYZ(const Domain &d, const int &xi, const int &yi, const int &zi, double &x, double &y,
            double &z)
{
	const int &n   = d.n;
	double     h_x = d.x_length / n;
	double     h_y = d.y_length / n;
	double     h_z = d.z_length / n;
	if (xi == -1) {
		x = d.x_start;
	} else if (xi == n) {
		x = d.x_start + d.x_length;
	} else {
		x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
	}
	if (yi == -1) {
		y = d.y_start;
	} else if (yi == n) {
		y = d.y_start + d.y_length;
	} else {
		y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
	}
	if (zi == -1) {
		z = d.z_start;
	} else if (zi == n) {
		z = d.z_start + d.z_length;
	} else {
		z = d.z_start + h_z / 2.0 + d.z_length * zi / n;
	}
}
TEST_CASE("GMGTriLinIntp works", "[GMG]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	SECTION("uniform mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevel(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		shared_ptr<GMG::Interpolator>   op(new GMG::TriLinIntp(coarse, fine, comm));
		PW<Vec>                         coarse_start  = coarse->getNewDomainVec();
		PW<Vec>                         fine_expected = fine->getNewDomainVec();
		double *                        ce_vec;
		VecGetArray(coarse_start, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						double x, y, z;
						getXYZ(d, xi, yi, zi, x, y, z);
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n] = x + 0.5 * y - z;
					}
				}
			}
		}
		VecRestoreArray(coarse_start, &ce_vec);

		VecGetArray(fine_expected, &ce_vec);
		for (auto data : fine->domains) {
			Domain &d = *data.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						double x, y, z;
						getXYZ(d, xi, yi, zi, x, y, z);
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n] = x + 0.5 * y - z;
					}
				}
			}
		}
		VecRestoreArray(fine_expected, &ce_vec);

		// check that interpolator works
		PW<Vec> fine_result = fine->getNewDomainVec();
		op->interpolate(coarse_start, fine_result);

#ifdef HAVE_VTK
		// print?
		VtkWriter writer(*fine, "fine_domain");
		writer.add(fine_expected, "Expected");
		writer.add(fine_result, "Result");
		writer.write();
#endif

		double expected_sum;
		VecSum(fine_expected, &expected_sum);
		VecAYPX(fine_result, -1, fine_expected);

		double two_norm;
		VecNorm(fine_result, NORM_2, &two_norm);
		two_norm /= expected_sum;
		REQUIRE(two_norm < 1e-10);

		double inf_norm;
		VecNorm(fine_result, NORM_INFINITY, &inf_norm);
		inf_norm /= expected_sum;
		REQUIRE(inf_norm < 1e-10);
	}
	SECTION("refined mesh")
	{
		shared_ptr<DomainCollection> coarse;
		shared_ptr<DomainCollection> fine;
		generateTwoLevelRefined(coarse, fine);
		shared_ptr<GMG::InterLevelComm> comm(new GMG::InterLevelComm(coarse, fine));
		shared_ptr<GMG::Interpolator>   op(new GMG::TriLinIntp(coarse, fine, comm));
		PW<Vec>                         coarse_start  = coarse->getNewDomainVec();
		PW<Vec>                         fine_expected = fine->getNewDomainVec();
		double *                        ce_vec;
		VecGetArray(coarse_start, &ce_vec);
		for (auto p : coarse->domains) {
			Domain &d = *p.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						double x, y, z;
						getXYZ(d, xi, yi, zi, x, y, z);
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n] = x + 0.5 * y - z;
					}
				}
			}
		}
		VecRestoreArray(coarse_start, &ce_vec);

		VecGetArray(fine_expected, &ce_vec);
		for (auto data : fine->domains) {
			Domain &d = *data.second;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						double x, y, z;
						getXYZ(d, xi, yi, zi, x, y, z);
						ce_vec[d.id_local * n * n * n + xi + yi * n + zi * n * n] = x + 0.5 * y - z;
					}
				}
			}
		}
		VecRestoreArray(fine_expected, &ce_vec);

		// check that interpolator works
		PW<Vec> fine_result = fine->getNewDomainVec();
		op->interpolate(coarse_start, fine_result);

#ifdef HAVE_VTK
		// print?
		VtkWriter writer(*fine, "fine_domain_refined");
		writer.add(fine_expected, "Expected");
		writer.add(fine_result, "Result");
		writer.write();
#endif

		double expected_sum;
		VecSum(fine_expected, &expected_sum);
		VecAYPX(fine_result, -1, fine_expected);

		double two_norm;
		VecNorm(fine_result, NORM_2, &two_norm);
		two_norm /= expected_sum;
		CHECK(two_norm < 1e-10);

		double inf_norm;
		VecNorm(fine_result, NORM_INFINITY, &inf_norm);
		inf_norm /= expected_sum;
		CHECK(inf_norm < 1e-10);
	}
}
