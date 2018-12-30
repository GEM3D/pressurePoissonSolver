#include "TriLinIntp.h"
#include <algorithm>
#include <functional>
using namespace std;
using namespace GMG;
/*
 * @brief Get an array of all pairs of sides that touch i.e. form a corner.
 *
 * @return The array.Cores: 1
 */
inline std::array<std::array<Side<3>, 2>, 12> getPairValues()
{
	return std::array<std::array<Side<3>, 2>, 12>({{{{Side<3>::west, Side<3>::south}},
	                                                {{Side<3>::west, Side<3>::north}},
	                                                {{Side<3>::west, Side<3>::bottom}},
	                                                {{Side<3>::west, Side<3>::top}},
	                                                {{Side<3>::east, Side<3>::south}},
	                                                {{Side<3>::east, Side<3>::north}},
	                                                {{Side<3>::east, Side<3>::bottom}},
	                                                {{Side<3>::east, Side<3>::top}},
	                                                {{Side<3>::south, Side<3>::bottom}},
	                                                {{Side<3>::south, Side<3>::top}},
	                                                {{Side<3>::north, Side<3>::bottom}},
	                                                {{Side<3>::north, Side<3>::top}}}});
}
TriLinIntp::TriLinIntp(shared_ptr<DomainCollection<3>> coarse_dc,
                       shared_ptr<DomainCollection<3>> fine_dc, shared_ptr<InterLevelComm<3>> ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
struct OctInfo {
	Orthant<3> oct;
	int        n;
	int        x_start;
	int        y_start;
	int        z_start;
	OctInfo(Orthant<3> oct, int n)
	{
		this->oct = oct;
		this->n   = n;
		x_start   = oct.isOnSide(Side<3>::west) ? 0 : n / 2;
		y_start   = oct.isOnSide(Side<3>::south) ? 0 : n / 2;
		z_start   = oct.isOnSide(Side<3>::bottom) ? 0 : n / 2;
	}
	inline void getCube(const double *coarse_vec, const int coarse_idx, const int xi, const int yi,
	                    const int zi, double *cube)
	{
		cube[0]
		= coarse_vec[coarse_idx + (x_start + xi) + (y_start + yi) * n + (z_start + zi) * n * n];
		cube[1]
		= coarse_vec[coarse_idx + (x_start + xi + 1) + (y_start + yi) * n + (z_start + zi) * n * n];
		cube[2]
		= coarse_vec[coarse_idx + (x_start + xi) + (y_start + yi + 1) * n + (z_start + zi) * n * n];
		cube[3] = coarse_vec[coarse_idx + (x_start + xi + 1) + (y_start + yi + 1) * n
		                     + (z_start + zi) * n * n];
		cube[4]
		= coarse_vec[coarse_idx + (x_start + xi) + (y_start + yi) * n + (z_start + zi + 1) * n * n];
		cube[5] = coarse_vec[coarse_idx + (x_start + xi + 1) + (y_start + yi) * n
		                     + (z_start + zi + 1) * n * n];
		cube[6] = coarse_vec[coarse_idx + (x_start + xi) + (y_start + yi + 1) * n
		                     + (z_start + zi + 1) * n * n];
		cube[7] = coarse_vec[coarse_idx + (x_start + xi + 1) + (y_start + yi + 1) * n
		                     + (z_start + zi + 1) * n * n];
	}
	inline void getFineCube(double *fine_vec, const int fine_idx, const int xi, const int yi,
	                        const int zi, double **cube)
	{
		cube[0] = fine_vec + fine_idx + (1 + 2 * xi) + (1 + 2 * yi) * n + (1 + 2 * zi) * n * n;
		cube[1] = fine_vec + fine_idx + (1 + 2 * xi + 1) + (1 + 2 * yi) * n + (1 + 2 * zi) * n * n;
		cube[2] = fine_vec + fine_idx + (1 + 2 * xi) + (1 + 2 * yi + 1) * n + (1 + 2 * zi) * n * n;
		cube[3]
		= fine_vec + fine_idx + (1 + 2 * xi + 1) + (1 + 2 * yi + 1) * n + (1 + 2 * zi) * n * n;
		cube[4] = fine_vec + fine_idx + (1 + 2 * xi) + (1 + 2 * yi) * n + (1 + 2 * zi + 1) * n * n;
		cube[5]
		= fine_vec + fine_idx + (1 + 2 * xi + 1) + (1 + 2 * yi) * n + (1 + 2 * zi + 1) * n * n;
		cube[6]
		= fine_vec + fine_idx + (1 + 2 * xi) + (1 + 2 * yi + 1) * n + (1 + 2 * zi + 1) * n * n;
		cube[7]
		= fine_vec + fine_idx + (1 + 2 * xi + 1) + (1 + 2 * yi + 1) * n + (1 + 2 * zi + 1) * n * n;
	}
};
class Helper
{
	public:
	virtual void apply(double *u_fine, double *u_coarse);
};
class InteriorHelper : public Helper
{
	private:
	int                  n;
	int                  x_start;
	int                  y_start;
	int                  z_start;
	static constexpr int center_coeffs[8][8]
	= {{27, 9, 9, 3, 9, 3, 3, 1}, {9, 27, 3, 9, 3, 9, 1, 3}, {9, 3, 27, 9, 3, 1, 9, 3},
	   {3, 9, 9, 27, 1, 3, 3, 9}, {9, 3, 3, 1, 27, 9, 9, 3}, {3, 9, 1, 3, 9, 27, 3, 9},
	   {3, 1, 9, 3, 9, 3, 27, 9}, {1, 3, 3, 9, 3, 9, 9, 27}};

	public:
	InteriorHelper(OctInfo info)
	{
		n       = info.n;
		x_start = info.oct.isOnSide(Side<3>::west) ? 0 : n / 2;
		y_start = info.oct.isOnSide(Side<3>::south) ? 0 : n / 2;
		z_start = info.oct.isOnSide(Side<3>::bottom) ? 0 : n / 2;
	}
	void apply(double *u_fine, double *u_coarse)
	{
		for (int zi = 0; zi < n / 2 - 1; zi++) {
			for (int yi = 0; yi < n / 2 - 1; yi++) {
				for (int xi = 0; xi < n / 2 - 1; xi++) {
					double cube[8];
					cube[0]
					= u_coarse[(x_start + xi) + (y_start + yi) * n + (z_start + zi) * n * n];
					cube[1]
					= u_coarse[(x_start + xi + 1) + (y_start + yi) * n + (z_start + zi) * n * n];
					cube[2]
					= u_coarse[(x_start + xi) + (y_start + yi + 1) * n + (z_start + zi) * n * n];
					cube[3] = u_coarse[(x_start + xi + 1) + (y_start + yi + 1) * n
					                   + (z_start + zi) * n * n];
					cube[4]
					= u_coarse[(x_start + xi) + (y_start + yi) * n + (z_start + zi + 1) * n * n];
					cube[5] = u_coarse[(x_start + xi + 1) + (y_start + yi) * n
					                   + (z_start + zi + 1) * n * n];
					cube[6] = u_coarse[(x_start + xi) + (y_start + yi + 1) * n
					                   + (z_start + zi + 1) * n * n];
					cube[7] = u_coarse[(x_start + xi + 1) + (y_start + yi + 1) * n
					                   + (z_start + zi + 1) * n * n];

					double *fine_cube[8];
					fine_cube[0] = u_fine + (1 + 2 * xi) + (1 + 2 * yi) * n + (1 + 2 * zi) * n * n;
					fine_cube[1]
					= u_fine + (1 + 2 * xi + 1) + (1 + 2 * yi) * n + (1 + 2 * zi) * n * n;
					fine_cube[2]
					= u_fine + (1 + 2 * xi) + (1 + 2 * yi + 1) * n + (1 + 2 * zi) * n * n;
					fine_cube[3]
					= u_fine + (1 + 2 * xi + 1) + (1 + 2 * yi + 1) * n + (1 + 2 * zi) * n * n;
					fine_cube[4]
					= u_fine + (1 + 2 * xi) + (1 + 2 * yi) * n + (1 + 2 * zi + 1) * n * n;
					fine_cube[5]
					= u_fine + (1 + 2 * xi + 1) + (1 + 2 * yi) * n + (1 + 2 * zi + 1) * n * n;
					fine_cube[6]
					= u_fine + (1 + 2 * xi) + (1 + 2 * yi + 1) * n + (1 + 2 * zi + 1) * n * n;
					fine_cube[7]
					= u_fine + (1 + 2 * xi + 1) + (1 + 2 * yi + 1) * n + (1 + 2 * zi + 1) * n * n;
					for (int f = 0; f < 8; f++) {
						double result = 0;
						for (int c = 0; c < 8; c++) {
							result += cube[c] * center_coeffs[f][c];
						}
						*fine_cube[f] += result / 64;
					}
				}
			}
		}
	}
};
class ExtFaceHelper : public Helper
{
	private:
	static constexpr int coeffs[4][8] = {{45, 15, 15, 5, -9, -3, -3, -1},
	                                     {15, 45, 5, 15, -3, -9, -1, -3},
	                                     {15, 5, 45, 15, -3, -1, -9, -3},
	                                     {5, 15, 15, 45, -1, -3, -3, -9}};
	static constexpr int divisor      = 64;
	int                  start        = -1;
	int                  start_inner  = -1;
	int                  start_coarse = -1;
	int                  stride_x     = -1;
	int                  stride_y     = -1;
	int                  n;

	public:
	ExtFaceHelper(OctInfo info, Side<3> s)
	{
		n = info.n;
		// set strides and starting indexes
		if (s == Side<3>::west || s == Side<3>::east) {
			stride_x = n;
			stride_y = n * n;
			if (s == Side<3>::west) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_inner  = info.x_start + 1 + info.y_start * n + info.z_start * n * n;
				start_coarse = 0;
			} else {
				start = (info.x_start + n / 2 - 1) + info.y_start * n + info.z_start * n * n;
				start_inner
				= (info.x_start + n / 2 - 1 - 1) + info.y_start * n + info.z_start * n * n;
				start_coarse = n - 1;
			}
		}
		if (s == Side<3>::south || s == Side<3>::north) {
			stride_x = 1;
			stride_y = n * n;
			if (s == Side<3>::south) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_inner  = info.x_start + (info.y_start + 1) * n + info.z_start * n * n;
				start_coarse = 0;
			} else {
				start = info.x_start + (info.y_start + n / 2 - 1) * n + info.z_start * n * n;
				start_inner
				= info.x_start + (info.y_start + n / 2 - 1 - 1) * n + info.z_start * n * n;
				start_coarse = (n - 1) * n;
			}
		}
		if (s == Side<3>::bottom || s == Side<3>::top) {
			stride_x = 1;
			stride_y = n;
			if (s == Side<3>::bottom) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_inner  = info.x_start + info.y_start * n + (info.z_start + 1) * n * n;
				start_coarse = 0;
			} else {
				start = info.x_start + info.y_start * n + (info.z_start + n / 2 - 1) * n * n;
				start_inner
				= info.x_start + info.y_start * n + (info.z_start + n / 2 - 1 - 1) * n * n;
				start_coarse = (n - 1) * n * n;
			}
		}
	}
	void apply(double *u_fine, double *u_coarse)
	{
		for (int yi = 0; yi < n / 2 - 1; yi++) {
			for (int xi = 0; xi < n / 2 - 1; xi++) {
				double cube[8];
				cube[0] = u_coarse[start + xi * stride_x + yi * stride_y];
				cube[1] = u_coarse[start + (xi + 1) * stride_x + yi * stride_y];
				cube[2] = u_coarse[start + xi * stride_x + (yi + 1) * stride_y];
				cube[3] = u_coarse[start + (xi + 1) * stride_x + (yi + 1) * stride_y];
				cube[4] = u_coarse[start_inner + xi * stride_x + yi * stride_y];
				cube[5] = u_coarse[start_inner + (xi + 1) * stride_x + yi * stride_y];
				cube[6] = u_coarse[start_inner + xi * stride_x + (yi + 1) * stride_y];
				cube[7] = u_coarse[start_inner + (xi + 1) * stride_x + (yi + 1) * stride_y];

				double fine[4] = {0, 0, 0, 0};
				for (int j = 0; j < 4; j++) {
					for (int i = 0; i < 8; i++) {
						fine[j] += cube[i] * coeffs[j][i];
					}
					fine[j] /= divisor;
				}

				u_fine[start_coarse + (1 + xi * 2) * stride_x + (1 + yi * 2) * stride_y] += fine[0];
				u_fine[start_coarse + (1 + xi * 2 + 1) * stride_x + (1 + yi * 2) * stride_y]
				+= fine[1];
				u_fine[start_coarse + (1 + xi * 2) * stride_x + (1 + yi * 2 + 1) * stride_y]
				+= fine[2];
				u_fine[start_coarse + (1 + xi * 2 + 1) * stride_x + (1 + yi * 2 + 1) * stride_y]
				+= fine[3];
			}
		}
	}
};
class IntFaceHelper : public Helper
{
	private:
	static constexpr int coeffs[4][8] = {{27, 9, 9, 3, 9, 3, 3, 1},
	                                     {9, 27, 3, 9, 3, 9, 1, 3},
	                                     {9, 3, 27, 9, 3, 1, 9, 3},
	                                     {3, 9, 9, 27, 1, 3, 3, 9}};
	static constexpr int divisor      = 64;
	int                  start        = -1;
	int                  start_outer  = -1;
	int                  start_coarse = -1;
	int                  stride_x     = -1;
	int                  stride_y     = -1;
	int                  n;

	public:
	IntFaceHelper(OctInfo info, Side<3> s)
	{
		n = info.n;
		// set strides and starting indexes
		if (s == Side<3>::west || s == Side<3>::east) {
			stride_x = n;
			stride_y = n * n;
			if (s == Side<3>::west) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_outer  = info.x_start - 1 + info.y_start * n + info.z_start * n * n;
				start_coarse = 0;
			} else {
				start = (info.x_start + n / 2 - 1) + info.y_start * n + info.z_start * n * n;
				start_outer
				= (info.x_start + n / 2 - 1 + 1) + info.y_start * n + info.z_start * n * n;
				start_coarse = n - 1;
			}
		}
		if (s == Side<3>::south || s == Side<3>::north) {
			stride_x = 1;
			stride_y = n * n;
			if (s == Side<3>::south) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_outer  = info.x_start + (info.y_start - 1) * n + info.z_start * n * n;
				start_coarse = 0;
			} else {
				start = info.x_start + (info.y_start + n / 2 - 1) * n + info.z_start * n * n;
				start_outer
				= info.x_start + (info.y_start + n / 2 - 1 + 1) * n + info.z_start * n * n;
				start_coarse = (n - 1) * n;
			}
		}
		if (s == Side<3>::bottom || s == Side<3>::top) {
			stride_x = 1;
			stride_y = n;
			if (s == Side<3>::bottom) {
				start        = info.x_start + info.y_start * n + info.z_start * n * n;
				start_outer  = info.x_start + info.y_start * n + (info.z_start - 1) * n * n;
				start_coarse = 0;
			} else {
				start = info.x_start + info.y_start * n + (info.z_start + n / 2 - 1) * n * n;
				start_outer
				= info.x_start + info.y_start * n + (info.z_start + n / 2 - 1 + 1) * n * n;
				start_coarse = (n - 1) * n * n;
			}
		}
	}
	void apply(double *u_fine, double *u_coarse)
	{
		for (int yi = 0; yi < n / 2 - 1; yi++) {
			for (int xi = 0; xi < n / 2 - 1; xi++) {
				double cube[8];
				cube[0] = u_coarse[start + xi * stride_x + yi * stride_y];
				cube[1] = u_coarse[start + (xi + 1) * stride_x + yi * stride_y];
				cube[2] = u_coarse[start + xi * stride_x + (yi + 1) * stride_y];
				cube[3] = u_coarse[start + (xi + 1) * stride_x + (yi + 1) * stride_y];
				cube[4] = u_coarse[start_outer + xi * stride_x + yi * stride_y];
				cube[5] = u_coarse[start_outer + (xi + 1) * stride_x + yi * stride_y];
				cube[6] = u_coarse[start_outer + xi * stride_x + (yi + 1) * stride_y];
				cube[7] = u_coarse[start_outer + (xi + 1) * stride_x + (yi + 1) * stride_y];

				double fine[4] = {0, 0, 0, 0};
				for (int j = 0; j < 4; j++) {
					for (int i = 0; i < 8; i++) {
						fine[j] += cube[i] * coeffs[j][i];
					}
					fine[j] /= divisor;
				}

				u_fine[start_coarse + (1 + xi * 2) * stride_x + (1 + yi * 2) * stride_y] += fine[0];
				u_fine[start_coarse + (1 + xi * 2 + 1) * stride_x + (1 + yi * 2) * stride_y]
				+= fine[1];
				u_fine[start_coarse + (1 + xi * 2) * stride_x + (1 + yi * 2 + 1) * stride_y]
				+= fine[2];
				u_fine[start_coarse + (1 + xi * 2 + 1) * stride_x + (1 + yi * 2 + 1) * stride_y]
				+= fine[3];
			}
		}
	}
};
class EdgeHelper : public Helper
{
	private:
	int                  start[4]   = {-1, -1, -1, -1};
	int                  start_fine = -1;
	int                  stride     = -1;
	int                  n;
	int                  num_int;
	static constexpr int coeffs[3][2][8]
	= {{{75, 25, -15, -5, -15, -5, 3, 1}, {25, 75, -5, -15, -5, -15, 1, 3}},
	   {{45, 15, 15, 5, -9, -3, -3, -1}, {15, 45, 5, 15, -3, -9, -1, -3}},
	   {{27, 9, 9, 3, 9, 3, 3, 1}, {9, 27, 3, 9, 3, 9, 1, 3}}};
	static constexpr int divisor = 64;

	public:
	EdgeHelper(int n, Orthant<3> o, std::array<Side<3>, 2> sides)
	{
		this->n = n;
		num_int = 0;
		for (Side<3> s : sides) {
			if (!o.isOnSide(s)) num_int++;
		}
		sort(sides.begin(), sides.end(),
		     [=](Side<3> s1, Side<3> s2) { return o.isOnSide(s1) < o.isOnSide(s2); });
		// get axis
		Side<3> orth_side;
		int     axis = (0x1 << (sides[0].toInt() / 2)) | (0x1 << (sides[1].toInt() / 2));
		switch (axis) {
			case 0b110:
				orth_side = Side<3>::west;
				break;
			case 0b101:
				orth_side = Side<3>::south;
				break;
			case 0b011:
				orth_side = Side<3>::bottom;
				break;
		}
		auto pow = [](int x, int n) {
			int retval = 1;
			for (int i = 0; i < n; i++) {
				retval *= x;
			}
			return retval;
		};
		int strides[2];
		strides[0] = pow(n, sides[0].toInt() / 2);
		strides[1] = pow(n, sides[1].toInt() / 2);
		stride     = pow(n, orth_side.toInt() / 2);
		int fine_coord[2];
		fine_coord[0] = sides[0].isLowerOnAxis() ? 0 : (n - 1);
		fine_coord[1] = sides[1].isLowerOnAxis() ? 0 : (n - 1);
		int first_inds[2];
		if (o.isOnSide(sides[0])) {
			if (sides[0].isLowerOnAxis()) {
				first_inds[0] = 0;
				first_inds[1] = 1;
			} else {
				first_inds[0] = n - 1;
				first_inds[1] = n - 2;
			}
		} else {
			if (sides[0].isLowerOnAxis()) {
				first_inds[0] = n / 2;
				first_inds[1] = n / 2 - 1;
			} else {
				first_inds[0] = n / 2 - 1;
				first_inds[1] = n / 2;
			}
		}
		int second_inds[2];
		if (o.isOnSide(sides[1])) {
			if (sides[1].isLowerOnAxis()) {
				second_inds[0] = 0;
				second_inds[1] = 1;
			} else {
				second_inds[0] = n - 1;
				second_inds[1] = n - 2;
			}
		} else {
			if (sides[1].isLowerOnAxis()) {
				second_inds[0] = n / 2;
				second_inds[1] = n / 2 - 1;
			} else {
				second_inds[0] = n / 2 - 1;
				second_inds[1] = n / 2;
			}
		}

		start_fine = strides[0] * fine_coord[0] + strides[1] * fine_coord[1];
		start[0]   = strides[0] * first_inds[0] + strides[1] * second_inds[0];
		start[1]   = strides[0] * first_inds[1] + strides[1] * second_inds[0];
		start[2]   = strides[0] * first_inds[0] + strides[1] * second_inds[1];
		start[3]   = strides[0] * first_inds[1] + strides[1] * second_inds[1];
		if (!o.isOnSide(orth_side)) {
			for (int i = 0; i < 4; i++) {
				start[i] += stride * n / 2;
			}
		}
	}
	void apply(double *u_fine, double *u_coarse)
	{
		for (int xi = 0; xi < n / 2 - 1; xi++) {
			double cube[8];
			cube[0] = u_coarse[start[0] + xi * stride];
			cube[1] = u_coarse[start[0] + (xi + 1) * stride];
			cube[2] = u_coarse[start[1] + xi * stride];
			cube[3] = u_coarse[start[1] + (xi + 1) * stride];
			cube[4] = u_coarse[start[2] + xi * stride];
			cube[5] = u_coarse[start[2] + (xi + 1) * stride];
			cube[6] = u_coarse[start[3] + xi * stride];
			cube[7] = u_coarse[start[3] + (xi + 1) * stride];

			double fine[2] = {0, 0};
			for (int j = 0; j < 2; j++) {
				for (int i = 0; i < 8; i++) {
					fine[j] += cube[i] * coeffs[num_int][j][i];
				}
				fine[j] /= divisor;
			}

			u_fine[start_fine + (1 + xi * 2) * stride] += fine[0];
			u_fine[start_fine + (1 + xi * 2 + 1) * stride] += fine[1];
		}
	}
};
constexpr int EdgeHelper::coeffs[3][2][8];
class CornerHelper : public Helper
{
	private:
	int                  cube[8]  = {-1, -1, -1, -1, -1, -1, -1, -1};
	int                  fine_idx = -1;
	int                  stride   = -1;
	int                  n;
	int                  num_int;
	static constexpr int coeffs[4][8] = {{125, -25, -25, 5, -25, 5, 5, -1},
	                                     {75, 25, -15, -5, -15, -5, 3, 1},
	                                     {45, 15, 15, 5, -9, -3, -3, -1},
	                                     {27, 9, 9, 3, 9, 3, 3, 1}};
	static constexpr int divisor      = 64;

	public:
	CornerHelper(int n, Orthant<3> o, Orthant<3> o_corner)
	{
		this->n                 = n;
		array<Side<3>, 3> sides = o_corner.getExteriorSides();
		num_int                 = 0;
		for (Side<3> s : sides) {
			if (!o.isOnSide(s)) num_int++;
		}
		sort(sides.begin(), sides.end(),
		     [=](Side<3> s1, Side<3> s2) { return o.isOnSide(s1) < o.isOnSide(s2); });
		auto pow = [](int x, int n) {
			int retval = 1;
			for (int i = 0; i < n; i++) {
				retval *= x;
			}
			return retval;
		};

		int strides[3];
		for (int i = 0; i < 3; i++) {
			strides[i] = pow(n, sides[i].toInt() / 2);
		}

		int fine_inds[3];
		for (int i = 0; i < 3; i++) {
			fine_inds[i] = sides[i].isLowerOnAxis() ? 0 : (n - 1);
		}

		fine_idx = 0;
		for (int i = 0; i < 3; i++) {
			fine_idx += fine_inds[i] * strides[i];
		}

		int inds[3][2];
		for (int i = 0; i < 3; i++) {
			if (o.isOnSide(sides[i])) {
				if (sides[i].isLowerOnAxis()) {
					inds[i][0] = 0;
					inds[i][1] = 1;
				} else {
					inds[i][0] = n - 1;
					inds[i][1] = n - 2;
				}
			} else {
				if (sides[i].isLowerOnAxis()) {
					inds[i][0] = n / 2;
					inds[i][1] = n / 2 - 1;
				} else {
					inds[i][0] = n / 2 - 1;
					inds[i][1] = n / 2;
				}
			}
		}
		cube[0] = strides[0] * inds[0][0] + strides[1] * inds[1][0] + strides[2] * inds[2][0];
		cube[1] = strides[0] * inds[0][1] + strides[1] * inds[1][0] + strides[2] * inds[2][0];
		cube[2] = strides[0] * inds[0][0] + strides[1] * inds[1][1] + strides[2] * inds[2][0];
		cube[3] = strides[0] * inds[0][1] + strides[1] * inds[1][1] + strides[2] * inds[2][0];
		cube[4] = strides[0] * inds[0][0] + strides[1] * inds[1][0] + strides[2] * inds[2][1];
		cube[5] = strides[0] * inds[0][1] + strides[1] * inds[1][0] + strides[2] * inds[2][1];
		cube[6] = strides[0] * inds[0][0] + strides[1] * inds[1][1] + strides[2] * inds[2][1];
		cube[7] = strides[0] * inds[0][1] + strides[1] * inds[1][1] + strides[2] * inds[2][1];
		/*
		int second_inds[2];
		if (o.isOnSide(sides[1])) {
		    if (sides[1].isLowerOnAxis()) {
		        second_inds[0] = 0;
		        second_inds[1] = 1;
		    } else {
		        second_inds[0] = n - 1;
		        second_inds[1] = n - 2;
		    }
		} else {
		    second_inds[0] = n / 2 - 1;
		    second_inds[1] = n / 2;
		}

		start_fine = strides[0] * fine_coord[0] + strides[1] * fine_coord[1];
		start[0]   = strides[0] * first_inds[0] + strides[1] * second_inds[0];
		start[1]   = strides[0] * first_inds[1] + strides[1] * second_inds[0];
		start[2]   = strides[0] * first_inds[0] + strides[1] * second_inds[1];
		start[3]   = strides[0] * first_inds[1] + strides[1] * second_inds[1];
		if (!o.isOnSide(orth_side)) {
		    for (int i = 0; i < 4; i++) {
		        start[i] += stride * n / 2;
		    }
		}
		*/
	}
	void apply(double *u_fine, double *u_coarse)
	{
		double fine = 0;
		for (int i = 0; i < 8; i++) {
			fine += u_coarse[cube[i]] * coeffs[num_int][i];
		}
		fine /= divisor;

		u_fine[fine_idx] += fine;
	}
};
constexpr int CornerHelper::coeffs[4][8];

void TriLinIntp::interpolate(PW<Vec> coarse, PW<Vec> fine) const
{
	// get vectors
	double *u_fine;
	double *u_coarse;
	PW<Vec> coarse_tmp = ilc->getNewCoarseDistVec();
	// scatter
	PW<VecScatter> scatter = ilc->getScatter();
	VecScatterBegin(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);

	VecGetArray(fine, &u_fine);
	VecGetArray(coarse_tmp, &u_coarse);
	for (auto p : ilc->getFineDomains()) {
		Domain<3> &d          = *p.d;
		int        n          = d.n;
		int        coarse_idx = p.local_index * n * n * n;

		Orthant<3> oct      = d.oct_on_parent;
		int        fine_idx = d.id_local * n * n * n;

		if (d.id == d.parent_id) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + xi + yi * n + zi * n * n];
					}
				}
			}
		} else {
			OctInfo info(oct, n);
			// interior points
			{
				InteriorHelper helper(info);
				helper.apply(u_fine + fine_idx, u_coarse + coarse_idx);
			}
			// faces
			for (Side<3> s : Side<3>::getValues()) {
				if (oct.isOnSide(s)) {
					ExtFaceHelper helper(info, s);
					helper.apply(u_fine + fine_idx, u_coarse + coarse_idx);
				} else {
					IntFaceHelper helper(info, s);
					helper.apply(u_fine + fine_idx, u_coarse + coarse_idx);
				}
			}
			// edges
			for (std::array<Side<3>, 2> sides : getPairValues()) {
				EdgeHelper helper(n, oct, sides);
				helper.apply(u_fine + fine_idx, u_coarse + coarse_idx);
			}
			for (Orthant<3> o : Orthant<3>::getValues()) {
				CornerHelper helper(n, oct, o);
				helper.apply(u_fine + fine_idx, u_coarse + coarse_idx);
			}
		}
	}
	VecRestoreArray(fine, &u_fine);
	VecRestoreArray(coarse_tmp, &u_coarse);
}
