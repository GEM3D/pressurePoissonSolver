#ifndef IFACE_H
#define IFACE_H
enum class BlockType {
	plain,
	fine,
	fine_out_left,
	fine_out_right,
	coarse,
	coarse_out_left,
	coarse_out_right
};
class ColIface
{
	public:
	const static int left_offset            = 4;
	const static int side_size              = 6 + 2;
	const static int right_offset           = 28 + 2 * 4;
	const static int domain_boundary_offset = 46 + 2 * 7;
	const static int zp_offset              = 46 + 2 * 7 + 2;
	const static int size                   = 46 + 2 * 7 + 2 + 2;
	bool             right;
	int              axis;
	int              refine_level = 1;
	int              main_i;
	Side             s;
	std::array<int, 4> global_i      = {{-1, -1, -1, -1}};
	std::array<int, 4> center_i      = {{-1, -1, -1, -1}};
	std::array<int, 4> refined_left  = {{-1, -1, -1, -1}};
	std::array<int, 4> refined_right = {{-1, -1, -1, -1}};

	std::array<std::bitset<4>, 9> domain_boundaries;

	bool           zero_patch = false;
	std::bitset<4> neumann;
	std::bitset<4> hasCoarseNbr;
	std::bitset<4> hasFineNbr;
	std::bitset<4> isCoarseLeft;
	void           setNeumann()
	{
		for (int q = 0; q < 4; q++) {
			neumann[q] = global_i[q] == -1;
		}
	}
	static void readIfaces(std::set<ColIface> &ifaces, int_vector_type &iface_info)
	{
		auto iface_view = iface_info.getLocalView<Kokkos::HostSpace>();
		for (size_t i = 0; i < iface_view.dimension(0); i += size) {
			ColIface left;
			ColIface right;

			left.zero_patch  = iface_view(i + zp_offset, 0);
			right.zero_patch = iface_view(i + zp_offset + 1, 0);
			int blr          = iface_view(i, 0);

			left.right            = false;
			left.axis             = iface_view(i + 1, 0);
			left.refine_level     = iface_view(i + 2, 0);
			left.global_i[0]      = iface_view(i + left_offset, 0);
			left.center_i[0]      = iface_view(i + left_offset + 1, 0);
			left.hasFineNbr[0]    = iface_view(i + left_offset + 2, 0);
			left.hasCoarseNbr[0]  = iface_view(i + left_offset + 3, 0);
			left.isCoarseLeft[0]  = iface_view(i + left_offset + 4, 0);
			left.refined_left[0]  = iface_view(i + left_offset + 5, 0);
			left.refined_right[0] = iface_view(i + left_offset + 6, 0);

			for (int q = 1; q < 4; q++) {
				left.global_i[q]      = iface_view(i + left_offset + q * side_size, 0);
				left.center_i[q]      = iface_view(i + left_offset + q * side_size + 1, 0);
				left.hasFineNbr[q]    = iface_view(i + left_offset + q * side_size + 2, 0);
				left.hasCoarseNbr[q]  = iface_view(i + left_offset + q * side_size + 3, 0);
				left.isCoarseLeft[q]  = iface_view(i + left_offset + q * side_size + 4, 0);
				left.refined_left[q]  = iface_view(i + left_offset + q * side_size + 5, 0);
				left.refined_right[q] = iface_view(i + left_offset + q * side_size + 6, 0);
			}

			right.right            = true;
			right.axis             = iface_view(i + 1, 0);
			right.refine_level     = iface_view(i + 2, 0);
			right.global_i[0]      = iface_view(i + left_offset, 0);
			right.center_i[0]      = iface_view(i + left_offset + 1, 0);
			right.hasFineNbr[0]    = iface_view(i + left_offset + 2, 0);
			right.hasCoarseNbr[0]  = iface_view(i + left_offset + 3, 0);
			right.isCoarseLeft[0]  = iface_view(i + left_offset + 4, 0);
			right.refined_left[0]  = iface_view(i + left_offset + 5, 0);
			right.refined_right[0] = iface_view(i + left_offset + 6, 0);
			for (int q = 1; q < 4; q++) {
				right.global_i[q]      = iface_view(i + right_offset + (q - 1) * side_size, 0);
				right.center_i[q]      = iface_view(i + right_offset + (q - 1) * side_size + 1, 0);
				right.hasFineNbr[q]    = iface_view(i + right_offset + (q - 1) * side_size + 2, 0);
				right.hasCoarseNbr[q]  = iface_view(i + right_offset + (q - 1) * side_size + 3, 0);
				right.isCoarseLeft[q]  = iface_view(i + right_offset + (q - 1) * side_size + 4, 0);
				right.refined_left[q]  = iface_view(i + right_offset + (q - 1) * side_size + 5, 0);
				right.refined_right[q] = iface_view(i + right_offset + (q - 1) * side_size + 6, 0);
			}

			if (blr == 0) {
				ifaces.insert(left);
				ifaces.insert(right);
			} else if (blr == 1) {
				ifaces.insert(left);
			} else if (blr == 2) {
				ifaces.insert(right);
			}
		}
	}
	friend bool operator<(const ColIface &l, const ColIface &r)
	{
		return std::tie(l.global_i[0], l.right) < std::tie(r.global_i[0], r.right);
	}
	friend bool operator==(const ColIface &l, const ColIface &r)
	{
		return std::tie(l.neumann, l.zero_patch) == std::tie(r.neumann, r.zero_patch);
	}
	/*friend bool operator!=(const ColIface &l, const ColIface &r)
	{
	    return std::tie(l.l_south, l.t_south, l.l_east, l.t_east, l.l_north, l.t_north, l.l_west,
	                    l.t_west)
	           == std::tie(r.l_south, r.t_south, r.l_west, r.t_west, r.l_north, r.t_north, r.l_east,
	                       r.t_east);
	}*/
	friend std::ostream &operator<<(std::ostream &os, const ColIface &iface)
	{
		os << "N:     " << iface.global_i[0] << std::endl;
		os << "E:     " << iface.global_i[1] << std::endl;
		os << "S:     " << iface.global_i[2] << std::endl;
		os << "W:     " << iface.global_i[3] << std::endl;
		os << "AXIS:  " << iface.axis << std::endl;
		os << "RIGHT: " << iface.right << std::endl;
		return os;
	}
	static void writeIfaces(Domain &d, int_vector_type &iface_info)
	{
		auto dist_view = iface_info.getLocalView<Kokkos::HostSpace>();
		// north
		Side iface_s = Side::north;
		do {
			if (d.hasNbr(iface_s)) {
				int iface_i = d.index(iface_s) * ColIface::size;
				int i       = 0;
				// left
				if (iface_s == Side::north || iface_s == Side::east) {
					if (d.hasCoarseNbr(iface_s) || d.hasFineNbr(iface_s)) {
						dist_view(iface_i, 0) = 1;
					}
					if (iface_s == Side::north || iface_s == Side::south) {
						dist_view(iface_i + 1, 0) = 0;
					} else {
						dist_view(iface_i + 1, 0) = 1;
					}
					dist_view(iface_i + 2, 0) = d.ds.refine_level;
					Side s = iface_s;
					dist_view(iface_i + zp_offset, 0) = d.zero_patch;
					i = ColIface::left_offset;
					do {
						dist_view(iface_i + i, 0)     = d.globalIndex(s);
						dist_view(iface_i + i + 1, 0) = d.globalIndexCenter(s);
						dist_view(iface_i + i + 2, 0) = d.hasFineNbr(s);
						dist_view(iface_i + i + 3, 0) = d.hasCoarseNbr(s);
						dist_view(iface_i + i + 4, 0) = d.isCoarseLeft(s);
						dist_view(iface_i + i + 5, 0) = d.globalIndexRefinedLeft(s);
						dist_view(iface_i + i + 6, 0) = d.globalIndexRefinedRight(s);
						i += ColIface::side_size;
						s++;
					} while (s != iface_s);
				}
				// right
				if (iface_s == Side::south || iface_s == Side::west) {
					if (d.hasCoarseNbr(iface_s) || d.hasFineNbr(iface_s)) {
						dist_view(iface_i, 0) = 2;
						if (iface_s == Side::north || iface_s == Side::south) {
							dist_view(iface_i + 1, 0) = 0;
						} else {
							dist_view(iface_i + 1, 0) = 1;
						}
						dist_view(iface_i + 2, 0)               = d.ds.refine_level;
						dist_view(iface_i + left_offset, 0)     = d.globalIndex(iface_s);
						dist_view(iface_i + left_offset + 1, 0) = d.globalIndexCenter(iface_s);
						dist_view(iface_i + left_offset + 2, 0) = d.hasFineNbr(iface_s);
						dist_view(iface_i + left_offset + 3, 0) = d.hasCoarseNbr(iface_s);
						dist_view(iface_i + left_offset + 4, 0) = d.isCoarseLeft(iface_s);
						dist_view(iface_i + left_offset + 5, 0) = d.globalIndexRefinedLeft(iface_s);
						dist_view(iface_i + left_offset + 6, 0)
						= d.globalIndexRefinedRight(iface_s);
					}
					Side s = iface_s;
					s++;
					dist_view(iface_i + zp_offset + 1, 0) = d.zero_patch;
					i = ColIface::right_offset;
					do {
						dist_view(iface_i + i, 0)     = d.globalIndex(s);
						dist_view(iface_i + i + 1, 0) = d.globalIndexCenter(s);
						dist_view(iface_i + i + 2, 0) = d.hasFineNbr(s);
						dist_view(iface_i + i + 3, 0) = d.hasCoarseNbr(s);
						dist_view(iface_i + i + 4, 0) = d.isCoarseLeft(s);
						dist_view(iface_i + i + 5, 0) = d.globalIndexRefinedLeft(s);
						dist_view(iface_i + i + 6, 0) = d.globalIndexRefinedRight(s);
						i += ColIface::side_size;
						s++;
					} while (s != iface_s);
				}
			}
			// add domain boundary information to array
			iface_s++;
		} while (iface_s != Side::north);
	}
};
enum class IfaceType { main, left, right };
class RowIface
{
	public:
	const static int main_offset         = 0;
	const static int side_offset         = 1;
	const static int type_offset         = 2;
	const static int hfn_offset          = 3;
	const static int hcn_offset          = 4;
	const static int icl_offset          = 5;
	const static int zp_offset           = 6;
	const static int global_offset       = 7;
	const static int right_offset        = 11;
	const static int coarse_right_offset = 22;
	const static int size                = 33;

	int       main_i;
	Side      s;
	IfaceType type;
	bool      hasFineNbr;
	bool      hasCoarseNbr;
	bool      isCoarseLeft;
	std::array<int, 4> global_i = {{-1, -1, -1, -1}};

	bool           zero_patch = false;
	std::bitset<4> neumann;
	void           setNeumann()
	{
		for (int q = 0; q < 4; q++) {
			neumann[q] = global_i[q] == -1;
		}
	}
	static void readIfaces(std::set<RowIface> &ifaces, int_vector_type &iface_info)
	{
		auto iface_view = iface_info.getLocalView<Kokkos::HostSpace>();
		for (size_t i = 0; i < iface_view.dimension(0); i += size) {
			RowIface left, right, coarse_right;

			int iface_i = i;

			left.main_i       = iface_view(iface_i + main_offset, 0);
			left.s            = static_cast<Side>(iface_view(iface_i + side_offset, 0));
			left.type         = static_cast<IfaceType>(iface_view(iface_i + type_offset, 0));
			left.hasFineNbr   = iface_view(iface_i + hfn_offset, 0);
			left.hasCoarseNbr = iface_view(iface_i + hcn_offset, 0);
			left.isCoarseLeft = iface_view(iface_i + icl_offset, 0);
			left.zero_patch   = iface_view(iface_i + zp_offset, 0);
			for (int q = 0; q < 4; q++) {
				left.global_i[q] = iface_view(iface_i + global_offset + q, 0);
			}

			iface_i = i + right_offset;

			right.main_i       = iface_view(iface_i + main_offset, 0);
			right.s            = static_cast<Side>(iface_view(iface_i + side_offset, 0));
			right.type         = static_cast<IfaceType>(iface_view(iface_i + type_offset, 0));
			right.hasFineNbr   = iface_view(iface_i + hfn_offset, 0);
			right.hasCoarseNbr = iface_view(iface_i + hcn_offset, 0);
			right.isCoarseLeft = iface_view(iface_i + icl_offset, 0);
			right.zero_patch   = iface_view(iface_i + zp_offset, 0);
			for (int q = 0; q < 4; q++) {
				right.global_i[q] = iface_view(iface_i + global_offset + q, 0);
			}

			if (iface_view(i + coarse_right_offset, 0) >= 0) {
				iface_i             = i + coarse_right_offset;
				coarse_right.main_i = iface_view(iface_i + main_offset, 0);
				coarse_right.s      = static_cast<Side>(iface_view(iface_i + side_offset, 0));
				coarse_right.type   = static_cast<IfaceType>(iface_view(iface_i + type_offset, 0));
				coarse_right.hasFineNbr   = iface_view(iface_i + hfn_offset, 0);
				coarse_right.hasCoarseNbr = iface_view(iface_i + hcn_offset, 0);
				coarse_right.isCoarseLeft = iface_view(iface_i + icl_offset, 0);
				coarse_right.zero_patch   = iface_view(iface_i + zp_offset, 0);
				for (int q = 0; q < 4; q++) {
					coarse_right.global_i[q] = iface_view(iface_i + global_offset + q, 0);
				}
				ifaces.insert(coarse_right);
			}
			ifaces.insert(left);
			ifaces.insert(right);
		}
	}
	friend bool operator<(const RowIface &l, const RowIface &r)
	{
		return std::tie(l.main_i, l.global_i, l.type, l.s)
		       < std::tie(r.main_i, r.global_i, r.type, r.s);
	}
	friend bool operator==(const RowIface &l, const RowIface &r)
	{
		return std::tie(l.neumann, l.zero_patch) == std::tie(r.neumann, r.zero_patch);
	}
	/*friend bool operator!=(const ColIface &l, const ColIface &r)
	{
	    return std::tie(l.l_south, l.t_south, l.l_east, l.t_east, l.l_north, l.t_north, l.l_west,
	                    l.t_west)
	           == std::tie(r.l_south, r.t_south, r.l_west, r.t_west, r.l_north, r.t_north, r.l_east,
	                       r.t_east);
	}*/
	static void writeIfaces(Domain &d, int_vector_type &iface_info)
	{
		auto dist_view = iface_info.getLocalView<Kokkos::HostSpace>();
		// north
		Side iface_s = Side::north;
		do {
			if (d.hasNbr(iface_s)) {
				int iface_i = d.index(iface_s) * RowIface::size;
				if (!d.hasFineNbr(iface_s)) {
					dist_view(iface_i + coarse_right_offset, 0) = -1;
				}
				if ((iface_s == Side::south || iface_s == Side::west) && !d.hasFineNbr(iface_s)) {
					iface_i += right_offset;
				}
				dist_view(iface_i + main_offset, 0) = d.globalIndex(iface_s);
				dist_view(iface_i + side_offset, 0) = static_cast<int>(iface_s);
				dist_view(iface_i + type_offset, 0) = static_cast<int>(IfaceType::main);
				dist_view(iface_i + hfn_offset, 0)  = d.hasFineNbr(iface_s);
				dist_view(iface_i + hcn_offset, 0)  = d.hasCoarseNbr(iface_s);
				dist_view(iface_i + icl_offset, 0)  = d.isCoarseLeft(iface_s);
				dist_view(iface_i + zp_offset, 0)   = d.zero_patch;
				for (int q = 0; q < 4; q++) {
					dist_view(iface_i + global_offset + q, 0) = d.globalIndex(iface_s + q);
				}
			}
			if (d.hasCoarseNbr(iface_s)) {
				int       iface_i = d.indexCenter(iface_s) * RowIface::size;
				IfaceType t;
				if (d.isCoarseLeft(iface_s)) {
					t = IfaceType::left;
					iface_i += right_offset;
				} else {
					t = IfaceType::right;
					iface_i += coarse_right_offset;
				}
				dist_view(iface_i + main_offset, 0) = d.globalIndexCenter(iface_s);
				dist_view(iface_i + side_offset, 0) = static_cast<int>(iface_s);
				dist_view(iface_i + type_offset, 0) = static_cast<int>(t);
				dist_view(iface_i + hfn_offset, 0)  = d.hasFineNbr(iface_s);
				dist_view(iface_i + hcn_offset, 0)  = d.hasCoarseNbr(iface_s);
				dist_view(iface_i + icl_offset, 0)  = d.isCoarseLeft(iface_s);
				dist_view(iface_i + zp_offset, 0)   = d.zero_patch;
				for (int q = 0; q < 4; q++) {
					dist_view(iface_i + global_offset + q, 0) = d.globalIndex(iface_s + q);
				}
			}
			if (d.hasFineNbr(iface_s)) {
				int iface_i = d.indexRefinedLeft(iface_s) * RowIface::size;
				if (iface_s == Side::south || iface_s == Side::west) {
					iface_i += right_offset;
				}
				dist_view(iface_i + main_offset, 0) = d.globalIndexRefinedLeft(iface_s);
				dist_view(iface_i + side_offset, 0) = static_cast<int>(iface_s);
				dist_view(iface_i + type_offset, 0) = static_cast<int>(IfaceType::left);
				dist_view(iface_i + hfn_offset, 0)  = d.hasFineNbr(iface_s);
				dist_view(iface_i + hcn_offset, 0)  = d.hasCoarseNbr(iface_s);
				dist_view(iface_i + icl_offset, 0)  = d.isCoarseLeft(iface_s);
				dist_view(iface_i + zp_offset, 0)   = d.zero_patch;
				for (int q = 0; q < 4; q++) {
					dist_view(iface_i + global_offset + q, 0) = d.globalIndex(iface_s + q);
				}
				// right
				iface_i = d.indexRefinedRight(iface_s) * RowIface::size;
				if (iface_s == Side::south || iface_s == Side::west) {
					iface_i += right_offset;
				}
				dist_view(iface_i + main_offset, 0) = d.globalIndexRefinedRight(iface_s);
				dist_view(iface_i + side_offset, 0) = static_cast<int>(iface_s);
				dist_view(iface_i + type_offset, 0) = static_cast<int>(IfaceType::right);
				dist_view(iface_i + hfn_offset, 0)  = d.hasFineNbr(iface_s);
				dist_view(iface_i + hcn_offset, 0)  = d.hasCoarseNbr(iface_s);
				dist_view(iface_i + icl_offset, 0)  = d.isCoarseLeft(iface_s);
				dist_view(iface_i + zp_offset, 0)   = d.zero_patch;
				for (int q = 0; q < 4; q++) {
					dist_view(iface_i + global_offset + q, 0) = d.globalIndex(iface_s + q);
				}
			}
			iface_s++;
		} while (iface_s != Side::north);
	}
	BlockType getBlockType()
	{
		BlockType bt = BlockType::plain;
		if (hasCoarseNbr) {
			if (main_i == global_i[0]) {
				bt = BlockType::fine;
			} else {
				if (isCoarseLeft) {
					bt = BlockType::fine_out_left;
				} else {
					bt = BlockType::fine_out_right;
				}
			}
		}
		if (hasFineNbr) {
			if (main_i == global_i[0]) {
				bt = BlockType::coarse;
			}
			if (type == IfaceType::left) {
				bt = BlockType::coarse_out_left;
			}
			if (type == IfaceType::right) {
				bt = BlockType::coarse_out_right;
			}
		}
		return bt;
	}
};

class MatrixBlock
{
	public:
	MatrixBlock(int i, int j, bool flip_i, bool flip_j, bool right, std::bitset<4> neumann,
	            bool zero_patch, Side s, BlockType type)
	{
		this->i          = i;
		this->j          = j;
		this->flip_i     = flip_i;
		this->flip_j     = flip_j;
		this->right      = right;
		this->neumann    = neumann;
		this->zero_patch = zero_patch;
		this->type       = type;
		this->s          = s;
	}
	int            i, j;
	bool           flip_i, flip_j, right;
	bool           zero_patch;
	std::bitset<4> neumann;
	Side           s;
	BlockType      type;
	friend bool operator==(const MatrixBlock &l, const MatrixBlock &r)
	{
		return std::tie(l.neumann, l.zero_patch) == std::tie(r.neumann, r.zero_patch);
	}
	friend bool operator<(const MatrixBlock &l, const MatrixBlock &r)
	{
		return std::tie(l.j, l.i, l.right) < std::tie(r.j, r.i, r.right);
	}

	static void getBlocks(std::set<MatrixBlock> &blocks, const std::set<RowIface> &ifaces)
	{
		for (RowIface iface : ifaces) {
			Side absolute_s = iface.s;
			auto getrevxy   = [&](Side s, bool &reverse_x, bool &reverse_y) {
				switch (s) {
					case Side::north:
						reverse_x = false;
						reverse_y = false;
						break;
					case Side::east:
						reverse_x = true;
						reverse_y = false;
						break;
					case Side::south:
						reverse_x = true;
						reverse_y = true;
						break;
					case Side::west:
						reverse_x = false;
						reverse_y = true;
				}
			};
			bool reverse_x = false;
			bool reverse_y = false;

			int i = iface.main_i;
			// north
			{
				Side s = Side::north;

				std::bitset<4> neumann = iface.neumann;

				getrevxy(s + absolute_s, reverse_x, reverse_y);
				bool        right = (absolute_s + s == Side::south || absolute_s + s == Side::west);
				bool        flip_i = reverse_x;
				bool        flip_j = reverse_x;
				int         j      = iface.global_i[0];
				MatrixBlock b(i, j, flip_i, flip_j, right, neumann, iface.zero_patch, s,
				              iface.getBlockType());
				blocks.insert(b);
			}
			// east
			if (iface.global_i[1] != -1) {
				Side s = Side::west;

				std::bitset<4> neumann;
				for (int q = 0; q < 4; q++) {
					neumann[q] = iface.neumann[(q + 1) % 4];
				}

				getrevxy(s + absolute_s, reverse_x, reverse_y);
				bool        right = (absolute_s + s == Side::south || absolute_s + s == Side::west);
				bool        flip_i = !reverse_y;
				bool        flip_j = !reverse_x;
				int         j      = iface.global_i[1];
				MatrixBlock b(i, j, flip_i, flip_j, right, neumann, iface.zero_patch, s,
				              iface.getBlockType());
				blocks.insert(b);
			}
			// south
			if (iface.global_i[2] != -1) {
				Side s = Side::south;

				std::bitset<4> neumann;
				for (int q = 0; q < 4; q++) {
					neumann[q] = iface.neumann[(q + 2) % 4];
				}

				getrevxy(s + absolute_s, reverse_x, reverse_y);
				bool        right = (absolute_s + s == Side::south || absolute_s + s == Side::west);
				bool        flip_i = reverse_x;
				bool        flip_j = reverse_x;
				int         j      = iface.global_i[2];
				MatrixBlock b(i, j, flip_i, flip_j, right, neumann, iface.zero_patch, s,
				              iface.getBlockType());
				blocks.insert(b);
			}
			// west
			if (iface.global_i[3] != -1) {
				Side s = Side::east;

				std::bitset<4> neumann;
				for (int q = 0; q < 4; q++) {
					neumann[q] = iface.neumann[(q + 3) % 4];
				}

				getrevxy(s + absolute_s, reverse_x, reverse_y);
				bool        right = (absolute_s + s == Side::south || absolute_s + s == Side::west);
				bool        flip_i = !reverse_y;
				bool        flip_j = !reverse_x;
				int         j      = iface.global_i[3];
				MatrixBlock b(i, j, flip_i, flip_j, right, neumann, iface.zero_patch, s,
				              iface.getBlockType());
				blocks.insert(b);
			}
		}
	}
};
#endif
