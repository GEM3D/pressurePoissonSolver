#ifndef IFACE_H
#define IFACE_H
enum class BlockType {
	north,
	east,
	south,
	west,
	// fine in
	north_fine_in,
	east_fine_in,
	south_fine_in,
	west_fine_in,
	// fine out
	// left
	north_fine_out_left,
	east_fine_out_left,
	south_fine_out_left,
	west_fine_out_left,
	// right
	north_fine_out_right,
	east_fine_out_right,
	south_fine_out_right,
	west_fine_out_right,
	// coarse in
	north_coarse_in,
	east_coarse_in,
	south_coarse_in,
	west_coarse_in,
	// coarse out
	// left
	north_coarse_out_left,
	east_coarse_out_left,
	south_coarse_out_left,
	west_coarse_out_left,
	// right
	north_coarse_out_right,
	east_coarse_out_right,
	south_coarse_out_right,
	west_coarse_out_right
};
class Blockk
{
	public:
	Blockk(int i, int j, bool flip_i, bool flip_j, bool right, std::bitset<4> neumann,
	       BlockType type)
	{
		this->i       = i;
		this->j       = j;
		this->flip_i  = flip_i;
		this->flip_j  = flip_j;
		this->right   = right;
		this->neumann = neumann;
		this->type    = type;
	}
	int            i, j;
	bool           flip_i, flip_j, right;
	std::bitset<4> neumann;
	BlockType      type;
	friend bool operator==(const Blockk &l, const Blockk &r)
	{
		return l.neumann == r.neumann;
	}
	friend bool operator<(const Blockk &l, const Blockk &r)
	{
		return std::tie(l.j, l.i, l.right) < std::tie(r.j, r.i, r.right);
	}
};
class Iface
{
	public:
	const static int left_offset  = 4;
	const static int side_size    = 6 + 2;
	const static int right_offset = 28 + 2 * 4;
	const static int size         = 46 + 2 * 7;
	bool             right;
	int              axis;
	int              refine_level = 1;
	std::array<int, 4> global_i      = {-1, -1, -1, -1};
	std::array<int, 4> center_i      = {-1, -1, -1, -1};
	std::array<int, 4> refined_left  = {-1, -1, -1, -1};
	std::array<int, 4> refined_right = {-1, -1, -1, -1};
	std::bitset<4> neumann;
	std::bitset<4> hasCoarseNbr;
	std::bitset<4> hasFineNbr;
	std::bitset<4> isCoarseLeft;
    void setNeumann(){
		for (int q = 0; q < 4; q++) {
			neumann[q] = global_i[q] == -1;
		}
    }
	static void readIfaces(std::set<Iface> &ifaces, int_vector_type &iface_info)
	{
		auto iface_view = iface_info.getLocalView<Kokkos::HostSpace>();
		for (size_t i = 0; i < iface_view.dimension(0); i += size) {
			Iface left;
			Iface right;
			int blr = iface_view(i, 0);

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
	friend bool operator<(const Iface &l, const Iface &r)
	{
		return std::tie(l.global_i[0], l.right) < std::tie(r.global_i[0], r.right);
	}
	friend bool operator==(const Iface &l, const Iface &r)
	{
		return l.neumann == r.neumann;
	}
	/*friend bool operator!=(const Iface &l, const Iface &r)
	{
		return std::tie(l.l_south, l.t_south, l.l_east, l.t_east, l.l_north, l.t_north, l.l_west,
		                l.t_west)
		       == std::tie(r.l_south, r.t_south, r.l_west, r.t_west, r.l_north, r.t_north, r.l_east,
		                   r.t_east);
	}*/
	friend std::ostream &operator<<(std::ostream &os, const Iface &iface)
	{
		os << "N:     " << iface.global_i[0] << std::endl;
		os << "E:     " << iface.global_i[1] << std::endl;
		os << "S:     " << iface.global_i[2] << std::endl;
		os << "W:     " << iface.global_i[3] << std::endl;
		os << "AXIS:  " << iface.axis << std::endl;
		os << "RIGHT: " << iface.right << std::endl;
		return os;
	}
    static void writeIfaces(Domain &d, int_vector_type &iface_info){
		auto dist_view = iface_info.getLocalView<Kokkos::HostSpace>();
		// north
		Side iface_s = Side::north;
		do {
			if (d.hasNbr(iface_s)) {
				int iface_i = d.globalIndex(iface_s) * Iface::size;
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
					i      = Iface::left_offset;
					do {
						dist_view(iface_i + i, 0)     = d.globalIndex(s);
						dist_view(iface_i + i + 1, 0) = d.globalIndexCenter(s);
						dist_view(iface_i + i + 2, 0) = d.hasFineNbr(s);
						dist_view(iface_i + i + 3, 0) = d.hasCoarseNbr(s);
						dist_view(iface_i + i + 4, 0) = d.isCoarseLeft(s);
						dist_view(iface_i + i + 5, 0) = d.globalIndexRefinedLeft(s);
						dist_view(iface_i + i + 6, 0) = d.globalIndexRefinedRight(s);
						i += Iface::side_size;
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
						dist_view(iface_i + 2, 0) = d.ds.refine_level;
						dist_view(iface_i + left_offset, 0)     = d.globalIndex(iface_s);
						dist_view(iface_i + left_offset + 1, 0) = d.globalIndexCenter(iface_s);
						dist_view(iface_i + left_offset + 2, 0) = d.hasFineNbr(iface_s);
						dist_view(iface_i + left_offset + 3, 0) = d.hasCoarseNbr(iface_s);
						dist_view(iface_i + left_offset + 4, 0) = d.isCoarseLeft(iface_s);
						dist_view(iface_i + left_offset + 5, 0) = d.globalIndexRefinedLeft(iface_s);
						dist_view(iface_i + left_offset + 6, 0) = d.globalIndexRefinedRight(iface_s);
					}
					Side s = iface_s;
					s++;
					i = Iface::right_offset;
					do {
						dist_view(iface_i + i, 0)     = d.globalIndex(s);
						dist_view(iface_i + i + 1, 0) = d.globalIndexCenter(s);
						dist_view(iface_i + i + 2, 0) = d.hasFineNbr(s);
						dist_view(iface_i + i + 3, 0) = d.hasCoarseNbr(s);
						dist_view(iface_i + i + 4, 0) = d.isCoarseLeft(s);
						dist_view(iface_i + i + 5, 0) = d.globalIndexRefinedLeft(s);
						dist_view(iface_i + i + 6, 0) = d.globalIndexRefinedRight(s);
						i += Iface::side_size;
						s++;
					} while (s != iface_s);
				}
			}
			iface_s++;
		} while (iface_s != Side::north);
	}
};
#endif
