#ifndef COEFFHELPER_H
#define COEFFHELPER_H
class CoeffHelper
{
	public:
	int    curr_index[2];
	int    size;
	double values[5];
	int    indices[5];
	int    off_size;
	double off_values[5];
	int    off_indices[5];
	CoeffHelper(Domain *d,Side s){
        this->s=s;
        this->d=d;
		n      = d->n;
		switch (s) {
			case Side::north:
				curr_index[0] = 0;
				curr_index[1] = n - 1;
				break;
			case Side::east:
				curr_index[0] = n - 1;
				curr_index[1] = n - 1;
				break;
			case Side::south:
				curr_index[0] = n - 1;
				curr_index[1] = 0;
				break;
			case Side::west:
				curr_index[0] = 0;
				curr_index[1] = 0;
		}
		updateValues();
		updateIndices();
	}
	CoeffHelper &operator++()
	{
		switch (s) {
			case Side::north:
				if (curr_index[0] < n - 1) {
					curr_index[0]++;
					i++;
				} else {
					is_done = true;
				}
				break;
			case Side::east:
				if (curr_index[1] > 0) {
					curr_index[1]--;
					i++;
				} else {
					is_done = true;
				}
				break;
			case Side::south:
				if (curr_index[0] > 0) {
					curr_index[0]--;
					i++;
				} else {
					is_done = true;
				}
				break;
			case Side::west:
				if (curr_index[1] < n - 1) {
					curr_index[1]++;
					i++;
				} else {
					is_done = true;
				}
		}
		updateValues();
		updateIndices();
		return *this;
	}

	bool done() { return is_done; }
	private:
	Side    s;
	Domain *d;
	int     n;
	int     i       = 0;
	bool    is_done = false;
	void    updateValues()
	{
		if (d->hasFineNbr(s)) {
			if (i == 0) {
				values[0] = 0;
				values[1] = -21.0 / 10.0;
				values[2] = 1.0;
				values[3] = 0;
				values[4] = 1.0 / 15.0;
				off_size  = 5;
			} else if (i == n - 1) {
				values[0] = 0;
				values[1] = -21.0 / 10.0;
				values[2] = 1.0;
				values[3] = 1.0 / 15.0;
				values[4] = 0;
				off_size  = 5;
			} else {
				values[0] = 0;
				values[1] = -2.0;
				values[2] = 1.0;
				values[3] = -1.0 / 30.0;
				values[4] = -1.0 / 30.0;
				off_size  = 4;
			}
			off_values[0] = 1.0 / 3.0;
			off_values[1] = 1.0 / 5.0;
			off_values[2] = 1.0 / 3.0;
			off_values[3] = 1.0 / 5.0;
			off_values[4] = -1.0 / 30.0;

			size     = 5;

		} else if (d->hasCoarseNbr(s)) {
			values[0] = 0;
			values[1] = -4.0 / 3.0;
			values[2] = 4.0 / 5.0;
			if (!d->isCoarseLeft(s) && i < 2) {
				if (i == 0) {
					off_values[0] = 3.0 / 4.0;
					off_values[1] = -3.0 / 10.0;
					off_values[2] = 1.0 / 12.0;
				} else {
					off_values[0] = 7.0 / 20.0;
					off_values[1] = 7.0 / 30.0;
					off_values[2] = -1.0 / 20.0;
				}
			} else if (d->isCoarseLeft(s) && i > n - 3) {
				if (i == n - 1) {
					off_values[0] = 1.0 / 12.0;
					off_values[1] = -3.0 / 10.0;
					off_values[2] = 3.0 / 4.0;
				} else {
					off_values[0] = -1.0 / 20.0;
					off_values[1] = 7.0 / 30.0;
					off_values[2] = 7.0 / 20.0;
				}
			} else {
				if (i % 2 == 0) {
					off_values[0] = 1.0 / 12.0;
					off_values[1] = 1.0 / 2.0;
					off_values[2] = -1.0 / 20.0;
				} else {
					off_values[0] = -1.0 / 20.0;
					off_values[1] = 1.0 / 2.0;
					off_values[2] = 1.0 / 12.0;
				}
			}
			size     = 3;
			off_size = 3;
		} else if (d->hasNbr(s)) {
			values[0] = 1;
			values[1] = -2;
			values[2] = 1;
			size      = 3;
			off_size  = 0;
		} else if (d->neumann) {
			values[0] = 0;
			values[1] = -1;
			values[2] = 1;
			size      = 3;
			off_size  = 0;
		} else {
			values[0] = 0;
			values[1] = -3;
			values[2] = 1;
			size      = 3;
			off_size  = 0;
		}
		for (int i = 0; i < 5; i++) {
			if (s == Side::north || s == Side::south) {
				values[i] /= (d->h_y * d->h_y);
				off_values[i] /= (d->h_y * d->h_y);
			} else {
				values[i] /= (d->h_x * d->h_x);
				off_values[i] /= (d->h_x * d->h_x);
			}
		}
        
	}
	void updateIndices()
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		switch (s) {
			case Side::north:
				indices[0] = 4;
				indices[1] = 0;
				indices[2] = 3;
				indices[3] = 1;
				indices[4] = 2;
				break;
			case Side::east:
				indices[0] = 2;
				indices[1] = 0;
				indices[2] = 1;
				indices[3] = 4;
				indices[4] = 3;
				break;
			case Side::south:
				indices[0] = 3;
				indices[1] = 0;
				indices[2] = 4;
				indices[3] = 2;
				indices[4] = 1;
				break;
			case Side::west:
				indices[0] = 1;
				indices[1] = 0;
				indices[2] = 2;
				indices[3] = 3;
				indices[4] = 4;
		}
		if (d->hasFineNbr(s)) {
			if (i == 0) {
				off_indices[0] = 6;
				off_indices[1] = 7;
				off_indices[2] = 8;
				off_indices[3] = 9;
				off_indices[4] = 5;
                //corner
				if (s != Side::north) {
					Side s2 = s;
					s2--;
					int shift = 0;
					if (d->hasCoarseNbr(s2)) {
						shift = 3;
					} else if (d->hasFineNbr(s2)) {
						shift = 5;
					}
					for (int i = 0; i < 5; i++) {
						off_indices[i] += shift;
					}
				}
			} else if (i == n - 1) {
				off_indices[0] = 5;
				off_indices[1] = 6;
				off_indices[2] = 7;
				off_indices[3] = 8;
				off_indices[4] = 9;
				// last corner
				if (s == Side::west) {
					int shift = 0;
					if (d->hasCoarseNbr(Side::north)) {
						shift = 3;
					} else if (d->hasFineNbr(Side::north)) {
						shift = 5;
					}
					for (int i = 0; i < 5; i++) {
						off_indices[i] += shift;
					}
				}
			} else {
				off_indices[0] = 5;
				off_indices[1] = 6;
				off_indices[2] = 7;
				off_indices[3] = 8;
			}

		} else if (d->hasCoarseNbr(s)) {
			off_indices[0] = 5;
			off_indices[1] = 6;
			off_indices[2] = 7;
			if (i == 0) {
				// corner
				if (s != Side::north) {
					Side s2 = s;
					s2--;
					int shift = 0;
					if (d->hasCoarseNbr(s2)) {
						shift = 3;
					} else if (d->hasFineNbr(s2)) {
						shift = 5;
					}
					for (int i = 0; i < 5; i++) {
						off_indices[i] += shift;
					}
				}
			}
			if (i == n - 1) {
				// last corner
				if (s == Side::west) {
					int shift = 0;
					if (d->hasCoarseNbr(Side::north)) {
						shift = 3;
					} else if (d->hasFineNbr(Side::north)) {
						shift = 5;
					}
					for (int i = 0; i < 5; i++) {
						off_indices[i] += shift;
					}
				}
			}
		}
	}
};
#endif
