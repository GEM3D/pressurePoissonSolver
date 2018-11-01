#include "PBMatrix.h"
#include "SchurHelper.h"
struct Block;
class SchurMatrixHelper
{
	private:
	std::shared_ptr<SchurHelper<3>> sh;
	int                             n;

	typedef std::function<void(Block *, std::shared_ptr<std::valarray<double>>)> inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	SchurMatrixHelper(std::shared_ptr<SchurHelper<3>> sh)
	{
		this->sh = sh;
	}
	PW_explicit<Mat> formCRSMatrix();
	PBMatrix *       formPBMatrix();
	void             getPBDiagInv(PC p);
	PW_explicit<Mat> getPBMatrix();
	PW_explicit<Mat> getPBDiagInv();
};
