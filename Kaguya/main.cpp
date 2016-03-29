#include "BlackBox.h"
//=============================================================================
int main(int argc, char **argv)
{
	BlackBox::initialize(&argc, argv);
	BlackBox::run();
	BlackBox::save();
	BlackBox::destroy();
	
	return 0;
}
//=============================================================================