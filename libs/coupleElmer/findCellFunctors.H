#ifndef findCellFunctors_H
#define findCellFunctors_H
namespace Foam
{
struct findCellFunctor
{
    const scalar xMin;
    const scalar xMax;
    const scalar yMin;
    const scalar yMax;
    const scalar zMin;
    const scalar zMax;
    const label cellNew;

	findCellFunctor(
			scalar _xMin,
			scalar _xMax,
			scalar _yMin,
			scalar _yMax,
			scalar _zMin,
			scalar _zMax,
			label _cellNew
		):
			xMin(_xMin),
			xMax(_xMax),
			yMin(_yMin),
			yMax(_yMax),
			zMin(_zMin),
			zMax(_zMax),
			cellNew(_cellNew)
		{}

	__HOST____DEVICE__
        //__device__
	scalar operator()(const vector& v1, const label& cellDefault)//const label& cellI, const vector& cC)
	{	//If point outside cell bounds then return list value, else return new index
		//bool outside = (v1[0] > xMax || v1[0] < xMin || v1[1] > yMax || v1[1] < yMin || v1[2] > zMax || v1[2] < zMin);
		//return cellDefault*outside+cellNew*(!outside);

		if(v1[0] > xMax || v1[0] < xMin || v1[1] > yMax || v1[1] < yMin || v1[2] > zMax || v1[2] < zMin)
		{
			//if(cellI1 > -1) return cellI1;
			//else return -1;
			return cellDefault;
		}
		else
		{
			return cellNew;
		}

	}
};
/*
struct copyLabelFunctor
{
    const scalar dummy;

	copyLabelFunctor() {}//scalar _dummy) : dummy(_dummy) {}

	__HOST____DEVICE__
        //__device__
	scalar operator()(const label& cellI, const label& cellI1)
	{
		if(cellI > -1)
		{
			return cellI;
		}
		else if(cellI1 > -1)
		{
			return cellI1;
		}
		else
		{
			return -1;
		}
	}
};
*/
}
#endif
