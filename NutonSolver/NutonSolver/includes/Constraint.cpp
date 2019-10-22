#include "Constraint.h"

static	double	ft_Stopper(Point A1, Point A2, Point B1, Point B2, double _value)
{
	return (0.0f);
}

void			Constraint::Fill_Free()
{
	pointA1 = &tmp_pointA1;
	pointA2 = &tmp_pointA2;
	pointB1 = &tmp_pointB1;
	pointB2 = &tmp_pointB2;
	value = 0;

	func = &ft_Stopper;

	func_dA1x = &ft_Stopper;
	func_dA1y = &ft_Stopper;
	func_dA2x = &ft_Stopper;
	func_dA2y = &ft_Stopper;
	func_dB1x = &ft_Stopper;
	func_dB1y = &ft_Stopper;
	func_dB2x = &ft_Stopper;
	func_dB2y = &ft_Stopper;

	//	A1x_
	func_dA1x_dA1x = &ft_Stopper;
	func_dA1x_dA1y = &ft_Stopper;
	func_dA1x_dA2x = &ft_Stopper;
	func_dA1x_dA2y = &ft_Stopper;
	func_dA1x_dB1x = &ft_Stopper;
	func_dA1x_dB1y = &ft_Stopper;
	func_dA1x_dB2x = &ft_Stopper;
	func_dA1x_dB2y = &ft_Stopper;
	// A1y_
	func_dA1y_dA1x = &ft_Stopper;
	func_dA1y_dA1y = &ft_Stopper;
	func_dA1y_dA2x = &ft_Stopper;
	func_dA1y_dA2y = &ft_Stopper;
	func_dA1y_dB1x = &ft_Stopper;
	func_dA1y_dB1y = &ft_Stopper;
	func_dA1y_dB2x = &ft_Stopper;
	func_dA1y_dB2y = &ft_Stopper;
	//
	//	A2x_
	func_dA2x_dA1x = &ft_Stopper;
	func_dA2x_dA1y = &ft_Stopper;
	func_dA2x_dA2x = &ft_Stopper;
	func_dA2x_dA2y = &ft_Stopper;
	func_dA2x_dB1x = &ft_Stopper;
	func_dA2x_dB1y = &ft_Stopper;
	func_dA2x_dB2x = &ft_Stopper;
	func_dA2x_dB2y = &ft_Stopper;
	// A1y_
	func_dA2y_dA1x = &ft_Stopper;
	func_dA2y_dA1y = &ft_Stopper;
	func_dA2y_dA2x = &ft_Stopper;
	func_dA2y_dA2y = &ft_Stopper;
	func_dA2y_dB1x = &ft_Stopper;
	func_dA2y_dB1y = &ft_Stopper;
	func_dA2y_dB2x = &ft_Stopper;
	func_dA2y_dB2y = &ft_Stopper;
	//
	//	B1x_
	func_dB1x_dA1x = &ft_Stopper;
	func_dB1x_dA1y = &ft_Stopper;
	func_dB1x_dA2x = &ft_Stopper;
	func_dB1x_dA2y = &ft_Stopper;
	func_dB1x_dB1x = &ft_Stopper;
	func_dB1x_dB1y = &ft_Stopper;
	func_dB1x_dB2x = &ft_Stopper;
	func_dB1x_dB2y = &ft_Stopper;
	// B1y_
	func_dB1y_dA1x = &ft_Stopper;
	func_dB1y_dA1y = &ft_Stopper;
	func_dB1y_dA2x = &ft_Stopper;
	func_dB1y_dA2y = &ft_Stopper;
	func_dB1y_dB1x = &ft_Stopper;
	func_dB1y_dB1y = &ft_Stopper;
	func_dB1y_dB2x = &ft_Stopper;
	func_dB1y_dB2y = &ft_Stopper;
	//
	//	B2x_
	func_dB2x_dA1x = &ft_Stopper;
	func_dB2x_dA1y = &ft_Stopper;
	func_dB2x_dA2x = &ft_Stopper;
	func_dB2x_dA2y = &ft_Stopper;
	func_dB2x_dB1x = &ft_Stopper;
	func_dB2x_dB1y = &ft_Stopper;
	func_dB2x_dB2x = &ft_Stopper;
	func_dB2x_dB2y = &ft_Stopper;
	// B2y_
	func_dB2y_dA1x = &ft_Stopper;
	func_dB2y_dA1y = &ft_Stopper;
	func_dB2y_dA2x = &ft_Stopper;
	func_dB2y_dA2y = &ft_Stopper;
	func_dB2y_dB1x = &ft_Stopper;
	func_dB2y_dB1y = &ft_Stopper;
	func_dB2y_dB2x = &ft_Stopper;
	func_dB2y_dB2y = &ft_Stopper;
}

Constraint::Constraint()
{
	Fill_Free();
}

Constraint::~Constraint()
{

}

double	Constraint::Function()
{
	return func(*pointA1, *pointA2, *pointB1, *pointB2, value);
}

double	Constraint::Derivative(Point* point, bool xy)
{
	if (point == pointA1)
	{
		if (xy)
			return (func_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
		return (func_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
	}
	if (point == pointA2)
	{
		if (xy)
			return (func_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
		return (func_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
	}
	if (point == pointB1)
	{
		if (xy)
			return (func_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
		return (func_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
	}
	if (point == pointB2)
	{
		if (xy)
			return (func_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
		return (func_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
	}

	return (0.0f);
}

double	Constraint::SecondDerivative(Point* point_d1, bool xy_d1, Point* point_d2, bool xy_d2)
{
	if (point_d1 == pointA1)
	{
		if (point_d2 == pointA1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA1x_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA1x_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA1y_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA1y_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointA2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA1x_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA1x_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA1y_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA1y_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA1x_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA1x_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA1y_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA1y_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA1x_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA1x_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA1y_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA1y_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
	}
	if (point_d1 == pointA2)
	{
		if (point_d2 == pointA1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA2x_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA2x_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA2y_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA2y_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointA2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA2x_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA2x_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA2y_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA2y_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA2x_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA2x_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA2y_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA2y_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dA2x_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dA2x_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dA2y_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dA2y_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
	}
	if (point_d1 == pointB1)
	{
		if (point_d2 == pointA1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB1x_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB1x_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB1y_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB1y_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointA2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB1x_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB1x_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB1y_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB1y_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB1x_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB1x_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB1y_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB1y_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB1x_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB1x_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB1y_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB1y_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
	}
	if (point_d1 == pointB2)
	{
		if (point_d2 == pointA1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB2x_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB2x_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB2y_dA1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB2y_dA1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointA2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB2x_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB2x_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB2y_dA2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB2y_dA2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB1)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB2x_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB2x_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB2y_dB1x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB2y_dB1y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
		if (point_d2 == pointB2)
		{
			if (xy_d1)
			{
				if (xy_d2)
					return (func_dB2x_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
				return (func_dB2x_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
			}
			if (xy_d2)
				return (func_dB2y_dB2x(*pointA1, *pointA2, *pointB1, *pointB2, value));
			return (func_dB2y_dB2y(*pointA1, *pointA2, *pointB1, *pointB2, value));
		}
	}

	return (0.0f);
}