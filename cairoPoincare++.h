/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file cairoPoincare.h
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2010/06/01
 **/

#include <complex>
#include "cairoPoincare.h"


/* geodesic midpoint in C++
    double cx,cy,c_r;
    computeCircleParameters(&a, &b, &cx, &cy, &c_r);
    std::complex<double> c = std::complex<double>(cx,cy);
    double D_theta = std::arg(c);
    double D_r = std::abs(c);
    //std::cout<<"theta = "<<D_theta<<", r = "<<D_r<<std::endl;
    drawPointCairo(cr, new Point(std::polar(D_r-c_r, D_theta)));
*/

struct Point : point {

  Point() {}

  Point(double xx, double yy) { X = xx; Y = yy; }

  Point(std::complex<double> p) { X = p.real(); Y = p.imag(); }

  double x() const
  {return X;}

  double y() const
  {return Y;}

};

void drawLine(const Point &a, const Point &b)
{
  drawLineCairo(cairo, &a, &b, true);
}
void drawLine(cairo_t *cr, const Point &a, const Point &b)
{
  drawLineCairo(cr, &a, &b, true);
}

void drawUnitCircle(cairo_t *cr)
{
  drawUnitCircleCairo(cr);
}
