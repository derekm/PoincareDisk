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

#include<cairo.h>
#include<cairo-pdf.h>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <float.h>
#include <stdbool.h>


//Size of the cairo board
#define SIZEX 1000
#define SIZEY  1000
//Radius of the "unit" disk
#define RAD 400.0

//#define EDGEWIDTH 0.5
#define EDGEWIDTH 2
#define POINTRADIUS 3


//Global drawing board
cairo_surface_t *CSglobal;
cairo_t *cairo;



struct point {
  double X, Y;
};


//Functions to map a point from the unit disk to the cairo board
double cX(double x)
{
  return (x*RAD+SIZEX/2.0);
}
double cY(double y)
{
  return (SIZEY-y*RAD) - SIZEY/2.0;
}


/**
 * Draw a point in the Poincaré Disk
 *
 * @param cr the Cairo context on which to draw
 * @param a the point to plot
 * @param r red value
 * @param g green value
 * @param b blue value
 * @param alpha transparency value
 */
void drawPointCairoColor(cairo_t *cr, struct point const *const a,
                         const double r, const double g, const double b, const double alpha)
{
  cairo_set_source_rgba (cr, r, g, b, alpha);
  cairo_set_line_width (cr, 0.5);
  cairo_arc(cr, cX(a->X),cY(a->Y), POINTRADIUS, 0, 2*M_PI);
  cairo_fill (cr);
  cairo_stroke(cr);
}

const double red=0.0, green=0.0, blue=0.0, alpha=0.5;

void drawPointCairo(cairo_t *cr, struct point const *const a)
{
  drawPointCairoColor(cr, a, red, green, blue, alpha);
}

void drawPoint(struct point const *const a)
{
  drawPointCairoColor(cairo, a, red, green, blue, alpha);
}



/**
 * Compute the circle containing points a and b, and orthogonal to the
 * unit circle (a.k.a shortest path between a and b in the hyperbolic sense).
 *
 * @param a input point a
 * @param b input point b
 * @param cx abscissa of the circle center
 * @param cy ordiniate of the circle center
 * @param radius radius of the circle
 *
 * @return true if the output is a circle. False if the shortest path
 * is a straight line (diameter of the unit circle).
 */
bool computeCircleParameters(struct point const *const a, struct point const *const b,
                             double *const cx, double *const cy, double *const radius)
{
  const double ax=a->X, ay=a->Y, bx=b->X, by=b->Y;

  if (fabs(ax*by - ay*bx) < DBL_EPSILON)
    return false;

  *cx = (-1/2.0*(ay*bx*bx + ay*by*by
		- (ax*ax + ay*ay + 1)*by + ay)
	/(ax*by -ay*bx));

  *cy = (1/2.0*(ax*bx*bx + ax*by*by
	       - (ax*ax + ay*ay + 1)*bx + ax)
	/(ax*by - ay*bx));


  *radius = (-1/2.0*sqrt(ax*ax*ax*ax*bx*bx + ax*ax*ax*ax*by*by - 2*ax*ax*ax*bx*bx*bx-
			2*ax*ax*ax*bx*by*by + 2*ax*ax*ay*ay*bx*bx + 2*ax*ax*ay*ay*by*by - 2*ax*ax*ay*bx*bx*by
			- 2*ax*ax*ay*by*by*by + ax*ax*bx*bx*bx*bx + 2*ax*ax*bx*bx*by*by + ax*ax*by*by*by*by -
			2*ax*ay*ay*bx*bx*bx - 2*ax*ay*ay*bx*by*by + ay*ay*ay*ay*bx*bx + ay*ay*ay*ay*by*by -
			2*ay*ay*ay*bx*bx*by - 2*ay*ay*ay*by*by*by + ay*ay*bx*bx*bx*bx + 2*ay*ay*bx*bx*by*by + ay*ay*by*by*by*by
			- 2*ax*ax*ax*bx - 2*ax*ax*ay*by + 4*ax*ax*bx*bx - 2*ax*ay*ay*bx + 8*ax*ay*bx*by
			- 2*ax*bx*bx*bx - 2*ax*bx*by*by - 2*ay*ay*ay*by + 4*ay*ay*by*by - 2*ay*bx*bx*by -
			2*ay*by*by*by + ax*ax - 2*ax*bx + ay*ay - 2*ay*by + bx*bx + by*by)
	    /(ax*by -	ay*bx) );

  if ((*radius) < 0.0)
    *radius= 1/2.0*sqrt(ax*ax*ax*ax*bx*bx + ax*ax*ax*ax*by*by - 2*ax*ax*ax*bx*bx*bx
		       - 2*ax*ax*ax*bx*by*by + 2*ax*ax*ay*ay*bx*bx + 2*ax*ax*ay*ay*by*by -
		       2*ax*ax*ay*bx*bx*by - 2*ax*ax*ay*by*by*by + ax*ax*bx*bx*bx*bx + 2*ax*ax*bx*bx*by*by +
		       ax*ax*by*by*by*by - 2*ax*ay*ay*bx*bx*bx - 2*ax*ay*ay*bx*by*by + ay*ay*ay*ay*bx*bx + ay*ay*ay*ay*by*by -
		       2*ay*ay*ay*bx*bx*by - 2*ay*ay*ay*by*by*by + ay*ay*bx*bx*bx*bx + 2*ay*ay*bx*bx*by*by + ay*ay*by*by*by*by
		       - 2*ax*ax*ax*bx - 2*ax*ax*ay*by + 4*ax*ax*bx*bx - 2*ax*ay*ay*bx + 8*ax*ay*bx*by
		       - 2*ax*bx*bx*bx - 2*ax*bx*by*by - 2*ay*ay*ay*by + 4*ay*ay*by*by - 2*ay*bx*bx*by -
		       2*ay*by*by*by + ax*ax - 2*ax*bx + ay*ay - 2*ay*by + bx*bx + by*by)/(ax*by - ay*bx);

  return true;

}


/**
 * Compute the angles of the points (a.k.a., endpoints, ideal points, omega
 * points, etc.) at which the circle meets the boundary of the unit circle,
 * using the circle parameters.
 *
 */
void computeOmegaPoints(double const *const cx, double const *const cy, double const *const r,
                        double *const theta1, double *const theta2)
{
  const double dtheta = atan2(*r, 1);
  double theta = atan2(*cy, *cx);
  if (theta<0) theta+=2*M_PI;
  *theta1 = theta - dtheta;
  *theta2 = theta + dtheta;
}


/**
 * Compute circle-line intersect for a given line and the unit circle.
 */
void computeCircleLineIntersect(const double r, struct point const *const p1, struct point const *const p2,
                                struct point *const meet1, struct point *const meet2)
{
	// Lua code
    //dx = l.p2.x - l.p1.x
    //dy = l.p2.y - l.p1.y
    //dr = math.sqrt(dx^2 + dy^2)
    //D = l.p1.x*l.p2.y - l.p2.x*l.p1.y
    //numer1 = D * dy
    //numer3 = math.sqrt(self.r^2 * dr^2 - D^2)
    //numer2 = sign(dy)*dx * numer3
    //denom = dr^2
    //x1 = (numer1 + numer2) / denom
    //x2 = (numer1 - numer2) / denom
    //numer1 = -(D * dx)
    //numer2 = math.abs(dy) * numer3
    //y1 = (numer1 + numer2) / denom
    //y2 = (numer1 - numer2) / denom
    //return { vec2( x1, y1 ), vec2( x2, y2 ) }

	const double dx = p2->X - p1->X;
	const double dy = p2->Y - p1->Y;
	const double dr = sqrt(pow(dx, 2.) + pow(dy, 2.));
	const double D = p1->X * p2->Y - p2->X * p1->Y;
	      double numer1 = D * dy;
	const double numer3 = sqrt(pow(r, 2.) * pow(dr, 2.) - pow(D, 2.));
	      double numer2 = (dy < 0 ? -1 : 1) * dx * numer3;
	const double denom = pow(dr, 2.);

	meet1->X = (numer1 + numer2) / denom;
	meet2->X = (numer1 - numer2) / denom;

	numer1 = -(D * dx);
	numer2 = abs(dy) * numer3;

	meet1->Y = (numer1 + numer2) / denom;
	meet2->Y = (numer1 - numer2) / denom;
}


/**
 * Compute the omega points from an arbitrary/assumed midpoint.
 */
void computeOmegaPointsFromMidpoint(struct point const *const m,
                                    struct point *const pole1, struct point *const pole2,
                                    struct point *const omega1, struct point *const omega2)
{
	const double phi = fmod(atan2(m->Y, m->X) + 2*M_PI, 2*M_PI);

	double polar1, polar2;
	if (phi < M_PI/2) {
		polar1 = phi + M_PI/2;
		polar2 = phi + M_PI + M_PI/2;
	} else if (phi > M_PI + M_PI/2) {
		polar1 = phi - (M_PI + M_PI/2);
		polar2 = phi - M_PI/2;
	} else {
		polar1 = phi + M_PI/2;
		polar2 = phi - M_PI/2;
	}

	pole1->X = cos(polar1);
	pole1->Y = sin(polar1);
	pole2->X = cos(polar2);
	pole2->Y = sin(polar2);

	struct point meet1a, meet1b;
	computeCircleLineIntersect(1., pole1, m, &meet1a, &meet1b);

	struct point meet2a, meet2b;
	computeCircleLineIntersect(1., pole2, m, &meet2a, &meet2b);

	const double meet1a_phi = fmod(atan2(meet1a.Y, meet1a.X) + 2*M_PI, 2*M_PI);
	if (meet1a_phi == polar1) {
		omega1->X = meet1b.X;
		omega1->Y = meet1b.Y;
	} else {
		omega1->X = meet1a.X;
		omega1->Y = meet1a.Y;
	}

	const double meet2a_phi = fmod(atan2(meet2a.Y, meet2a.X) + 2*M_PI, 2*M_PI);
	if (meet2a_phi == polar2) {
		omega2->X = meet2b.X;
		omega2->Y = meet2b.Y;
	} else {
		omega2->X = meet2a.X;
		omega2->Y = meet2a.Y;
	}
}


/**
 * Compute the angles of the intersect points from the perspective of the
 * intersecting circle, since that is the perspective from which cairo_arc()
 * draws.
 *
 */
void computeIntersectAngles(double const *const cx, double const *const cy, double const *const r,
                            double *const phi1, double *const phi2)
{
  const double dphi = atan2(1, *r);
  double itheta = atan2(*cy, *cx);
  if (itheta < 0) itheta += 2*M_PI;
  if (itheta >= M_PI) itheta -= M_PI;
  else if (itheta < M_PI) itheta += M_PI;
  itheta = 2*M_PI - itheta; // convert to Cairo's polar orientation
  *phi1 = itheta - dphi;
  *phi2 = itheta + dphi;
}

void computeMidpoint(double const *const cx, double const *const cy, double const *const r,
                     double *const R, double *const theta)
{
  *R = sqrt(pow(*cx, 2.) + pow(*cy, 2.)) - *r;
  *theta = atan2(*cy, *cx);
}

/**
 * Draw the hyperbolic line through a and b.
 *
 * @param cr the Cairo context on which to draw
 * @param a first point.
 * @param b second point.
 * @param withPoint true means that  points are also displayed.
 */
void drawLineCairo(cairo_t *cr, struct point const *const a, struct point const *const b, bool withPoint)
{
  if (withPoint)
    {
      drawPointCairo(cr, a);
      drawPointCairo(cr, b);
    }

  const double ax=a->X, ay=a->Y, bx=b->X, by=b->Y;
  double cx,cy,r;

  bool result = computeCircleParameters(a,b,&cx,&cy,&r);

  if (not(result))
    {
      cairo_set_source_rgba (cr, 0, 0.6, 0, 0.5);
      cairo_set_line_width (cr, 0.5);

      //We project a and b points onto the unit circle.
      const double theta = atan2(ay,ax);
      const double theta2 = atan2(by,bx);

      cairo_move_to (cr, cX(cos(theta)),cY(sin(theta)));
      cairo_line_to (cr, cX(cos(theta2)),cY(sin(theta2)));
      cairo_stroke(cr);
    }
  else
    {
      double theta1, theta2;
      computeIntersectAngles(&cx, &cy, &r, &theta1, &theta2);
      cairo_set_source_rgba (cr, 0, 0.6, 0, 0.5);
      cairo_set_line_width (cr, EDGEWIDTH);
      cairo_arc(cr, cX(cx),cY(cy), r*RAD, theta1, theta2);
      cairo_stroke(cr);
    }
}

void drawLine(struct point const *const a, struct point const *const b, bool withPoint)
{
  drawLineCairo(cairo, a, b, withPoint);
}


/**
 * Internal method for the edge drawing.
 *
 * @param a first point.
 * @param b second point.
 * @param withPoint if true, points are displayed.
 */
void internaldrawEdgeCairo(cairo_t *cr, struct point const *const a, struct point const *const b, bool withPoint)
{
  if (withPoint)
    {
      drawPointCairo(cr, a);
      drawPointCairo(cr, b);
    }

  const double ax=a->X, ay=a->Y, bx=b->X, by=b->Y;
  double cx,cy,r;
  bool result  = computeCircleParameters(a,b,&cx,&cy,&r);

  //Near-aligned points -> straight segment.
  if (not(result))
    {
      cairo_move_to (cr, cX(ax),cY(ay));
      cairo_line_to (cr, cX(bx),cY(by));
    }
  else
    {
      double theta  = atan2(-ay+cy,ax-cx);
      double theta2 = atan2(-by+cy,bx-cx);

      //We recenter the angles to [0, 2Pi]
      while (theta<0)
	theta += 2*M_PI;
      while (theta2<0)
	theta2 += 2*M_PI;

      if (theta > theta2)
	if ( theta - theta2 > M_PI)
	  cairo_arc(cr, cX(cx),cY(cy), r*RAD ,theta, theta2);
	else
	  cairo_arc_negative(cr, cX(cx),cY(cy), r*RAD ,theta, theta2);
      else
	if ( theta2 - theta < M_PI)
	  cairo_arc(cr, cX(cx),cY(cy), r*RAD ,theta, theta2);
	else
	  cairo_arc_negative(cr, cX(cx),cY(cy), r*RAD ,theta, theta2);
    }
}
void internaldrawEdge(struct point const *const a, struct point const *const b, bool withPoint)
{
  internaldrawEdgeCairo(cairo, a, b, withPoint);
}

/**
 * Main method to draw an hyperbolic edge between a and b
 *
 * @param a first point
 * @param b second point
 * @param withPoints if true, points are displayed.
 */
void drawEdgeCairo(cairo_t *cr, struct point const *const a, struct point const *const b, const bool withPoints)
{
  cairo_set_source_rgba (cr, 0, 0.0, 1, 0.8);
  cairo_set_line_width (cr, EDGEWIDTH);

  internaldrawEdgeCairo(cr,a,b,withPoints);
  cairo_stroke(cr);
}
void drawEdge(struct point const *const a, struct point const *const b, const bool withPoints)
{
  drawEdgeCairo(cairo, a, b, withPoints);
}
/**
 * Draw an hyperbolic triangle with vertices (a,b,c)..
 *
 * @param a first point
 * @param b second point
 * @param c third point
 * @param r red color for the triangle
 * @param g green color
 * @param bl blue color
 * @param withPoint if true, points are displayed.
 * @param withSupportLines if true, supporting line are displayed.
 */
void drawTriangle(struct point const *const a, struct point const *const b, struct point const *const c,
		  const double r=0, const double g=0, const double bl=1.0,
		  const bool withPoint=true, const bool withSupportLines=true)
{
  if (withPoint)
    {
      drawPoint(a);
      drawPoint(b);
      drawPoint(c);
    }

  if (withSupportLines)
    {
      drawLine(a,b,false);
      drawLine(b,c,false);
      drawLine(c,a,false);
    }

  cairo_set_source_rgba (cairo, 0, 0.0, 1, 0.8);
  cairo_set_line_width (cairo, EDGEWIDTH);

  cairo_new_path (cairo);
  internaldrawEdge(a,b,false);
  internaldrawEdge(b,c,false);
  internaldrawEdge(c,a,false);


  cairo_set_source_rgba (cairo, r, g, bl, 0.2);
  cairo_fill_preserve(cairo);

  cairo_set_source_rgba (cairo, 0.0, 0.0, 1, 0.8);
  cairo_set_line_width (cairo, EDGEWIDTH);
  cairo_stroke(cairo);
}

/**
 * Draw the unit circle and clip every cairo drawning to it.
 *
 */
void drawUnitCircleCairo(cairo_t *cr)
{
  cairo_arc(cr, cX(0),cY(0), RAD, 0, 2*M_PI);
  cairo_clip (cr);
  cairo_new_path (cr); /* path not consumed by clip()*/

  cairo_set_source_rgba (cr, 1.0, 0, 0, 0.7);
  cairo_set_line_width (cr, 10.0);
  cairo_arc(cr, cX(0),cY(0), RAD, 0, 2*M_PI);

  cairo_stroke(cr);
}
void drawUnitCircle()
{
  drawUnitCircleCairo(cairo);
}

/**
 * Init the PDF layout.
 *
 * @param fname
 */
void initPDF(char const *const fname)
{
  CSglobal=cairo_pdf_surface_create(fname, SIZEX+20, SIZEY+20);
  cairo=cairo_create(CSglobal);
}


/**
 * Terminate the PDF file.
 *
 */
void flushPDF()
{
  cairo_show_page(cairo);
  cairo_destroy(cairo);
  cairo_surface_flush(CSglobal);
  cairo_surface_destroy(CSglobal);
}
