/* COMPILE USING:  g++ -std=c++11 -Wextra -o testAnimation `pkg-config --cflags --libs gtk+-3.0` testAnimation.cpp */
#include <gtk/gtk.h>
#include <complex>
//#include <limits>
#include <iostream>
#include <vector>
#include "cairoPoincare++.h"

using namespace std;

#define WINDOW_WIDTH  1000
#define WINDOW_HEIGHT 1000
static void draw_inner_rings(cairo_t*, int, int, vector<Point>*, double&, double&, const double&, const bool = true);

static gboolean draw_cb(GtkWidget *widget, cairo_t *cr, gpointer data)
{
  // R code:
  // five2 <- 2*pi / 5
  // angles <- 1:5
  // for (ii in seq(1, 0.2, -0.05)) {
  //   poincare.plot(list(r=ii, rad=five2), list(r=ii, rad=five2*2),T);
  //   for (i in angles) {
  //     poincare.plot(list(r=ii, rad=five2*i), list(r=ii, rad=five2*(i+1)),F);
  //     poincare.plot(list(r=ii, rad=five2*i), list(r=ii, rad=five2*(i+2)),F);
  //   }
  // }
  static const int points = 9;
  //static const double epsilon = std::numeric_limits<double>::epsilon();
  static const double delta = 2*M_PI / points;
  static double r = 1.0;
  static double r_inc = 0.005;
  static const double r_upper = r_inc;
  static const double r_lower = -(r_inc*2);
  static double rot = M_PI/2;
  //static const double rot_inc = 0.0025;

  cairo_set_source_rgb(cr, 255, 255, 255);
  cairo_paint(cr);
  drawUnitCircle(cr);

  //if (r == 1.0 || r < 0.2) r_inc = r_inc * -1;
  if (r == 1.0 || r < -1.0) r_inc = r_inc * -1;
  //if (rot > 2*M_PI) rot -= 2*M_PI;

  draw_inner_rings(cr, points, 2, nullptr, r, rot, delta, true);
/*
  for (int i = 1; i <= points; ++i) {
    Point a = Point(std::polar(r,delta*i+rot));
    Point b = Point(std::polar(r,delta*(i+1)+rot));
    drawEdgeCairo(cr, &a, &b, true);
    double cx, cy, c_r;
    double R, theta;
    computeCircleParameters(&a, &b, &cx, &cy, &c_r);
    computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    inner_ring.push_back(Point(std::polar(R, theta)));

    //b = Point(std::polar(r,delta*(i+2)+rot));
    //drawEdgeCairo(cr, &a, &b, true);
    //computeCircleParameters(&a, &b, &cx, &cy, &c_r);
    //computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    //drawPointCairo(cr, new Point(std::polar(R, theta)));

    //b = Point(std::polar(r,delta*(i+3)+rot));
    //drawEdgeCairo(cr, &a, &b, true);
    //computeCircleParameters(&a, &b, &cx, &cy, &c_r);
    //computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    //drawPointCairo(cr, new Point(std::polar(R, theta)));
    //drawLine(cr, Point(std::polar(r,delta*i+rot)), Point(std::polar(r,delta*(i+4)+rot)));
    //drawLine(cr, Point(std::polar(r,delta*i+rot)), Point(std::polar(r,delta*(i+5)+rot)));
    //drawLine(cr, Point(std::polar(r,delta*i+rot)), Point(std::polar(r,delta*(i+6)+rot)));
  }

  vector<Point> inner_ring2 = vector<Point>();

  for (vector<Point>::const_iterator i = inner_ring.begin(); i != inner_ring.end(); ++i) {
    const Point *a = &*i;
    const Point *b;
    if (i+1 != inner_ring.end())
      b = &*(i + 1);
    else
      b = &*(inner_ring.begin());
    drawEdgeCairo(cr, a, b, true);

    double cx, cy, c_r;
    double R, theta;
    computeCircleParameters(a, b, &cx, &cy, &c_r);
    computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    inner_ring2.push_back(Point(std::polar(R, theta)));
  }

  for (vector<Point>::const_iterator i = inner_ring2.begin(); i != inner_ring2.end(); ++i) {
    const Point *a = &*i;
    const Point *b;
    if (i+1 != inner_ring2.end())
      b = &*(i + 1);
    else
      b = &*(inner_ring2.begin());
    drawEdgeCairo(cr, a, b, true);
  }
*/
//std::cout<<"epsilon="<<epsilon;
//std::cout<<" r_inc="<<r_inc<<" r="<<r<<" diff="<<fabs(r-0.2)<<" bool="<<(bool)(fabs(r-0.2)>fabs(r_inc))<<" rot="<<rot;
//std::cout<<" abs(r)="<<fabs(r);
//std::cout<<std::endl;

  if (r > r_lower && r <= r_upper) r += r_inc; // skip r ~= 0
  r += r_inc;
  //rot += rot_inc;

  return FALSE;
}

static void draw_inner_rings(cairo_t* cr, int points, int count, vector<Point> *inner_ring, double &r, double &rot, const double &delta, const bool lines)
{
  vector<Point> inner_ring2 = vector<Point>();
  const Point *a, *b;
  double cx, cy, c_r;
  double R, theta;

  if (inner_ring == nullptr || inner_ring->empty())
  for (int i = 1; i <= points; ++i) {
    a = new Point(std::polar(r,delta*i+rot));
    b = new Point(std::polar(r,delta*(i+1)+rot));
    if (!lines)
      drawEdgeCairo(cr, a, b, true);
    else
      drawLineCairo(cr, a, b, false);

    computeCircleParameters(a, b, &cx, &cy, &c_r);
    computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    const Point m = Point(std::polar(R, theta));
    drawPointCairo(cr, &m);
    Point polarX, idealX, polarY, idealY;
    computeOmegaPointsFromMidpoint(&m, &polarX, &polarY, &idealX, &idealY);
    double theta1, theta2;
    computeOmegaPoints(&cx, &cy, &c_r, &theta1, &theta2);
    idealX = Point(std::polar(1., theta1));
    idealY = Point(std::polar(1., theta2));
cout<<std::abs(std::complex<double>(polarX.X, polarX.Y))<<","
    <<std::fmod(std::arg(std::complex<double>(polarX.X, polarX.Y))+2*M_PI, 2*M_PI)<<"  "
    <<std::abs(std::complex<double>(polarY.X, polarY.Y))<<","
    <<std::fmod(std::arg(std::complex<double>(polarY.X, polarY.Y))+2*M_PI, 2*M_PI)
    <<endl
    <<std::abs(std::complex<double>(idealX.X, idealX.Y))<<","
    <<std::fmod(std::arg(std::complex<double>(idealX.X, idealX.Y))+2*M_PI, 2*M_PI)<<"  "
    <<std::abs(std::complex<double>(idealY.X, idealY.Y))<<","
    <<std::fmod(std::arg(std::complex<double>(idealY.X, idealY.Y))+2*M_PI, 2*M_PI)
    <<endl;
      cairo_set_source_rgba (cr, 0, 0.6, 0, 0.5);
      cairo_set_line_width (cr, 0.5);
      //We project a and b points onto the unit circle.
      //double theta = atan2(polarX.Y, polarX.X);
      //double theta2 = atan2(idealX.Y, idealX.X);
      //cairo_move_to (cr, cX(cos(theta)),cY(sin(theta)));
      //cairo_line_to (cr, cX(cos(theta2)),cY(sin(theta2)));
      //cairo_stroke(cr);

      //theta = atan2(polarY.Y, polarY.X);
      //theta2 = atan2(idealY.Y, idealY.X);
      //cairo_move_to (cr, cX(cos(theta)),cY(sin(theta)));
      //cairo_line_to (cr, cX(cos(theta2)),cY(sin(theta2)));
      //cairo_stroke(cr);
      
      cairo_move_to (cr, cX(polarX.X),cY(polarX.Y));
      cairo_line_to (cr, cX(idealX.X),cY(idealX.Y));
      cairo_stroke(cr);
      cairo_move_to (cr, cX(polarY.X),cY(polarY.Y));
      cairo_line_to (cr, cX(idealY.X),cY(idealY.Y));
      cairo_stroke(cr);

    inner_ring2.push_back(Point(std::polar(R, theta)));
  }

  else
  for (vector<Point>::const_iterator i = inner_ring->begin(); i != inner_ring->end(); ++i) {
    a = &*i;
    if (i+1 != inner_ring->end())
      b = &*(i + 1);
    else
      b = &*(inner_ring->begin());
    if (!lines)
      drawEdgeCairo(cr, a, b, true);
    else
      drawLineCairo(cr, a, b, false);

    computeCircleParameters(a, b, &cx, &cy, &c_r);
    computeMidpoint(&cx, &cy, &c_r, &R, &theta);
    inner_ring2.push_back(Point(std::polar(R, theta)));
  }

  if (count != 0)
    draw_inner_rings(cr, inner_ring2.size(), --count, &inner_ring2, r, rot, delta, lines);

}

static gboolean draw_trigger(GtkWidget *widget)
{
  if (widget == NULL) return FALSE;

  gtk_widget_queue_draw(widget);

  return TRUE;
}

int main (int argc, char *argv[])
{
   GtkWidget *window;
   GtkWidget *da;

   gtk_init (&argc, &argv);

   window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
   g_signal_connect (window, "destroy", G_CALLBACK (gtk_main_quit), NULL);

   da = gtk_drawing_area_new();
   gtk_widget_set_size_request (da, WINDOW_WIDTH, WINDOW_HEIGHT);
   g_signal_connect (da, "draw", G_CALLBACK(draw_cb),  NULL);

   gtk_container_add (GTK_CONTAINER (window), da);
   gtk_widget_show (da);
   gtk_widget_show (window);

   g_timeout_add (34, (GSourceFunc)draw_trigger, (gpointer)da);

   gtk_main ();

   return 0;
}
