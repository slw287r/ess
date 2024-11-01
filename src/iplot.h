#include <stdio.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <ctype.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <inttypes.h>

#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include "utils.h"

typedef struct
{
	int pk;
	double lsd, rsd;
} sd_t;

// plot dimensions
#define MARGIN (35*3)
#define DIM_X (605*2.5)
#define DIM_Y (165*2.5)
#define WIDTH (DIM_X + MARGIN * 2)
#define HEIGHT (DIM_Y + MARGIN * 2)

void draw_box(cairo_t *cr, double x, double y, double width, double height);
void draw_rect(cairo_t *cr, double x, double y, double width, double height);
void draw_rrect(cairo_t *cr);
void draw_arrow(cairo_t *cr, double start_x, double start_y, double end_x, double end_y);
void draw_xlab(cairo_t *cr, const char *xlab);
void draw_ylab(cairo_t *cr, const char *ylab);
void draw_y2lab(cairo_t *cr, const char *ylab);
void draw_xticks(cairo_t *cr, const double xmax);
void draw_yticks(cairo_t *cr, const sp_t *sp);
void draw_is(cairo_t *cr, const int *is, const double *cis, const double pk,
		const int n);
void do_drawing(cairo_t *cr, const int *is, const double *cis, const int n,
		const sd_t *sd, const char *title, const char *sub);
