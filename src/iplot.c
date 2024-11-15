#include "iplot.h"

void draw_box(cairo_t *cr, double x, double y, double width, double height)
{
	cairo_save(cr);
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	cairo_rectangle(cr, x, y, width, height);
	cairo_stroke(cr);
	cairo_restore(cr);
}

void draw_rect(cairo_t *cr, double x, double y, double width, double height)
{
	cairo_save(cr);
	cairo_new_sub_path (cr);
	cairo_rectangle(cr, x, y, width, height);
	cairo_set_source_rgba (cr, .96, .96, .96, 0.86);
	cairo_fill(cr);
	cairo_restore(cr);
}

void draw_rrect(cairo_t *cr)
{
	// a custom shape that could be wrapped in a function
	double x         = 0,        // parameters like cairo_rectangle
	       y         = 0,
	       width         = WIDTH,
	       height        = HEIGHT,
	       aspect        = 1.0,     // aspect ratio
	       corner_radius = height / 60.0;   // and corner curvature radius
	double radius = corner_radius / aspect;
	double degrees = M_PI / 180.0;
	cairo_new_sub_path (cr);
	cairo_arc (cr, x + width - radius, y + radius, radius, -90 * degrees, 0 * degrees);
	cairo_arc (cr, x + width - radius, y + height - radius, radius, 0 * degrees, 90 * degrees);
	cairo_arc (cr, x + radius, y + height - radius, radius, 90 * degrees, 180 * degrees);
	cairo_arc (cr, x + radius, y + radius, radius, 180 * degrees, 270 * degrees);
	cairo_close_path (cr);
	cairo_set_source_rgba (cr, .96, .96, .96, .5);
	cairo_fill(cr);
}

void draw_arrow(cairo_t *cr, double start_x, double start_y, double end_x, double end_y)
{
	double angle = atan2(end_y - start_y, end_x - start_x) + M_PI;
	double arrow_degrees_ = M_PI / 15;
	double arrow_length_ = 10;
	double x1 = end_x + arrow_length_ * cos(angle - arrow_degrees_);
	double y1 = end_y + arrow_length_ * sin(angle - arrow_degrees_);
	double x2 = end_x + arrow_length_ * cos(angle + arrow_degrees_);
	double y2 = end_y + arrow_length_ * sin(angle + arrow_degrees_);
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, start_x, start_y);
	cairo_line_to(cr, (x1 + x2) / 2, (y1 + y2) / 2);
	cairo_stroke(cr);
	cairo_set_line_width(cr, fmin(w1, w2) / 2);
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER);
	cairo_move_to(cr, x1, y1);
	cairo_line_to(cr, x2, y2);
	cairo_line_to(cr, end_x, end_y);
	cairo_line_to(cr, x1, y1);
	cairo_close_path(cr);
	cairo_fill(cr);
}

void draw_xlab(cairo_t *cr, const char *xlab)
{
	double x, y;
	cairo_text_extents_t ext;
	cairo_set_font_size(cr, 18.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = DIM_Y + ext.height * 3;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, xlab);
}

void draw_ylab(cairo_t *cr, const char *ylab)
{
	cairo_save(cr);
	cairo_set_font_size(cr, 18.0);
	cairo_text_extents_t ext;
	cairo_set_source_rgb(cr, 87 / 255.0, 62 / 255.0, 166 / 255.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_translate(cr, MARGIN / 2.0, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, ylab, &ext);
	cairo_move_to(cr, MARGIN / 2.0 + ext.height / 10, -MARGIN * 1.25);
	cairo_show_text(cr, ylab);
	cairo_restore(cr);
}

void draw_y2lab(cairo_t *cr, const char *y2lab)
{
	cairo_save(cr);
	cairo_set_font_size(cr, 18.0);
	cairo_text_extents_t ext;
	cairo_set_source_rgb(cr, 187 / 255.0, 12 / 255.0, 16 / 255.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_translate(cr, MARGIN / 2.0, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, y2lab, &ext);
	cairo_move_to(cr, MARGIN / 2.0, WIDTH - MARGIN * 1.75);
	cairo_show_text(cr, y2lab);
	cairo_restore(cr);
}

void draw_xticks(cairo_t *cr, const double xmax)
{
	double x, m;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	cairo_set_font_size(cr, 16.0);
	const double dashes[] = {0.75, 5.0, 0.75, 5.0};
	int ndash = sizeof(dashes) / sizeof(dashes[0]);
	char buf[sizeof(uint64_t) * 8 + 1];
	cairo_text_extents(cr, "m", &ext);
	double y_offset = ext.height;
	double major = nice_interval(xmax, 15), minor = major / 10;
	for (x = 0; x <= xmax; x += major)
	{
		sprintf(buf, "%g", x);
		cairo_set_source_rgb(cr, 0.16, 0.16, 0.16);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, buf, &ext);
		cairo_move_to(cr, DIM_X * x / xmax - ext.width / 2, DIM_Y + y_offset * 2);
		cairo_show_text(cr, buf);
		// major ticks
		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		cairo_move_to(cr, DIM_X * x / xmax, DIM_Y);
		cairo_line_to(cr, DIM_X * x / xmax, DIM_Y - y_offset);
		// minor ticks
		for (m = minor; m <= minor * 9; m += minor)
		{
			if ((x + m) / xmax > 1)
				break;
			cairo_set_line_width(cr, fmin(w1, w2) / 3.0);
			cairo_move_to(cr, DIM_X * (x + m) / xmax, DIM_Y);
			cairo_line_to(cr, DIM_X * (x + m) / xmax, DIM_Y - y_offset / (m == minor * 5 ? 1.25 : 1.75));
			cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
		}
		cairo_stroke(cr);
		if (x != 0 && x != xmax)
		{
			cairo_set_dash(cr, dashes, ndash, 0);
			cairo_move_to(cr, DIM_X * x / xmax, DIM_Y - y_offset);
			cairo_line_to(cr, DIM_X * x / xmax, 0);
			cairo_stroke(cr);
			cairo_set_dash(cr, dashes, 0, 0);
		}
	}
	cairo_stroke(cr);
}

void draw_yticks(cairo_t *cr, const sp_t *sp, double *scale)
{
	int i;
	double x, y;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	cairo_set_font_size(cr, 16.0);
	const double dashes[] = {0.75, 5.0, 0.75, 5.0};
	int ndash = sizeof(dashes) / sizeof(dashes[0]);
	if (sp->step >= 1e9)
		*scale = 1e9;
	else if (sp->step >= 1e6)
		*scale = 1e6;
	else if (sp->step >= 1e3)
		*scale = 1e3;
	char buf[sizeof(uint64_t) * 8 + 1];
	{
		cairo_text_extents(cr, "m", &ext);
		double x_offset = ext.width;
		// get precision of step
		int p = 0;
		// check precision
		int pc = 0;
		char *pp = NULL;
		for (i = 0; i <= sp->peak / sp->step; ++i)
		{
			sprintf(buf, "%g", i * 100 * sp->step / sp->peak);
			if ((pp = strchr(buf, '.')))
				pc = fmax(pc, strlen(pp + 1));
		}
		for (i = 0; i <= sp->peak / sp->step; ++i)
		{
			sprintf(buf, "%.*f", p, i * sp->step / (*scale));
			cairo_text_extents(cr, buf, &ext);
			x = -ext.width - x_offset / 2.5;
			y = i * sp->step / sp->peak;
			cairo_set_source_rgb(cr, 87 / 255.0, 62 / 255.0, 166 / 255.0);
			cairo_move_to(cr, x, DIM_Y - y * DIM_Y + ext.height / 2);
			cairo_show_text(cr, buf);
			// mirror y
			sprintf(buf, "%.*f", i == sp->peak / sp->step || !i ? 0 : pc, i * 100 * sp->step / sp->peak);
			cairo_text_extents(cr, "m", &ext);
			x = -ext.width / 4;
			y = i * sp->step / sp->peak;
			cairo_set_source_rgb(cr, 187 / 255.0, 12 / 255.0, 16 / 255.0);
			cairo_move_to(cr, DIM_X - x, DIM_Y - y * DIM_Y + ext.height / 2);
			cairo_show_text(cr, buf);
			// major ticks
			cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
			cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
			cairo_move_to(cr, DIM_X, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, DIM_X - x_offset * .75, DIM_Y - y * DIM_Y);
			cairo_stroke(cr);
			if (i != 0 && i != sp->peak / sp->step)
			{
				cairo_set_dash(cr, dashes, ndash, 0);
				cairo_move_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
				cairo_line_to(cr, DIM_X, DIM_Y - y * DIM_Y);
				cairo_stroke(cr);
				// reset dashes
				cairo_set_dash(cr, dashes, 0, 0);
			}
		}
		cairo_stroke(cr);
	}
}

void draw_is(cairo_t *cr, const int *is, const double *cis, const double pk,
		const int n)
{
	int i;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) * 2);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	int xmin = 0, xmax = n;
	for (xmin = 0; xmin <= n; ++xmin)
		if (is[xmin])
			break;
	// isize frequency
	cairo_move_to(cr, (double)xmin / xmax, 1 - (double)is[xmin] / pk);
	cairo_set_source_rgb(cr, 87 / 255.0, 62 / 255.0, 166 / 255.0);
	for (i = xmin + 1; i <= xmax; ++i)
		cairo_line_to(cr, (double)i / xmax, 1 - (double)is[i] / pk);
	cairo_save(cr);
	cairo_scale(cr, 1.0, (double)DIM_X / DIM_Y);
	cairo_stroke(cr);
	cairo_restore(cr);
	// isize cumulative
	cairo_move_to(cr, (double)xmin / xmax, 1 - cis[xmin]);
	cairo_set_source_rgb(cr, 187 / 255.0, 12 / 255.0, 16 / 255.0);
	for (i = xmin + 1; i <= xmax; ++i)
		if (cis[i] >= 0)
			cairo_line_to(cr, (double)i / xmax, 1 - cis[i]);
	cairo_save(cr);
	cairo_scale(cr, 1.0, (double)DIM_X / DIM_Y);
	cairo_stroke(cr);
	cairo_restore(cr);
}

void do_drawing(cairo_t *cr, const int *is, const double *cis, const int n,
		const sd_t *sd, const char *title, const char *sub)
{
	int i = 0;
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
	draw_rrect(cr);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_translate(cr, MARGIN, MARGIN / 1);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	// title
	cairo_set_font_size(cr, 24.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, title, &ext);
	if (!sub)
	{
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.height / 2 + ext.y_bearing * 1.5;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, title);
	}
	else
	{
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = -ext.height * 2;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, title);
		// subtitle
		cairo_set_font_size(cr, 20.0);
		cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_ITALIC, CAIRO_FONT_WEIGHT_BOLD);
		cairo_text_extents(cr, sub, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, sub);
	}
	// zlab
	char zlab[NAME_MAX];
	sprintf(zlab, "Mode: %d SD (-%.0f/+%.0f)", sd->pk, sd->lsd, sd->rsd);
	cairo_set_font_size(cr, 18.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, zlab, &ext);
	if (sd->pk >= n / 2.0)
		x = ext.x_bearing * 5;
	else
		x = DIM_X - ext.width - ext.x_bearing * 5;
	y = ext.height / 2 - ext.y_bearing;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, zlab);
	// draw isize
	int xmax = n, ymax = 0;
	for (i = 0; i <= n; ++i)
		ymax = fmax(ymax, is[i]);
	// axis
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 1.25);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 87 / 255.0, 62 / 255.0, 166 / 255.0);
	cairo_move_to(cr, 0, 0);
	cairo_line_to(cr, 0, DIM_Y); // yaxis
	cairo_set_source_rgb(cr, 0, 0, 0);
	// ticks
	double scale = 1.0f;
	sp_t sp = {0.0f, 0.0f};
	step_and_peak(ymax, &sp);
	draw_xticks(cr, xmax);
	draw_yticks(cr, &sp, &scale);
	draw_box(cr, 0, 0, DIM_X, DIM_Y);
	// xlab
	char xlab[] = "Insert Sizes (bp)";
	draw_xlab(cr, xlab);
	// ylab
	char ylab[NAME_MAX] = {'\0'};
	if (scale == 1)
		strncpy(ylab, "Frequency", NAME_MAX);
	else
		snprintf(ylab, NAME_MAX, "Frequency (%c)", scale == 1e9 ? 'G' :
				(scale == 1e6 ? 'M' : 'K'));
	draw_ylab(cr, ylab);
	// y2lab
	char y2lab[] = "Cumulative (%)";
	draw_y2lab(cr, y2lab);
	cairo_save(cr);
	cairo_scale(cr, DIM_X, DIM_Y);
	draw_is(cr, is, cis, sp.peak, xmax);
	cairo_restore(cr);
}
