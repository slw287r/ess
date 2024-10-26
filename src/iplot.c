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
	cairo_set_font_size(cr, 12.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = DIM_Y + MARGIN / 2.5 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, xlab);
}

void draw_ylab(cairo_t *cr, const char *ylab)
{
	cairo_save(cr);
	cairo_set_font_size(cr, 12.0);
	cairo_text_extents_t ext;
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_translate(cr, MARGIN / 2.0, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, ylab, &ext);
	cairo_move_to(cr, MARGIN / 2.0 - ext.width / 2, -MARGIN * 1.25);
	cairo_show_text(cr, ylab);
	cairo_restore(cr);
}

void draw_y2lab(cairo_t *cr, const char *ylab)
{
	cairo_save(cr);
	cairo_text_extents_t ext;
	cairo_set_source_rgb(cr, 166 / 255.0, 122 / 255.0, 87 / 255.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_translate(cr, MARGIN / 2, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, M_PI / 2.0); // was 270
	cairo_text_extents(cr, ylab, &ext);
	//cairo_move_to(cr, MARGIN / 2.5 + WIDTH, -MARGIN * 1.25);
	cairo_move_to(cr, -MARGIN, MARGIN * 1.25 - WIDTH);
	cairo_show_text(cr, ylab);
	cairo_restore(cr);
}

void draw_xticks(cairo_t *cr, const double xmax)
{
	int x, m;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	cairo_set_font_size(cr, 10.0);
	const double dashes[] = {0.75, 5.0, 0.75, 5.0};
	int ndash = sizeof(dashes) / sizeof(dashes[0]);
	char buf[sizeof(uint64_t) * 8 + 1];
	cairo_text_extents(cr, "m", &ext);
	double y_offset = ext.height;
	for (x = 0; x <= xmax; x += 50)
	{
		sprintf(buf, "%d", x);
		cairo_set_source_rgb(cr, 0.16, 0.16, 0.16);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, buf, &ext);
		cairo_move_to(cr, DIM_X * x / xmax - ext.width / 2, DIM_Y + y_offset * 2);
		cairo_show_text(cr, buf);
		// major ticks
		cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
		cairo_move_to(cr, DIM_X * x / xmax, DIM_Y);
		cairo_line_to(cr, DIM_X * x / xmax, DIM_Y - y_offset / 2);
		for (m = 10; m <= 40 && m <= xmax - 10; m += 10)
		{
			cairo_set_line_width(cr, fmin(w1, w2) / 3.0);
			cairo_move_to(cr, DIM_X * (x + m) / xmax, DIM_Y);
			cairo_line_to(cr, DIM_X * (x + m) / xmax, DIM_Y - y_offset / 3);
			cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
		}
		cairo_stroke(cr);
		if (x != 0 && x != xmax)
		{
			cairo_set_dash(cr, dashes, ndash, 0);
			cairo_move_to(cr, DIM_X * x / xmax, DIM_Y - y_offset / 2);
			cairo_line_to(cr, DIM_X * x / xmax, 0);
			cairo_stroke(cr);
			cairo_set_dash(cr, dashes, 0, 0);
		}
	}
	cairo_stroke(cr);
}

void draw_yticks(cairo_t *cr, const sp_t *sp, const bool logscale)
{
	int i, j;
	double x, y;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	cairo_set_font_size(cr, 10.0);
	const double dashes[] = {0.75, 5.0, 0.75, 5.0};
	int ndash = sizeof(dashes) / sizeof(dashes[0]);
	char buf[sizeof(uint64_t) * 8 + 1];
	if (logscale)
	{
		double h = ceil(log10(sp->peak));
		cairo_text_extents(cr, "m", &ext);
		double x_offset = ext.width;
		for (i = 0; i <= h; ++i)
		{
			sprintf(buf, "%d", (int)pow(10, h - i));
			cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_text_extents(cr, buf, &ext);
			x = -ext.width - x_offset / 2.5;
			y = 1 - (double)i / (h + 1);
			cairo_move_to(cr, x, DIM_Y - y * DIM_Y + ext.height / 2);
			cairo_show_text(cr, buf);
			// major ticks
			cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
			// minor ticks
			for (j = 2; j <= 9 && i < h; ++j)
			{
				y = (log10((11 - j) * pow(10, i)) + 1)  / (h + 1);
				cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
				cairo_line_to(cr, x_offset * .375, DIM_Y - y * DIM_Y);
			}
		}
		cairo_stroke(cr);
	}
	else
	{
		cairo_text_extents(cr, "m", &ext);
		double x_offset = ext.width;
		for (i = 0; i <= sp->peak / sp->step; ++i)
		{
			sprintf(buf, "%.0f", i * sp->step);
			cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
			cairo_text_extents(cr, buf, &ext);
			x = -ext.width - x_offset / 2.5;
			y = i * sp->step / sp->peak;
			cairo_set_source_rgb(cr, 0.16, 0.16, 0.16);
			cairo_move_to(cr, x, DIM_Y - y * DIM_Y + ext.height / 2);
			cairo_show_text(cr, buf);
			// major ticks
			cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
			cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, x_offset * .5, DIM_Y - y * DIM_Y);
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

void draw_y2ticks(cairo_t *cr, const double ymax)
{
	int i, j;
	double x, y;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	char buf[sizeof(uint64_t) * 8 + 1];
	double h = ceil(log10(ymax));
	cairo_text_extents(cr, "m", &ext);
	double x_offset = ext.width;
	for (i = 0; i <= h; ++i)
	{
		sprintf(buf, "%d", (int)pow(10, h - i));
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, buf, &ext);
		x = -ext.width - x_offset / 2.5;
		y = 1 - (double)i / (h + 1);
		cairo_move_to(cr, x, DIM_Y - y * DIM_Y + ext.height / 2);
		cairo_show_text(cr, buf);
		// major ticks
		cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
		cairo_line_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
		// minor ticks
		for (j = 2; j <= 9 && i < h; ++j)
		{
			y = (log10((11 - j) * pow(10, i)) + 1)  / (h + 1);
			cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, x_offset * .375, DIM_Y - y * DIM_Y);
		}
	}
	cairo_stroke(cr);
}

void draw_is(cairo_t *cr, const int *is, const double peak, const int n)
{
	int i;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	int xmin = 0, xmax = n, ymax = 0;
	for (xmin = 0; xmin <= n; ++xmin)
		if (is[xmin])
			break;
	for (i = 0; i <= n; ++i)
		ymax = (int)fmax(ymax, is[i]);
	double h = pow(10, floor(log10(ymax)));
	h = (int)(ymax / h + 1) * h;
	cairo_move_to(cr, (double)xmin / xmax, 1 - (double)is[xmin] / peak);
	for (i = xmin + 1; i < xmax; ++i)
		cairo_line_to(cr, (double)i / xmax, 1 - (double)is[i] / peak);
	cairo_save(cr);
	cairo_scale(cr, 1.0, (double)DIM_X / DIM_Y);
	cairo_stroke(cr);
	cairo_restore(cr);
}

void do_drawing(cairo_t *cr, const int *is, const int n, const sd_t *sd, const char *sname)
{
	int i = 0;
	cairo_set_source_rgb (cr, 0, 0, 0);
	//cairo_translate(cr, MARGIN / 2, MARGIN / 2.0);
	cairo_translate(cr, MARGIN, MARGIN / 2.0);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	// title
	cairo_set_font_size(cr, 13.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, sname, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = ext.height / 2 + ext.y_bearing * 1.5;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, sname);
	// zlab
	char zlab[NAME_MAX];
	sprintf(zlab, "Mode: %d SD (-%.0f/+%.0f)", sd->pk, sd->lsd, sd->rsd);
	cairo_set_font_size(cr, 10.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, zlab, &ext);
	x = ext.x_bearing * 5;
	y = ext.height / 2 - ext.y_bearing;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, zlab);
	// xlab
	char xlab[] = "Insert Sizes (bp)";
	draw_xlab(cr, xlab);
	// ylab
	char ylab[] = "Frequency";
	draw_ylab(cr, ylab);
	// draw isize
	int xmax = n, ymax = 0;
	for (i = 0; i <= n; ++i)
		ymax = fmax(ymax, is[i]);
	// axis
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 1.25);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, 0, 0);
	cairo_line_to(cr, 0, DIM_Y); // yaxis
	// yticks
	sp_t sp = {0, 0};
	step_and_peak(ymax, &sp);
	draw_xticks(cr, xmax);
	draw_yticks(cr, &sp, false);
	draw_box(cr, 0, 0, DIM_X, DIM_Y);
	cairo_save(cr);
	cairo_scale(cr, DIM_X, DIM_Y);
	draw_is(cr, is, sp.peak, xmax);
	cairo_restore(cr);
}
