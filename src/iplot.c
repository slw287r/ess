#include "iplot.h"

void draw_box(cairo_t *cr, double x, double y, double width, double height)
{
	cairo_save(cr);
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
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
	cairo_set_font_size(cr, 10.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, xlab);
}

void draw_ylab(cairo_t *cr, const char *ylab)
{
	cairo_save(cr);
	cairo_text_extents_t ext;
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	//cairo_translate(cr, MARGIN / 1.25, HEIGHT / 2.0); // translate origin to the center
	cairo_translate(cr, MARGIN / 2.0, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, ylab, &ext);
	//cairo_move_to(cr, MARGIN / 5, -MARGIN);
	cairo_move_to(cr, MARGIN / 5.0, -MARGIN / 1.5);
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

void draw_yticks(cairo_t *cr, const double ymax, const bool logscale)
{
	int i, j;
	double x, y;
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	char buf[sizeof(uint64_t) * 8 + 1];
	if (logscale)
	{
		double h = ceil(log10(ymax));
		cairo_text_extents(cr, "m", &ext);
		double x_offset = ext.width;
		for (i = 0; i <= h; ++i)
		{
			sprintf(buf, "%d", (int)pow(10, h - i));
			cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
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
	}
	else
	{
		; // TODO
	}
	cairo_stroke(cr);
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

void draw_is(cairo_t *cr, const int *is, const int n)
{
	int i;
	double w1 = 1.0, w2 = 1.0, x = 0, y = 0, x1, y1;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	// get max cpu usage
	for (i = n; i > 0; --i)
		if (is[i])
			break;
	int xmax = i, ymax = 0;
	for (i = xmax; i >= 0; --i)
		ymax = fmax(ymax, is[i]);
	if (xmax)
	{
		// draw cpu history
		for (i = 0; i <= xmax; ++i)
		{
			x1 = (double)(i + 1) / (xmax + 1);
			y1 = (double)is[i] / ymax;
			cairo_move_to(cr, x, 1 - y);
			cairo_line_to(cr, x1, 1 - y1);
			x = x1;
			y = y1;
		}
		cairo_save(cr);
		cairo_scale(cr, 1.0, (double)DIM_X / DIM_Y);
		cairo_stroke(cr);
		cairo_restore(cr);
	}
}

void do_drawing(cairo_t *cr, const int *is, const int n, const sd_t *sd)
{
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_translate(cr, MARGIN / 2, MARGIN / 2.0);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	//cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
	char xlab[] = "Insert size (bp)";
	char ylab[] = "Frequency";
	// title
	char title[] = "Insert size distribution";
	cairo_set_font_size(cr, 10.0);
	cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, title, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = ext.height / 2 + ext.y_bearing * 1.5;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, title);
	// zlab
	char zlab[NAME_MAX];
	sprintf(zlab, "Mode: %d SD (-%.0f/+%.0f)", sd->pk, sd->lsd, sd->rsd);
	cairo_set_font_size(cr, 10.0);
	cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, zlab, &ext);
	x = ext.x_bearing * 35;
	y = ext.height / 2 - ext.y_bearing;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, zlab);
	// xlab
	draw_xlab(cr, xlab);
	// ylab
	draw_ylab(cr, ylab);
	// get max cpu and mem
	int i = 0;
	draw_box(cr, 0, 0, DIM_X, DIM_Y);
	// draw isize
	cairo_save(cr);
	cairo_scale(cr, DIM_X, DIM_Y);
	draw_is(cr, is, n);
	cairo_restore(cr);
	// axis
	//draw_arrow(cr, 0, DIM_Y*1.005, DIM_X, DIM_Y*1.005); // xaxis
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 1.25);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, 0, 0);
	cairo_line_to(cr, 0, DIM_Y); // yaxis
	// yticks
	int ymax = 0;
	for (i = 0; i < n; ++i)
		ymax = fmax(ymax, is[i]);
	draw_yticks(cr, ymax, false);
	// y2ticks
	//draw_y2ticks(cr, mem_max);
	char *a = NULL;
	/* draw runtime at bottom right
	stoa(mns[n - 1]->ts - mns[0]->ts, &a);
	asprintf(&a, "Runtime: %s", a);
	cairo_set_font_size(cr, 8.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, a, &ext);
	x = DIM_X - ext.width - ext.x_bearing;
	y = DIM_Y + ext.height * 3 + ext.y_bearing; // bottom right
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, a);
	*/
	/* start time
	asprintf(&a, "Start: %s", st);
	cairo_text_extents(cr, a, &ext);
	x = ext.x_bearing;
	y = DIM_Y + ext.height * 3 + ext.y_bearing; // bottom left
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, a);
	*/
	// legend
	/* max cpu usage
	asprintf(&a, "%.*f", cpu_max < 1 ? 2 : 0, cpu_max);
	cairo_text_extents(cr, a, &ext);
	cairo_set_font_size(cr, 8.0);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	x = -ext.x_bearing * 4 - ext.width;
	y = ext.height;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, a);
	// max memory
	if (mem_max < 1)
	{
		double mem_max_m = mem_max * 1000;
		if (mem_max_m < 1)
			asprintf(&a, "<1M");
		else
			asprintf(&a, "%.0fM", mem_max_m);
	}
	else
		asprintf(&a, "%.0fG", mem_max);
	cairo_text_extents(cr, a, &ext);
	cairo_set_source_rgb(cr, 166 / 255.0, 122 / 255.0, 87 / 255.0);
	//x = DIM_X - ext.width + ext.x_bearing;
	x = DIM_X + ext.x_bearing;
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, a);
	*/
	/* draw legend
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	asprintf(&a, "—CPU");
	cairo_text_extents(cr, a, &ext);
	x = DIM_X - ext.width * 1.25;
	y = ext.height - ext.y_bearing * 1.5;
	cairo_move_to(cr, x, y);
	cairo_text_path(cr, a);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_fill_preserve (cr);
	cairo_set_source_rgb (cr, 1, 1, 1); // white border
	cairo_set_line_width (cr, .35);
	cairo_stroke (cr);
	
	asprintf(&a, "—RSS");
	x = DIM_X - ext.width * 1.25;
	y = ext.height - ext.y_bearing * 3.0;
	cairo_move_to(cr, x, y);
	cairo_text_path(cr, a);
	cairo_set_source_rgb(cr, 166 / 255.0, 122 / 255.0, 87 / 255.0);
	cairo_fill_preserve (cr);
	cairo_set_source_rgb (cr, 1, 1, 1); // white border
	cairo_set_line_width (cr, .35);
	cairo_stroke (cr);

	asprintf(&a, "—SHR");
	x = DIM_X - ext.width * 1.25;
	y = ext.height - ext.y_bearing * 4.5;
	cairo_move_to(cr, x, y);
	cairo_text_path(cr, a);
	cairo_set_source_rgb(cr, 218 / 255.0, 165 / 255.0, 32 / 255.0);
	cairo_fill_preserve (cr);
	cairo_set_source_rgb (cr, 1, 1, 1); // white border
	cairo_set_line_width (cr, .35);
	cairo_stroke (cr);
	*/
	if (a)
		free(a);
}
