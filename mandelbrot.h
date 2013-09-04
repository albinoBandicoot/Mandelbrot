#ifndef MANDELBROT_H
#define MANDELBROT_H

#include <cmplx.h>
#include <color.h>
struct tile {
	int	llx;
	int 	lly;
	int 	xsize;
	int	ysize;
};

struct mand_opts {
	int	xsize;
	int	ysize;
	int	tilesize;	// square tiles, except for edges & corner.
	Complex ll;
	Complex	size;

	// coloring options
	double	exponent;
	double	linear;
	double	constant;

	// edge refinement and iteration count options
	int	iter;	// base iteration count
	int	refinement_passes;

	// supersampling options
	double	qcx_thresh;
	double	grid_thresh;
	double	grid_size;	// in pixels
	int	grid_fac;
	double	grid_wt_tot;
	double	qcx_wt;	// per sample

};

struct render_thread_options {
	int	nthreads;
	int	thr_num;
	int 	tile_num;	// which tile we are currently working on; useful for passing to message printing methods.
	mand_opts*	mopts;
	color*	img;
};

void init_cramp ();

void* render (void *opt);

void render_tile (tile t, mand_opts opts, color* img, render_thread_options rto);

void get_color (color *pix, double v, mand_opts *opts);

double getV (Complex loc, int iters);

int map (double v, mand_opts *opts);

int main (int argc, const char *argv[]);

#endif
