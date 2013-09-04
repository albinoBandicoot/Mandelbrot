#include <color.h>
#include <image.h>
#include <cmplx.h>
#include <mandelbrot.h>
#include <pthread.h>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <cstring>
using namespace std;

/* First we define the color ramp to use. We start with a small number of control points
 * and then interpolate to fill out a list of RPOINTS colors. This way we just look up in
 * the list instead of compute the colors on the fly.
 */

#define	CPOINTS	10
#define	RPOINTS	2048

color *ramp;	// the main ramp
color *control;	// the control points. For now, they are equally spaced.

unsigned long ICT = 0;
int QCX_TOT_PX = 0;
int GRID_TOT_PX = 0;
int ncols = 80;

void fill_line (int len){
	for (int i=len; i<ncols; i++){
		printf(" ");
	}
	fflush(stdout);
}

void clear_line (){
	for (int i=0; i<ncols; i++){
		putc((char) 8, stdout);
	}
	fflush(stdout);
}

void init_cramp (){	// call this at the beginning to set up the color ramp.
	control = (color *)(malloc(CPOINTS * sizeof (color)));
	ramp    = (color *)(malloc(RPOINTS * sizeof (color)));
	/* initialize the colors to something nice */
	col_init_uc(&control[0], 235, 255, 255);
	col_init_uc(&control[1], 250, 242, 166);
	col_init_uc(&control[2], 247, 168,   2);
	col_init_uc(&control[3], 169,  64,  14);
	col_init_uc(&control[4], 84,  11,  35);
	col_init_uc(&control[5], 27,   2,  58);
	col_init_uc(&control[6], 7,  38, 153);
	col_init_uc(&control[7], 33, 104, 195);
	col_init_uc(&control[8], 78, 160, 227);
	col_init_uc(&control[9], 171, 228, 251);

	/* now fill out the main ramp by interpolation */
	for (int i=0; i<RPOINTS; i++){
		int pcol = (int) (((double) CPOINTS*i)/RPOINTS);	// the index of the color to the 'left'
		int ncol = (pcol+1)%CPOINTS;				// the index of the color to the 'right'
		double d = (((double) (CPOINTS*i))/RPOINTS)-pcol;	// the blending factor
		color temp1, temp2;
		col_mul(&temp1, &control[pcol], 1-d);
		col_mul(&temp2, &control[ncol], d);
		col_add(&ramp[i], &temp1, &temp2);
	}
}

int min (int a, int b){
	if (a < b) return a;
	return b;
}

int max (int a, int b){
	if (a > b) return a;
	return b;
}

/* Each pthread will call its own version of render. It is passed a render_thread_options*,
 * from which it will pull all of its options. Included in this is the thread number and the
 * # of render threads. It will then loop over tiles, selecting each one whose index is 
 * congruent to thr_num % nthreads. Then it renders it.
 */

pthread_mutex_t	print_mutex;
char **messages;

void print_messages (int nthreads){
	pthread_mutex_lock(&print_mutex);
	for (int i=0; i<nthreads; i++){
		clear_line();
	}
	for (int i=0; i<nthreads; i++){
		int len = printf("%s", messages[i]);
		fill_line(len);
	}
	pthread_mutex_unlock(&print_mutex);
}

void message (render_thread_options rto, char *str){
	sprintf(messages[rto.thr_num], "Thread %i working on tile #%i: %s", rto.thr_num, rto.tile_num, str);
	print_messages(rto.nthreads);
}

pthread_mutex_t tile_manager;
int tllx = 0;
int tlly = 0;
bool DONE = false;

tile get_tile (mand_opts opts){
	pthread_mutex_lock(&tile_manager);
	tile t;
	if (DONE){
		t.llx = -1;
		t.xsize = 0;
		t.ysize = 0;
		pthread_mutex_unlock(&tile_manager);
		return t;
	}
	t.llx = tllx;
	t.lly = tlly;
	t.xsize = min(opts.tilesize, opts.xsize-t.llx);
	t.ysize = min(opts.tilesize, opts.ysize-t.lly);
	tllx += opts.tilesize;
	if (tllx >= opts.xsize){
		tllx = 0;
		tlly += opts.tilesize;
		DONE = tlly >= opts.ysize;
	}
	pthread_mutex_unlock(&tile_manager);
	return t;
}

void* render (void *options){	// this is the method that we use in the pthread creation
	/* extract options */
	render_thread_options rto = *((render_thread_options*) options);
	mand_opts opts = *(rto.mopts);
	while (true){
		tile t = get_tile(opts);
		rto.tile_num++;
		if (t.llx == -1){
			sprintf(messages[rto.thr_num], "Thread %i has stopped.", rto.thr_num);
			print_messages(rto.nthreads);
			return NULL;
		}
		render_tile(t, opts, rto.img, rto);
	}
	/* set up the tile */
	/*
	tile t;
	t.llx = 0;
	t.lly = 0;
	int tnum = 0;
	// The main loop - looping over the possible tiles 
//	printf("Inside the thread, opts.xsize = %i, opts.ysize = %i\n", opts.xsize, opts.ysize);
	while (t.llx < opts.xsize){
		while (t.lly < opts.ysize){
//			printf("THREAD %i; looking at llx = %i; lly = %i.\n", rto.thr_num, t.llx, t.lly);
			if (tnum%rto.nthreads == rto.thr_num){	// it's ours
				// We have to be a little careful on the size computation with the border tiles 
				t.xsize = min(opts.tilesize, opts.xsize-t.llx);
				t.ysize = min(opts.tilesize, opts.ysize-t.lly);
//				printf("\tTHREAD %i; will render llx = %i; lly = %i.\n", rto.thr_num, t.llx, t.lly);
				render_tile (t, opts, rto.img, rto);
				rto.tile_num++;
//				printf("\t after render, llx = %i, lly = %i, tilesize = %i\n", t.llx, t.lly, opts.tilesize);
			}
			tnum++;
			t.lly += opts.tilesize;
		}
		t.lly = 0;
		t.llx += opts.tilesize;
	}
	*/
}

void clear (bool *ptr, int size){
	for (int i=0; i<size; i++){
		ptr[i] = false;
	}
}

void fill (char *ptr, int size, char val){
	for (int i=0; i<size; i++){
		ptr[i] = val;
	}
}

double *filter_kernel;
double kernel_sum;
const double B = 1;	// B and C are the two parameters for the family of cubic splines
const double C = 0;
const double FK_BIAS = 0.5;

double kernel (double x){
	x = abs(x);
	if (x >= 2) return 0;
	if (x >= 1){
		return (-B-6*C)*x*x*x + (6*B+30*C)*x*x + (-12*B-48*C)*x + (8*B+24*C);
	}
	return (12-9*B-6*C)*x*x*x + (-18+12*B+6*C)*x*x + (6-2*B);
}

double kernel_avg (double lo, double hi){	// compute an approximation of the average value of the kernel over the interval [lo..hi]
	double step = (hi-lo)/500;
	double pos = lo;
	double res = 0;
	for (int i=0; i<500; i++){
		res += kernel(pos);
		pos += step;
	}
	return res/500;
}

void compute_kernel (int size){
	printf("Filtering kernel:\n");
	filter_kernel = (double *)(malloc(size*size*sizeof(double)));
	double k_1d[size];
	double ppos = -2;
	for (int i=0; i<size; i++){
		double x = 4*((i+1.0)/(size) - 0.5);	// map to [-2..2]
		printf("k_1d[%i] = avg(%f, %f)\n", i, ppos, x);
		k_1d[i] = kernel_avg(ppos, x) + FK_BIAS;
		ppos = x;
	}
	for (int i=0; i<size; i++){
		for (int j=0; j<size; j++){
			filter_kernel[i*size+j] = k_1d[i] * k_1d[j];
			kernel_sum += filter_kernel[i*size+j];
		}
	}
	for (int i=0; i<size; i++){
		for (int j=0; j<size; j++){
			filter_kernel[i*size+j] /= kernel_sum;
			printf("%.4f\t", filter_kernel[i*size+j]);
		}
		printf("\n");
	}
	kernel_sum = 1;
}



#define EDGE4	// it turns out 4-neighbor edges is more efficient than 8-neighbor edges (4% speedup on one example)

/* The guts: loops over pixels, and does any edge refinement and supersampling. */
void render_tile (tile t, mand_opts opts, color* img, render_thread_options rto){
	Complex pxstep = {opts.size.a/opts.xsize, opts.size.b/opts.ysize};	// the size of a pixel
	double 	icts[t.xsize+2][t.ysize+2];	// the border is for the edge refinement stuff.
	message(rto, "initial set detection");
	for (int x=t.llx-1; x<t.llx+t.xsize+1; x++){
		for (int y=t.lly-1; y<t.lly+t.ysize+1; y++){
			if (x < t.llx || x > t.llx+t.xsize || y < t.lly || y > t.lly+t.ysize){
//				icts[x-t.llx+1][y-t.lly+1] = 0;	// anything but -1.
			} else {
				Complex px = {opts.ll.a + pxstep.a*x, opts.ll.b + pxstep.b*y};	// coords of the current pixel
				icts[x-t.llx+1][y-t.lly+1] = getV(px, opts.iter);
			}
		}
	}
	for (int x=0; x<t.xsize+2; x++){
		icts[x][0] = max(0,icts[x][1]);	// we use max(0, ...) here so we never set it to -1.
		icts[x][t.ysize+1] = max(0,icts[x][t.ysize]);
	}
	for (int y=0; y<t.ysize+2; y++){
		icts[0][y] = max(0,icts[1][y]);
		icts[t.xsize+1][y] = max(0,icts[t.xsize][y]);
	}

	/* now we proceed with edge refinement. The idea is that adding iterations shaves off pixels
	 * from the edge (boundary) of the set. So, we look only at the in-set pixels that border
	 * out-of-set pixels. Then we repeat until only a few pixels are being taken off. Then we
	 * increase the iterations again. 
	 *
	 * we keep a taboo list of pixels already searched at the current iteration count, so we don't
	 * waste time re-evaluating them later.
	 */

	char str[60];

	bool taboo[t.xsize][t.ysize];
	char itersave[t.xsize+2][t.ysize+2];	// with a border.
	fill (&itersave[0][0], (t.xsize+2)*(t.ysize+2), 0);

	const int REMCT_THRESH = 1;
	int iterations = opts.iter;
	int tot_remct = REMCT_THRESH+1;
	for (int i=0; i<opts.refinement_passes; i++){
		tot_remct = 0;
		sprintf(str, "refining edge: pass %i with %i iterations...", i+1, iterations*2);
		message(rto, str);
//		printf("Starting refinement pass %i with %i iterations...", i+1, iterations*2);
		clear(&taboo[0][0], t.xsize*t.ysize);
		iterations *= 2;
		int remct = REMCT_THRESH+1;	// the count of the number of pixels removed
		int passct = 0;
		while (remct > REMCT_THRESH){
			remct = 0;
			// now we loop over the pixels.
			for (int x=1; x<t.xsize+1; x++){
				for (int y=1; y<t.ysize+1; y++){
					if (icts[x][y] == -1 && !taboo[x-1][y-1]){	// if we're in the set and not taboo'd, check if we're on the edge.
						// now the question - is it more efficient to define the edge by 4-way or 8-way surrouning checks?
#ifdef EDGE4
						bool on_edge = (icts[x][y-1] + icts[x][y+1] + icts[x-1][y] + icts[x+1][y]) != -4;
#else
						bool on_edge = (icts[x][y-1] + icts[x][y+1] + icts[x-1][y] + icts[x+1][y] + icts[x-1][y-1] + icts[x-1][y+1] + icts[x+1][y-1] + icts[x+1][y+1]) != -8;
#endif
						if (on_edge){
							/* note that modifying icts will influence the evaluation of subsequent pixels in this pass. It is perfectly correct, just different than the way it was before. */
							Complex px = {opts.ll.a + pxstep.a*(x-1+t.llx), opts.ll.b + pxstep.b*(y-1+t.lly)};
							icts[x][y] = getV(px, iterations);
							if (icts[x][y] != -1) remct ++;	// we removed a pixel.
							taboo[x-1][y-1] = true;	// and now we set taboo so we don't do this work again.
							itersave[x][y] = i+1;
						}
					}
				}
			}
			passct++;
			tot_remct += remct;
		}
//		printf("removed %i pixels in %i passes\n", tot_remct, passct);
	}
supersample:
	message(rto, "initial coloring");
	char ssfacs[t.xsize+2][t.ysize+2];
	float divs[t.xsize][t.ysize];
	color* tempimg = (color *)(malloc((t.xsize+2)*(t.ysize+2)*sizeof(color)));	// this also has a border

	/* we re-set the border of icts so that we allow -1, since we want to use this for SS detection now. */
	for (int x=0; x<t.xsize+2; x++){
		icts[x][0] = icts[x][1];	
		icts[x][t.ysize+1] = icts[x][t.ysize];
	}
	for (int y=0; y<t.ysize+2; y++){
		icts[0][y] = icts[1][y];
		icts[t.xsize+1][y] = icts[t.xsize][y];
	}

	// first we do an initial coloring of the image.
//	printf("starting initial coloring...\n");
	for (int x=0; x<t.xsize+2; x++){
		for (int y=0; y<t.ysize+2; y++){
			get_color(&tempimg[x*(t.ysize+2)+y], icts[x][y], &opts);	// note that we color the border from icts
			if (x > 0 && x < t.xsize+1 && y > 0 && y < t.ysize+1){	// and we copy the colors to the final image
				col_copy(&img[(x+t.llx-1)*opts.ysize + (y+t.lly-1)], &tempimg[x*(t.ysize+2)+y]);
			}
		}
	}
	/* now we make a pass to evaluate the pixels to see to what degree SS is necessary there. We compare 
	 * computed float values with the thresholds in opts, to get either 0 = no SS, 1 = QCX, 2 = full grid.
	 */
//	printf("evaluating SS thresholds...\n");
	message(rto, "selecting pixels for supersampling");
	int qcx_tot = 0;
	int grid_tot = 0;
	for (int x=1; x<t.xsize+1; x++){
		for (int y=1; y<t.ysize+1; y++){
			double cdiff = 0.0;
			divs[x-1][y-1] = 1.0;	// for the center sample
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x-1)*(t.ysize+2)+(y-1)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x-1)*(t.ysize+2)+(y)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x-1)*(t.ysize+2)+(y+1)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x)  *(t.ysize+2)+(y-1)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x)  *(t.ysize+2)+(y+1)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x+1)*(t.ysize+2)+(y+1)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x+1)*(t.ysize+2)+(y)]);
			cdiff += col_diff (tempimg[x*(t.ysize+2)+y], tempimg[(x+1)*(t.ysize+2)+(y+1)]);
			if (cdiff > opts.grid_thresh){
				ssfacs[x][y] = 2;
				grid_tot ++;
			} else if (cdiff > opts.qcx_thresh){
				ssfacs[x][y] = 1;
				qcx_tot ++;
			} else {
				ssfacs[x][y] = 0;
			}
		}
	}
	/* Finally we are ready to start supersampling. First we do the QCX pass. */
	sprintf(str, "supsersampling... quincunx: %i pixels marked", qcx_tot);
	message(rto, str);
	int qcx_done = 0;
	color temp;
	for (int x=t.llx; x<t.llx+t.xsize; x++){
		for (int y=t.lly; y<t.lly+t.ysize; y++){
			if (ssfacs[x-t.llx+1][y-t.lly+1] == 1 || ssfacs[x-t.llx][y-t.lly+1] == 1 || ssfacs[x-t.llx][y-t.lly] == 1 || ssfacs[x-t.llx+1][y-t.lly] == 1){	// then this sample needs to be evaluated. 
				Complex px = {opts.ll.a + pxstep.a*(x-0.5), opts.ll.b + pxstep.b*(y-0.5)};
				int iters = max(max(itersave[x-t.llx+1][x-t.lly+1], itersave[x-t.llx+2][y-t.lly+1]), max(itersave[x-t.llx+1][y-t.lly+2], itersave[x-t.llx+2][y-t.lly+2]));
				iters = (int) (opts.iter * pow(2.0, iters));
				get_color(&temp, getV(px, iters), &opts);
				col_mul(&temp, &temp, opts.qcx_wt);
				col_add(&img[x*opts.ysize+y], &img[x*opts.ysize+y], &temp);
				divs[x-t.llx][y-t.lly] += opts.qcx_wt;
				// since we have the samples, why not give them to everybody whether they requested it or not?
				if (x > t.llx){
					col_add(&img[(x-1)*opts.ysize+y], &img[(x-1)*opts.ysize+y], &temp);
					divs[x-t.llx-1][y-t.lly] += opts.qcx_wt;
				}
				if (y > t.lly){
					col_add(&img[x*opts.ysize+y-1], &img[x*opts.ysize+y-1], &temp);
					divs[x-t.llx][y-t.lly-1] += opts.qcx_wt;
				}
				if (x > t.llx && y > t.lly){
					col_add(&img[(x-1)*opts.ysize+y-1], &img[(x-1)*opts.ysize+y-1], &temp);
					divs[x-t.llx-1][y-t.lly-1] += opts.qcx_wt;
				}
				qcx_done++;
			}
		}
		if ((x-t.llx)%(t.xsize/4) == 0){
			sprintf(str, "supersampling... quincunx: %i samples; %i pixels marked", qcx_done, qcx_tot);
		}
	}
	QCX_TOT_PX += qcx_tot;
	
	/* Now we do the full grid pass. No attempt is made to share samples here. */

	sprintf(str, "supersampling... 0 of %i pixels", grid_tot);
	message(rto, str);
	int grid_done = 0;
	for (int x=t.llx; x<t.llx+t.xsize; x++){
		for (int y=t.lly; y<t.lly+t.ysize; y++){
			if (ssfacs[x-t.llx+1][y-t.lly+1] == 2){
				grid_done++;
				Complex cent = {opts.ll.a + pxstep.a*x, opts.ll.b + pxstep.b*y};	// coords of the current pixel center
				double gxsize = pxstep.a*opts.grid_size;
				double gysize = pxstep.b*opts.grid_size;
				double xstep  = gxsize/opts.grid_fac;
				double ystep  = gysize/opts.grid_fac;
//				printf("(x,y) = (%i, %i); cent = (%f, %f); gxsize = %f; gysize = %f, xstep = %f, ystep = %f\n", x, y, cent.a, cent.b, gxsize, gysize, xstep, ystep);
				double a = cent.a-gxsize/2;
				double b = cent.b-gysize/2;
				col_init_fl(&img[x*opts.ysize+y], 0.0f, 0.0f, 0.0f);	// clear the color
				for (int i=0; i<opts.grid_fac; i++){
					for (int j=0; j<opts.grid_fac; j++){
//						printf("a = %f, b = %f, grid fac = %i\n", a, b, opts.grid_fac);
						Complex px = {a,b};
						int iters = itersave[x-t.llx+1][y-t.lly+1];
						iters = max(iters, max(itersave[x-t.llx+1][y-t.lly], itersave[x-t.llx+1][y-t.lly+2]));
						iters = max(max(iters, itersave[x-t.llx][y-t.lly+1]), max(itersave[x-t.llx][y-t.lly], itersave[x-t.llx][y-t.lly+2]));
						iters = max(max(iters, itersave[x-t.llx+2][y-t.lly+1]), max(itersave[x-t.llx+2][y-t.lly], itersave[x-t.llx+2][y-t.lly+2]));
						iters = (int) (opts.iter * pow(2.0, iters));
						get_color(&temp, getV(px, iters), &opts);
//						col_init_uc(&temp, 0, 255, 0);
//						col_mul(&temp, &temp, opts.grid_wt_tot/(opts.grid_size*opts.grid_size));
						col_mul(&temp, &temp, filter_kernel[i*opts.grid_fac+j]);
						col_add(&img[x*opts.ysize+y], &img[x*opts.ysize+y], &temp);
//						divs[x-t.llx][y-t.lly] += opts.grid_wt_tot/(opts.grid_size*opts.grid_size);
						b += ystep;
					}
					b = cent.b-gysize/2;
					a += xstep;
				}
				col_clamp(&img[x*opts.ysize+y]);
			}
		}
		if ((x-t.llx)%(t.xsize/6) == 0){
			sprintf(str, "supersampling... %i of %i pixels", grid_done, grid_tot);
			message(rto, str);
		}
	}
	GRID_TOT_PX += grid_tot;
//	printf("dividing...\n");
	/* now we have to do the divisions */
	for (int x=t.llx; x<t.llx+t.xsize; x++){
		for (int y=t.lly; y<t.lly+t.ysize; y++){
			col_mul(&img[x*opts.ysize+y], &img[x*opts.ysize+y], 1.0/divs[x-t.llx][y-t.lly]);
		}
	}
	free(tempimg);
//	printf("done\n");
}

void get_color (color *col, double v, mand_opts *opts){	// given the result of getV and options, sticks the color in col.
	int pos = map (v, opts);
	if (pos == -1){
		col_init_fl(col, 0.0f, 0.0f, 0.0f);
	} else {
		col_copy(col, &ramp[pos]);
	}
}

const double log_2 = log(2.0);
const double bailout_sqr = 1e80;
const double log_b = log(1e40);

double getV (Complex loc, int iters){	// the main iteration loop
	int ct = 0;
	double za = 0.0;
	double zb = 0.0;
	double tmp = 0.0;
	double zasq = 0.0;
	double zbsq = 0.0;
	while (zasq + zbsq < bailout_sqr && ct < iters){
		tmp = zasq - zbsq + loc.a;
		zb = 2*za*zb + loc.b; //zb*za + za*zb;
		za = tmp;
		zasq = za*za;
		zbsq = zb*zb;
		
		ct++;
	}
	/*
	Complex z = {0.0,0.0};
	Complex temp;
	while (dist_sqr(z) < 1e80 && ct < iters){
		z = mul(z, z);
		z = add(z,loc);
		ct++;
	}
	*/
	ICT += ct;
	if (ct == iters) return -1;
	return ct - log(log(sqrt(za*za+zb*zb)/log_b))/log_2;	// for smooth coloring
}

int map (double v, mand_opts *opts){	// maps onto 0..512 for access into the color ramp.
	if (v == -1) return -1;
	v = pow(v, opts->exponent);
	v = v*opts->linear + opts->constant;
	v *= RPOINTS;
	return ((int) (v))%RPOINTS;	// for now
}

const char *version = "1.0";
const bool print_comp_time = true;

void print_usage (){
	printf("Usage: \n \
		Argument 1: # of threads to use \n \
		Argument 2: # of columns on terminal window \n \
		Argument 3: name of options file. \n");
}

int main (int argc, const char *argv[]){
	pthread_mutex_init(&print_mutex, NULL);
	pthread_mutex_init(&tile_manager, NULL);
	init_cramp();

	printf("Fast Mandelbrot Set Renderer - version %s - compiled %s %s\n", version, __DATE__, print_comp_time?__TIME__:"");

	if (argc < 3 || !strcmp(argv[1], "-h")){
		print_usage();
		exit(0);
	}

	int nthreads = atoi(argv[1]);
	printf("Using %i threads\n", nthreads);

	ncols = atoi(argv[2]);

	mand_opts opts;	// opts is the same across all threads.

	/* initialize defaults */
	opts.xsize = 700;
	opts.ysize = 700;
	opts.tilesize = 64;
	Complex corner = {-2.0, -1.5};
	Complex size   = {3.0, 3.0};
	opts.ll = corner;
	opts.size = size;
	opts.iter = 500;
	opts.exponent = 0.6;
	opts.linear = 0.11;
	opts.constant = 0;
	opts.refinement_passes = 4;
	opts.qcx_thresh = 0;
	opts.grid_thresh = 0.2;
	opts.grid_size = 1.3;
	opts.grid_fac = 4;
	opts.grid_wt_tot = 3;
	opts.qcx_wt = 0.5;

	/* Here we read in the options from a file. Anything specified here will override the defaults. */ 
	FILE *inp = fopen(argv[3], "r");
	char name[20];
	while (!feof(inp)){
		fscanf(inp, "%s", name);
		if (feof(inp)){
			break;
		}
		printf("Found name %s\n", name);
		if (strcmp(name, "xsize") == 0){
			fscanf(inp, "%i", &opts.xsize);
		} else if (strcmp(name, "ysize") == 0){
			fscanf(inp, "%i", &opts.ysize);
		} else if (strcmp(name, "corner") == 0){
			fscanf(inp, "%lf %lf", &opts.ll.a, &opts.ll.b);
		} else if (strcmp(name, "size") == 0){
			fscanf(inp, "%lf %lf", &opts.size.a, &opts.size.b);
		} else if (strcmp(name, "iter") == 0){
			fscanf(inp, "%i", &opts.iter);
		} else if (strcmp(name, "exp") == 0){
			fscanf(inp, "%lf", &opts.exponent);
		} else if (strcmp(name, "linear") == 0){
			fscanf(inp, "%lf", &opts.linear);
		} else if (strcmp(name, "offset") == 0){
			fscanf(inp, "%lf", &opts.constant);
		} else if (strcmp(name, "ref_passes") == 0){
			fscanf(inp, "%i", &opts.refinement_passes);
		} else if (strcmp(name, "qcx_thresh") == 0){
			fscanf(inp, "%lf", &opts.qcx_thresh);
		} else if (strcmp(name, "grid_thresh") == 0){
			fscanf(inp, "%lf", &opts.grid_thresh);
		} else if (strcmp(name, "grid_size") == 0){
			fscanf(inp, "%lf", &opts.grid_size);
		} else if (strcmp(name, "grid_fac") == 0){
			fscanf(inp, "%i", &opts.grid_fac);
		} else if (strcmp(name, "grid_wt_tot") == 0){
			fscanf(inp, "%lf", &opts.grid_wt_tot);
		} else if (strcmp(name, "qcx_wt") == 0){
			fscanf(inp, "%lf", &opts.qcx_wt);
		} else if (strcmp(name, "tilesize") == 0){
			fscanf(inp, "%i", &opts.tilesize);
		} else {
			printf("Unrecognized name: %s; aborting.\n", name);
			exit(1);
		}
	}
	fclose(inp);

	compute_kernel(opts.grid_fac);

	int ntiles = (opts.xsize%opts.tilesize==0?0:1)+(opts.xsize/opts.tilesize);
	ntiles *= (opts.ysize%opts.tilesize==0?0:1)+(opts.ysize/opts.tilesize);
	printf("Will render a total of %i tiles (about %i per thread)\n", ntiles, ntiles/nthreads);

	color *img = (color *)(malloc(opts.xsize*opts.ysize*sizeof(color)));
	if (img == NULL){
		printf("bad image - could not malloc\n");
		exit(1);
	}

	pthread_t *threads = (pthread_t *)(malloc(nthreads*sizeof(pthread_t)));
	render_thread_options *thread_opts = (render_thread_options *)(malloc(nthreads*sizeof(render_thread_options)));
	messages = (char **)(malloc(nthreads*sizeof(char*)));
	long time = clock();
	for (int i=0; i<nthreads; i++){
		thread_opts[i].nthreads = nthreads;
		thread_opts[i].thr_num = i;
		thread_opts[i].mopts = &opts;
		thread_opts[i].img = img;
		messages[i] = (char *)(malloc(128*sizeof(char)));	// messages should be shorter than 128 chars.

		pthread_create(&threads[i], NULL, render, (void *) &thread_opts[i]);
		
	}

	for (int i=0; i<nthreads; i++){
		pthread_join(threads[i], NULL);
	}


	printf("Render threads have stopped; user time = %i ms. ", (int) ((clock()-time)/1000));
	time = clock();

	FILE *out = fopen("test.bmp", "w");
	printf("Proceeding with image writing...\n");
	write_img(img, opts.xsize, opts.ysize, out);

	printf("Image writing took %i ms.\n", (int) ((clock()-time)/1000));
	printf("ICT = %lu\n", ICT);
	printf("Quincunx supersampling: %i samples.\n", QCX_TOT_PX);
	printf("Full grid supersampling: %i pixels.\n", GRID_TOT_PX);
	return 0;
}
