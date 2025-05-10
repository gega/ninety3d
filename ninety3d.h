/*
 *  nighty3d.h
 *
 *  3d software renderer
 *
 *  - scanline rendering
 *  - edge drawing
 *  - arbitrary lightsource location
 *  - ambient light
 *  - objects must be normalized and transformed to origo
 *  - no math library required
 *  - no stdlib required
 *  - zero dependencies
 *  - define N3D_IMPLEMENTATION before inclusion
 *
 */
#ifndef N3D_H
#define N3D_H

#include <stdint.h>

/* config area, override these values here or outside
 */

#ifndef N3D_CAM_X
#define N3D_CAM_X (1.5f)
#endif
#ifndef N3D_OBJECT_RED
#define N3D_OBJECT_RED (100)
#endif
#ifndef N3D_OBJECT_GREEN
#define N3D_OBJECT_GREEN (100)
#endif
#ifndef N3D_OBJECT_BLUE
#define N3D_OBJECT_BLUE (150)
#endif
#ifndef N3D_LIGHT_X
#define N3D_LIGHT_X (10.5f)
#endif
#ifndef N3D_LIGHT_Y
#define N3D_LIGHT_Y (11.0f)
#endif
#ifndef N3D_LIGHT_Z
#define N3D_LIGHT_Z (11.0f)
#endif
#ifndef N3D_AMBIENT_LIGHT
#define N3D_AMBIENT_LIGHT (32)
#endif
#ifndef N3D_SW
#define N3D_SW (240) /* screen width */
#endif
#ifndef N3D_SH
#define N3D_SH (240) /* screen height */
#endif
#ifndef N3D_LINE_R
#define N3D_LINE_R (255)
#endif
#ifndef N3D_LINE_G
#define N3D_LINE_G (255)
#endif
#ifndef N3D_LINE_B
#define N3D_LINE_B (255)
#endif

/* end of config area
 */

#define N3D_NOFLIP (-1)


typedef struct
{
    float x[3], y[3], z[3];
    float depth;
} tri;

typedef void draw_line_type(int,int,int,int,int,int,int);
typedef void draw_horizontal_line_type(int,int,int,int,int,int);

typedef struct
{
  tri *triangles;
  int triangle_count;
  draw_line_type *draw_line;
  draw_horizontal_line_type *draw_horizontal_line;
} n3d;


#ifdef N3D_IMPLEMENTATION


#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) ((sizeof x) / (sizeof *x))
#endif

#ifndef NULL
#define NULL ((void *)0L)
#endif


static float n3d_rot_x[3][3], n3d_rot_y[3][3], n3d_rot_z[3][3];

/*
 * void gen_trig_tables(int res)
 * {
 *   for(int i=0; i<res; i++)
 *   {
 *     float angle = ((float)i/res) * 2*M_PI;
 *     sin_table[i] = sinf(angle);
 *     cos_table[i] = cosf(angle);
 *   }
 * }
 */
static const float n3d_sin_table[]=
{
   0.000, 0.025, 0.049, 0.074, 0.098, 0.122, 0.147, 0.171, 0.195, 0.219, 0.243, 0.267, 0.290, 0.314, 0.337, 0.360,
   0.383, 0.405, 0.428, 0.450, 0.471, 0.493, 0.514, 0.535, 0.556, 0.576, 0.596, 0.615, 0.634, 0.653, 0.672, 0.690,
   0.707, 0.724, 0.741, 0.757, 0.773, 0.788, 0.803, 0.818, 0.831, 0.845, 0.858, 0.870, 0.882, 0.893, 0.904, 0.914,
   0.924, 0.933, 0.942, 0.950, 0.957, 0.964, 0.970, 0.976, 0.981, 0.985, 0.989, 0.992, 0.995, 0.997, 0.999, 1.000,
   1.000, 1.000, 0.999, 0.997, 0.995, 0.992, 0.989, 0.985, 0.981, 0.976, 0.970, 0.964, 0.957, 0.950, 0.942, 0.933,
   0.924, 0.914, 0.904, 0.893, 0.882, 0.870, 0.858, 0.845, 0.831, 0.818, 0.803, 0.788, 0.773, 0.757, 0.741, 0.724,
   0.707, 0.690, 0.672, 0.653, 0.634, 0.615, 0.596, 0.576, 0.556, 0.535, 0.514, 0.493, 0.471, 0.450, 0.428, 0.405,
   0.383, 0.360, 0.337, 0.314, 0.290, 0.267, 0.243, 0.219, 0.195, 0.171, 0.147, 0.122, 0.098, 0.074, 0.049, 0.025,
  -0.000,-0.025,-0.049,-0.074,-0.098,-0.122,-0.147,-0.171,-0.195,-0.219,-0.243,-0.267,-0.290,-0.314,-0.337,-0.360,
  -0.383,-0.405,-0.428,-0.450,-0.471,-0.493,-0.514,-0.535,-0.556,-0.576,-0.596,-0.615,-0.634,-0.653,-0.672,-0.690,
  -0.707,-0.724,-0.741,-0.757,-0.773,-0.788,-0.803,-0.818,-0.831,-0.845,-0.858,-0.870,-0.882,-0.893,-0.904,-0.914,
  -0.924,-0.933,-0.942,-0.950,-0.957,-0.964,-0.970,-0.976,-0.981,-0.985,-0.989,-0.992,-0.995,-0.997,-0.999,-1.000,
  -1.000,-1.000,-0.999,-0.997,-0.995,-0.992,-0.989,-0.985,-0.981,-0.976,-0.970,-0.964,-0.957,-0.950,-0.942,-0.933,
  -0.924,-0.914,-0.904,-0.893,-0.882,-0.870,-0.858,-0.845,-0.831,-0.818,-0.803,-0.788,-0.773,-0.757,-0.741,-0.724,
  -0.707,-0.690,-0.672,-0.653,-0.634,-0.615,-0.596,-0.576,-0.556,-0.535,-0.514,-0.493,-0.471,-0.450,-0.428,-0.405,
  -0.383,-0.360,-0.337,-0.314,-0.290,-0.267,-0.243,-0.219,-0.195,-0.171,-0.147,-0.122,-0.098,-0.074,-0.049,-0.025,
};
static const float n3d_cos_table[]=
{
   1.000, 1.000, 0.999, 0.997, 0.995, 0.992, 0.989, 0.985, 0.981, 0.976, 0.970, 0.964, 0.957, 0.950, 0.942, 0.933,
   0.924, 0.914, 0.904, 0.893, 0.882, 0.870, 0.858, 0.845, 0.831, 0.818, 0.803, 0.788, 0.773, 0.757, 0.741, 0.724,
   0.707, 0.690, 0.672, 0.653, 0.634, 0.615, 0.596, 0.576, 0.556, 0.535, 0.514, 0.493, 0.471, 0.450, 0.428, 0.405,
   0.383, 0.360, 0.337, 0.314, 0.290, 0.267, 0.243, 0.219, 0.195, 0.171, 0.147, 0.122, 0.098, 0.074, 0.049, 0.025,
  -0.000,-0.025,-0.049,-0.074,-0.098,-0.122,-0.147,-0.171,-0.195,-0.219,-0.243,-0.267,-0.290,-0.314,-0.337,-0.360,
  -0.383,-0.405,-0.428,-0.450,-0.471,-0.493,-0.514,-0.535,-0.556,-0.576,-0.596,-0.615,-0.634,-0.653,-0.672,-0.690,
  -0.707,-0.724,-0.741,-0.757,-0.773,-0.788,-0.803,-0.818,-0.831,-0.845,-0.858,-0.870,-0.882,-0.893,-0.904,-0.914,
  -0.924,-0.933,-0.942,-0.950,-0.957,-0.964,-0.970,-0.976,-0.981,-0.985,-0.989,-0.992,-0.995,-0.997,-0.999,-1.000,
  -1.000,-1.000,-0.999,-0.997,-0.995,-0.992,-0.989,-0.985,-0.981,-0.976,-0.970,-0.964,-0.957,-0.950,-0.942,-0.933,
  -0.924,-0.914,-0.904,-0.893,-0.882,-0.870,-0.858,-0.845,-0.831,-0.818,-0.803,-0.788,-0.773,-0.757,-0.741,-0.724,
  -0.707,-0.690,-0.672,-0.653,-0.634,-0.615,-0.596,-0.576,-0.556,-0.535,-0.514,-0.493,-0.471,-0.450,-0.428,-0.405,
  -0.383,-0.360,-0.337,-0.314,-0.290,-0.267,-0.243,-0.219,-0.195,-0.171,-0.147,-0.122,-0.098,-0.074,-0.049,-0.025,
   0.000, 0.025, 0.049, 0.074, 0.098, 0.122, 0.147, 0.171, 0.195, 0.219, 0.243, 0.267, 0.290, 0.314, 0.337, 0.360,
   0.383, 0.405, 0.428, 0.450, 0.471, 0.493, 0.514, 0.535, 0.556, 0.576, 0.596, 0.615, 0.634, 0.653, 0.672, 0.690,
   0.707, 0.724, 0.741, 0.757, 0.773, 0.788, 0.803, 0.818, 0.831, 0.845, 0.858, 0.870, 0.882, 0.893, 0.904, 0.914,
   0.924, 0.933, 0.942, 0.950, 0.957, 0.964, 0.970, 0.976, 0.981, 0.985, 0.989, 0.992, 0.995, 0.997, 0.999, 1.000,
};
#define n3d_sin_fast(angle) n3d_sin_table[(int)((angle) * (ARRAY_SIZE(n3d_sin_table) / 6.2831853f)) % ARRAY_SIZE(n3d_sin_table)]
#define n3d_cos_fast(angle) n3d_cos_table[(int)((angle) * (ARRAY_SIZE(n3d_cos_table) / 6.2831853f)) % ARRAY_SIZE(n3d_cos_table)]
#define n3d_tan_fast(angle) (n3d_sin_fast((angle))/n3d_cos_fast((angle)))


static inline void n3d_mini_ssort_swap(void * const a_, void * const b_, int size)
{
  char * const a = a_;
  char * const b = b_;
  char tmp;

  for(int i=0; i<size; i++)
  {
    tmp = a[i];
    a[i] = b[i];
    b[i] = tmp;
  }
}

/*
 * Marcin Ciura / 2001 / Best Increments for the Average Case of Shellsort
 * https://link.springer.com/chapter/10.1007/3-540-44669-9_12
 */
static void n3d_mini_ssort(void *base, int nmemb, int size, int (*cmp)(const void *, const void *))
{
  const static unsigned short gaps[] = { 44861, 19938, 8861, 3938, 1750, 701, 301, 132, 57, 23, 10, 4, 1 };
  typedef char type[size];
  type *array = base;

  for(const unsigned short *gap = gaps; gap < gaps + ARRAY_SIZE(gaps); gap++)
  {
    for(int i = *gap; i < nmemb; i++)
    {
      for(int j = i; j >= *gap && cmp(array[j - *gap], array[j]) > 0; j -= *gap)
      {
        n3d_mini_ssort_swap(array[j], array[j - *gap], size);
      }
    }
  }
}

static void n3d_compute_rotation(float ax, float ay, float az)
{
  float cx = n3d_cos_fast(ax), sx = n3d_sin_fast(ax);
  float cy = n3d_cos_fast(ay), sy = n3d_sin_fast(ay);
  float cz = n3d_cos_fast(az), sz = n3d_sin_fast(az);

  n3d_rot_x[0][0] = 1;   n3d_rot_x[0][1] = 0;   n3d_rot_x[0][2] = 0;
  n3d_rot_x[1][0] = 0;   n3d_rot_x[1][1] = cx;  n3d_rot_x[1][2] = -sx;
  n3d_rot_x[2][0] = 0;   n3d_rot_x[2][1] = sx;  n3d_rot_x[2][2] = cx;

  n3d_rot_y[0][0] = cy;  n3d_rot_y[0][1] = 0;   n3d_rot_y[0][2] = sy;
  n3d_rot_y[1][0] = 0;   n3d_rot_y[1][1] = 1;   n3d_rot_y[1][2] = 0;
  n3d_rot_y[2][0] = -sy; n3d_rot_y[2][1] = 0;   n3d_rot_y[2][2] = cy;

  n3d_rot_z[0][0] = cz;  n3d_rot_z[0][1] = -sz; n3d_rot_z[0][2] = 0;
  n3d_rot_z[1][0] = sz;  n3d_rot_z[1][1] = cz;  n3d_rot_z[1][2] = 0;
  n3d_rot_z[2][0] = 0;   n3d_rot_z[2][1] = 0;   n3d_rot_z[2][2] = 1;
}

static void n3d_rotate(float in[3], float out[3])
{
  float tmp1[3], tmp2[3];
  int i;

  for(i = 0; i < 3; i++) tmp1[i] = n3d_rot_x[i][0]*in[0]   + n3d_rot_x[i][1]*in[1]   + n3d_rot_x[i][2]*in[2];
  for(i = 0; i < 3; i++) tmp2[i] = n3d_rot_y[i][0]*tmp1[0] + n3d_rot_y[i][1]*tmp1[1] + n3d_rot_y[i][2]*tmp1[2];
  for(i = 0; i < 3; i++) out[i]  = n3d_rot_z[i][0]*tmp2[0] + n3d_rot_z[i][1]*tmp2[1] + n3d_rot_z[i][2]*tmp2[2];
}

static void n3d_project(float in[3], int *x, int *y, float *depth)
{
  static const float fov_deg = 70.0f;
  static const float fov_rad = fov_deg * (3.1415926f / 180.0f);
  static const float scale = 0.5f * N3D_SW / n3d_tan_fast(fov_rad * 0.5f);

  float px = scale * in[0] / (N3D_CAM_X - in[2]);
  float py = scale * in[1] / (N3D_CAM_X - in[2]);
  *x = (int)(px + N3D_SW / 2);
  *y = (int)(-py + N3D_SH / 2);
  *depth = in[2];
}

/*
 * https://en.wikipedia.org/wiki/Fast_inverse_square_root
 */
static float n3d_fast_rsqrt(float number)
{
  union { float f; uint32_t i; } conv = { number };

  conv.i = 0x5f3759df - (conv.i >> 1);
  float y = conv.f;
  y = y * (1.5f - 0.5f * number * y * y);

  return(y);
}

static int n3d_cmp_depth(const void *a, const void *b)
{
  float d1 = ((tri*)a)->depth;
  float d2 = ((tri*)b)->depth;
  return (d1 > d2) - (d1 < d2);
}

void n3d_scene(n3d *n3, float angle_x, float angle_y, float angle_z, int flip)
{
  int i,v;

  n3d_compute_rotation(angle_x, angle_y, angle_z);

  float transformed_tri[n3->triangle_count][3][3];
  for(i=0; i<n3->triangle_count; i++)
  {
    float avg_z = 0;
    for (int v = 0; v < 3; v++)
    {
      float p[3] = { n3->triangles[i].x[v], n3->triangles[i].y[v], n3->triangles[i].z[v] };
      n3d_rotate(p, transformed_tri[i][v]);
      avg_z += transformed_tri[i][v][2];
    }
    n3->triangles[i].depth=avg_z / 3.0f;
  }

  n3d_mini_ssort(n3->triangles, n3->triangle_count, sizeof(tri), n3d_cmp_depth);

  for(i = 0; i < n3->triangle_count; i++)
  {
    float (*transformed)[3] = transformed_tri[i];

    float ux = transformed[1][0] - transformed[0][0];
    float uy = transformed[1][1] - transformed[0][1];
    float uz = transformed[1][2] - transformed[0][2];
    float vx = transformed[2][0] - transformed[0][0];
    float vy = transformed[2][1] - transformed[0][1];
    float vz = transformed[2][2] - transformed[0][2];

    float nx = uy * vz - uz * vy;
    float ny = uz * vx - ux * vz;
    float nz = ux * vy - uy * vx;

    float view_dir[3] = { 0.0f     - transformed[0][0],
                          0.0f     - transformed[0][1],
                         N3D_CAM_X - transformed[0][2] };

    float view_dot = nx * view_dir[0] + ny * view_dir[1] + nz * view_dir[2];
    if (view_dot < 0) continue;

    float inv_len = n3d_fast_rsqrt(nx * nx + ny * ny + nz * nz); // float nl = sqrtf(nx * nx + ny * ny + nz * nz);
    nx *= inv_len;						 // nx /= nl;
    ny *= inv_len;						 // ny /= nl;
    nz *= inv_len;						 // nz /= nl;

    float lx = N3D_LIGHT_X, ly = N3D_LIGHT_Y, lz = N3D_LIGHT_Z;
    float inv_ll = n3d_fast_rsqrt(lx*lx + ly*ly + lz*lz);       // float ll = sqrtf(lx*lx + ly*ly + lz*lz);
    lx *= inv_ll;						// lx /= ll;
    ly *= inv_ll;						// ly /= ll;
    lz *= inv_ll;						// lz /= ll;

    float dot = nx * lx + ny * ly + nz * lz;
    if (dot < 0) dot = 0;

    int r = (N3D_AMBIENT_LIGHT + dot * (N3D_OBJECT_RED   - N3D_AMBIENT_LIGHT));
    int g = (N3D_AMBIENT_LIGHT + dot * (N3D_OBJECT_GREEN - N3D_AMBIENT_LIGHT));
    int b = (N3D_AMBIENT_LIGHT + dot * (N3D_OBJECT_BLUE  - N3D_AMBIENT_LIGHT));

    int screen[3][2];
    float z[3];
    for(v = 0; v < 3; v++) n3d_project(transformed[v], &screen[v][0], &screen[v][1], &z[v]);

    int miny = N3D_SH, maxy = 0;
    for(v = 0; v < 3; v++)
    {
        if (screen[v][1] < miny) miny = screen[v][1];
        if (screen[v][1] > maxy) maxy = screen[v][1];
    }

    if(miny < 0) miny = 0;
    if(maxy >= N3D_SH) maxy = N3D_SH-1;

    for(int y = miny; y <= maxy; y++)
    {
      if((y%2)==flip) continue;
      float xint[3];
      int n=0;
      for(int i0=0; i0<3 && n<2; i0++)
      {
        int i1 = (i0 + 1) % 3;
        int y0 = screen[i0][1], y1 = screen[i1][1];
        int x0 = screen[i0][0], x1 = screen[i1][0];
        if ((y >= y0 && y < y1) || (y >= y1 && y < y0))
        {
          float f = (float)(y - y0) / (float)(y1 - y0);
          xint[n++] = x0 + f * (x1 - x0);
        }
      }
      if(n>=2)
      {
        int x0 = (int)xint[0];
        int x1 = (int)xint[1];
        if (x0 > x1) { int tmp = x0; x0 = x1; x1 = tmp; }
        if (x0 < 0) x0 = 0;
        if (x1 >= N3D_SW) x1 = N3D_SW - 1;
        if(NULL!=n3->draw_horizontal_line) n3->draw_horizontal_line(r,g,b, y, x0, x1);
      }
    }

    for(v = 0; v < 3; v++)
    {
      int v1 = (v + 1) % 3;
      if(NULL!=n3->draw_line) n3->draw_line(N3D_LINE_R,N3D_LINE_G,N3D_LINE_B, screen[v][0], screen[v][1], screen[v1][0], screen[v1][1]);
    }
  }
}

#endif
#endif
