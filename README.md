# ninety3d
The missing 3d renderer from the 90s

## Highlights

- scanline rendering
- edge drawing
- arbitrary lightsource location
- ambient lighting
- objects must be normalized
- no math library required
- no stdlib required
- zero dependencies
- define N3D_IMPLEMENTATION before inclusion

## Configuration

name              | default | type  | note
------------------|---------|-------|--------------
N3D_CAM_X         | 1.5     | float | camera location on X axis
N3D_OBJECT_RED    | 100     | int   | object color
N3D_OBJECT_GREEN  | 100     | int   |
N3D_OBJECT_BLUE   | 150     | int   |
N3D_LIGHT_X       | 10.5    | float | lightsource location
N3D_LIGHT_Y       | 11.0    | float |
N3D_LIGHT_Z       | 11.0    | float |
N3D_AMBIENT_LIGHT | 32      | float | ambient light
N3D_SW            | 240     | int   | window width
N3D_SH            | 240     | int   | window height
N3D_LINE_R        | 255     | int   | edge line color
N3D_LINE_G        | 255     | int   |
N3D_LINE_B        | 255     | int   |

## Usage

1. Keep the object as triangle mesh in the following struct:

   ```c
   typedef struct
   {
     float x[3], y[3], z[3];
     float depth;            // keep empty for internal use
   } tri;
   ```

   Reading an STL file directly to this structure will work.

2. Implement the optional callback interface functions

   ```c
   void draw_line(int r, int g, int b, int x0, int y0, int x1, int y1)
   {
   }
   
   void draw_horizontal_line(int r, int g, int b, int y, int x0, int x1)
   {
   }
   ```

4. Initialize the n3d struct

   ```c
   typedef struct
   {
     tri *triangles;
     int triangle_count;
     draw_line_type *draw_line;
     draw_horizontal_line_type *draw_horizontal_line;
   } n3d;

   n3d ninety3d;
   
   ninety3d.triangles=triangles;
   ninety3d.triangle_count=number_of_triangles;
   ninety3d.draw_line=NULL;
   ninety3d.draw_horizontal_line=draw_horizontal_line;
   ```

5. Call the interation function from your main loop
   ```c
   n3d_scene( &ninety3d, angle_x, angle_y, angle_z, N3D_NOFLIP);
   ```
## Demo

Pinetime ST7789v demo

https://github.com/user-attachments/assets/e198a00d-b2b0-4b1a-beb6-43dda7f9212d
