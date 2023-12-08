#define STRMAX 512

typedef struct c_surface_type {
    double *elevation;
    int elevation_shape[1];
    double *node_x;
    int node_x_shape[1];
    double *node_y;
    int node_y_shape[1];
    double *node_z;
    int node_z_shape[1];
    double *face_x;
    int face_x_shape[1];
    double *face_y;
    int face_y_shape[1];
    double *face_z;
    int face_z_shape[1];
}surface_type;
extern struct c_surface_type* bind_surface_init(int ny, int nx);
extern void bind_surface_final(struct c_surface_type *obj);

typedef struct c_PerlinArguments {
      double damp;
      double damp0;
      double damp_scale;
      double freq;
      double gain;
      double gain0;
      double lacunarity;
      double noise_height;
      double pers;
      double slope;
      double warp;
      double warp0;
}PerlinArguments;

double bind_perlin_noise(const char *model, double *x, double *y, double *z, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0);