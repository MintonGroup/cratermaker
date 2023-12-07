#define STRMAX 512

typedef struct c_surface_type {
    double *elevation;
    int elevation_len;
    double *node_x;
    int node_x_len;
    double *node_y;
    int node_y_len;
    double *node_z;
    int node_z_len;
    double *face_x;
    int face_x_len;
    double *face_y;
    int face_y_len;
    double *face_z;
    int face_z_len;
}surface_type;
extern struct c_surface_type* bind_surface_init(int ny, int nx);
extern void bind_surface_final(struct c_surface_type *obj);
extern void bind_surface_set_name(struct c_surface_type *obj, const char *c_string);
extern char* bind_surface_get_name(struct c_surface_type *obj);

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

double bind_perlin_noise(const char *model, double x, double y, double z, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0);